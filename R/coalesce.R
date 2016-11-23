#Christopher Barry, started on 16/03/2015 at University of Birmingham

# functions for coalescing a set of particles, clumping groups that are close together

if(!require("plyr")){install.packages("plyr"); require("plyr")}
if(!require("data.table")){install.packages("data.table"); require("data.table")}

# a random field of points, for testing
tstpts <- as.data.table(matrix(c(sample(0:50000, 100, T)/500, sample(0:50000, 100, T)/500,
                   sample((1:100)/100, 100, T)),
                 ncol = 3, dimnames = list(NULL, c("x", "y", "m"))))
tstpts3D <- cbind(tstpts[, 1:2, with = F], z = sample((0:100)/10, 100, T), m = tstpts[, 3, with = F])
big.tst <- function(n){
  tstpts <<- as.data.table(matrix(c(sample(0:50000, n, T)/500, sample(0:50000, n, T)/500,
                                    sample((1:100)/100, n, T)),
                                  ncol = 3, dimnames = list(NULL, c("x", "y", "m"))))
  tstpts3D <<- cbind(tstpts[, 1:2, with = F], z = sample((0:100)/10, n, T), m = tstpts[, m])
}

# xyzm should be a matrix or table of four (or three if TwoD = TRUE) columns: x, y, z, mass in that order
# cdh and cdv are the coalescing radii in the horizontal and vertical directions (forming an ellipsoid)
# mm is the minimum mass for resulting particles: any particles below this mass are coalesced with the nearest particle regardless of distance
# subregions = TRUE is recommended, even though some potential joins may be missed, because it greatly speeds up the code for large collections of particles: it analyses the inter-particle distances in subregions rather than as the whole set
# TwoD = TRUE tells the function that the particles are distributed in two dimensions only.  xyzm should be interpreted as xym in this case.
# nppsr is the maximum number of particles per subregion: the default 256 is the most efficient for large numbers of particles
# na.rep specifies whether rows containing NA in either the spatial or mass columns (which are removed during the function's calculations) should be reaffixed after the calculation, so they are not lost
# The final two arguments allow the user to set a desired maximum number of particles for the result: maxnp is that maximum (limitless by default), maxbymm = TRUE indicates that the threshold should be achieved by increasing mm, otherwise cd would be increased if necessary.  (note: 1. maxbymm = FALSE is not yet developed - will have no effect currently; 2. if mm = 0, then this will be adjusted to min(xyzm$m) so no fear of infinite loop)
coalesce3Dm <- function(xyzm, cdh, cdv = cdh/10, mm = 0, subregions = TRUE, TwoD = FALSE, nppsr = 256L, na.rep = FALSE, maxnp = Inf, maxbymm = TRUE){
  st.time <- Sys.time()

  cns <- copy(colnames(xyzm)) # copy makes a deep copy, removing the by reference relationship to the xyzm data table, which would cause it to retain a link with xyzm and therefore be modified along with it
  xyzm <- as.data.table(xyzm)
  
  # Remove rows containing NA for the calculation. Optionally save them for reattachment at the end
  nas <- Reduce(`|`, lapply(xyzm, is.na))
  if(na.rep && any(nas)) xyzmna <- xyzm[nas,]
  if(any(nas)) xyzm <- xyzm[!nas,]
  
  if(!is.finite(cdv)) cdv <- 0
  if(cdh == 0) cdv <- 1 #because otherwise get 0/0: doesn't change result in this case
  nan.flag <- cdv == 0

  setnames(xyzm, c("x", "y", if(!TwoD) "z", "m"))

  xyzm <- xyzm[m != 0,] #remove zero-mass particles
  np <- nrow(xyzm)

  # number of subregions: the largest square number that will allow at least 256 (or nppsr) points per subregion
  if(subregions){
    sqrtnsr <- as.integer(sqrt(np)%/%sqrt(nppsr))
    nsr <- as.integer(sqrtnsr^2)
    if(nsr < 4L) subregions <- F #not enough particles to justify subregions (or effectively one subregion)
  }
  
  #groups points into 16 subgroups to reduce the number of distance calculations that must be performed (no point if small number of particles, and should not be done if cd is large compared to bounding box)
  #grouping is by x and y only
  if(subregions){
    # divide into roughly equally populated sub-y-sections
    setkey(xyzm, y)
    xyzm[, sry := as.integer(.I%/%ceiling((np + 1L)/sqrtnsr))]
    
    #divide the sub-y-sections into roughly equally populated sub-x-sections
    setkey(xyzm, sry, x)
    xyzm[, srx := as.integer(.I%/%ceiling((.N + 1L)/sqrtnsr)), by = sry]
    
    #give each subregion a single unique identifier
    xyzm[, sr := .GRP, by = c("sry", "srx")]
    xyzm[, c("srx", "sry") := NULL]
  }
  
  sr.time <- Sys.time()
  
  if(cdh > 0){
    setkey(xyzm, m)
    
    grlabs <- integer(np) #group to which assigned
    
    xyzm[, group := {
      grvec <- integer(.N)
      asd <- logical(.N) # vector of FALSEs
      
      #those particle numbers within particle p's subregion to which p is close, including self
      cl <- l_ply(.N:1, function(p){
        #if this particle is already assigned to a group, then skip
        if(asd[p]) return(integer(0L))
        
        #find all particles (within this subregion) that are within the coalescing radius of this particle
        wcl <- which({
          zsepsq <- if(TwoD) 0 else ((z[p] - z)^2)*(cdh/cdv)^2
          if(nan.flag && !TwoD) zsepsq[is.nan(zsepsq)] <- 0 # 0*Inf = NaN, but would like 0 here (when cdv = 0 and equal z values are compared)
          (x[p] - x)^2 + (y[p] - y)^2 + zsepsq <= cdh^2
        })
        
        #mark the close-by particles as assigned
        asd[wcl] <<- TRUE
        
        #assign group labels - doesn't matter what the labels are as long as they are distinct
        #only need to be distinct within subregion because by = c("group", "sr") is used for the grouping
        grvec[wcl] <<- p
      })
      
      #return value
      grvec
    }, by = if(subregions) sr]
  }
  
  cl.time <- Sys.time()

  # aggregate with data table by group label
  # don't bother if no coalescing actually requested
  # `[.data.table` parses its arguments unusually: "by" must be specified
  #  outside the statement if not simple
  if(cdh > 0){
    bys <- if(subregions) "group,sr" else "group"
    xyzm <- if(TwoD){
      xyzm[, if(identical(.N, 1L)){
        list(x = x, y = y, m = m)
      }else{
        xy <- llply(list(x = x, y = y), weighted.mean, m)
        c(xy, list(m = sum(m)))
      }, by = bys]
    }else{
      xyzm[, if(identical(.N, 1L)){
        list(x = x, y = y, z = z, m = m)
      }else{
        xyz <- llply(list(x = x, y = y, z = z), weighted.mean, m)
        c(xyz, list(m = sum(m)))
      }, by = bys]
    }
  }
  
  gr.time <- Sys.time()
  
  np <- nrow(xyzm)
  
  #which particles are smaller than mm, or which are the smallest particles if there is a set maximum count
  if(maxbymm && mm <= 0) mm <- min(xyzm$m)
  
  #when maxnp != Inf or mm > 0 and subregions = TRUE, there is the possibility of subregions with one very small particle
  #in this case, an infinite loop occurs because there is no particle for these small particles to merge with
  #in order to avoid this, a variable called mnp.modifier is used that will modify the used value for maxnp
  #it is updated with each loop to be the number of isolated small particles which are the sole occupants of their subregions
  #the result is that maxnp is effectively reduced when finding which particles to flag as small, so that other subregions have reduced particle counts
  #mm > 0 is dealt with in a different way: groups with one small particle are given a special -1 group label which are then ignored
  #thus is is possible for some particles with m < mm to be retained
  if(subregions && is.finite(maxnp)) mnp.modifier <- 0L #this will be reset with each loop
  
  smallf <- expression({
    if(maxbymm && np > maxnp){
      pflags <- logical(np)
      msrt <- sort.list(xyzm$m, method = "quick", na.last = NA)
      # the lightest np - maxnp particles are labelled as small
      pflags[msrt[1:(np - (maxnp - mnp.modifier))]] <- TRUE; pflags
    }else{
      smtmp <- xyzm$m < mm
      smtmp[xyzm$group == -1L] <- FALSE #no point trying to group a subregion with only one particle - ignore
      smtmp
    }
  })
  
  #coalesce particles with mass less than mm with nearest particles until no particles have mass less than mm
  while(mm > 0 && any(sm <- eval(smallf))){
    if(identical(np, 1L)) break #infinite loop would be possible here
    if(all(sm)){
      #all particles are small - return one super-particle
      xyzm <- xyzm[, c(list(weighted.mean(x, m), weighted.mean(y, m)),
                       if(!TwoD) list(weighted.mean(z, m)), list(sum(m)))]
      break
    }
    if(subregions && is.finite(maxnp)) mnp.modifier <- 0L # reset
    
    xyzm[, small := sm]; rm(sm)
    
    #make groups to coalesce each small particle with one or more others
    # `[.data.table` parses its arguments unusually: "by" must be specified
    #  outside the statement if not simple
    bys <- if(subregions) "sr"
    xyzm[, group := {
      if(all(small)){
        #this step avoids infinite loops - see above
        if(subregions && is.finite(maxnp) && .N == 1L) mnp.modifier <<- mnp.modifier + 1L
        
        #all small: all to be coalesced together
        -1L
      }else{
        grvec <- 1:.N
        asd <- logical(.N)
        
        grvec[small] <- vapply(which(small), function(p){
          if(asd[p]) return(grvec[p])
          
          #find closest particle (within subregion) to this small particle
          #the closest particle's group label is assigned as this small particle's group label
          zsepsq <- if(TwoD) 0 else ((z[p] - z)^2)*(cdh/cdv)^2
          if(nan.flag && !TwoD) zsepsq[is.nan(zsepsq)] <- 0
          dsq <- (x[p] - x)^2 + (y[p] - y)^2 + zsepsq
          dsq[p] <- NA # don't want it to find itself; which.min ignores NAs
          closest <- which.min(dsq)
          
          asd[closest] <<- TRUE # mark as assigned to group already
          grvec[closest]
        }, integer(1L))
        
        grvec
      }
    }, by = bys]
    
    #coalesce groups due to smallness
    # `[.data.table` parses its arguments unusually: "by" must be specified
    #  outside the statement if not simple
    bys <- if(subregion) "group,sr" else "group"
    xyzm <- if(TwoD){
      xyzm[, if(identical(.N, 1L)){
        list(x = x, y = y, m = m)
      }else{
        xy <- llply(list(x = x, y = y), weighted.mean, m)
        c(xy, list(m = sum(m)))
      }, by = bys]
    }else{
      xyzm[, if(identical(.N, 1L)){
        list(x = x, y = y, z = z, m = m)
      }else{
        xyz <- llply(list(x = x, y = y, z = z), weighted.mean, m)
        c(xyz, list(m = sum(m)))
      }, by = bys]
    }

    #small particles are now moved to new homes, so to speak, so can delete their original entries
    #mass and centre of mass are conserved
    np <- nrow(xyzm)
  }

  mm.time <- Sys.time()

  timings <<- diff(c(st.time, subregions = sr.time, close = cl.time, coalescing = gr.time, minmass = mm.time))
  
  suppressWarnings(xyzm[, c("group", "sr", "small") := NULL])
  if(!is.null(cns)) setnames(xyzm, cns)
  if(na.rep && any(nas)) xyzm <- rbind(xyzm, xyzmna)
  
  return(xyzm)
}
