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
# nppsr is the maximum number of particles per subregion
# na.rep specifies whether rows containing NA in either the spatial or mass columns (which are removed during the function's calculations) should be reaffixed after the calculation, so they are not lost
# The final two arguments allow the user to set a desired maximum number of particles for the result: maxnp is that maximum (limitless by default), maxbymm = TRUE indicates that the threshold should be achieved by increasing mm, otherwise cd would be increased if necessary.  (note: 1. maxbymm = FALSE is not yet developed - will have no effect currently; 2. if mm = 0, then this will be adjusted to min(xyzm$m) so no fear of infinite loop)
coalesce3Dm <- function(xyzm, cdh, cdv = cdh/10, mm = 0, subregions = TRUE, TwoD = FALSE, nppsr = 256L, na.rep = FALSE, maxnp = Inf, maxbymm = TRUE){
  st.time <- Sys.time()

  xyzm <- as.data.table(xyzm)
  
  # Remove rows containing NA for the calculation. Optionally save them for reattachment at the end
  nas <- is.na(rowSums(xyzm, na.rm = FALSE))
  if(na.rep && any(nas)) xyzmna <- xyzm[nas,]
  if(any(nas)) xyzm <- xyzm[!nas,]
  
  if(!is.finite(cdv)) cdv <- 0
  if(cdh == 0) cdv <- 1 #because otherwise get 0/0: doesn't change result in this case
  nan.flag <- cdv == 0
  cns <- copy(colnames(xyzm)) # copy is necessary because of a quirk of setnames, which operates by reference

  setnames(xyzm, c(letters[23L + 1:ifelse(TwoD, 2L, 3L)], "m"))

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
    #defines four roughly equal groups by x; should be length 4 and not include np
    xsri <- seq(1L, np, as.integer(ceiling(np/sqrtnsr))) #start index of each x group
    dxsrinp <- diff(c(xsri, np + 1L)) #number of particles within each x group
    
    #defines sixteen groups: four y groups within the four x groups
    sri <- do.call(c, Map(function(start, n) seq(start, start + n - 1L, as.integer(ceiling(n/sqrtnsr))), xsri, dxsrinp)) #start index
    srend <- c(sri[-1], np + 1L) #end index + 1: used later on
    
    #re-order xym by x
    xyzm <- xyzm[sort.list(xyzm[, x], method = "quick", na.last = NA),]
    
    #re-order xym within each x group by y
    for(xsr in seq_len(sqrtnsr)){
      xyzm[xsri[xsr] + 0:(dxsrinp[xsr] - 1L),] <-
        xyzm[xsri[xsr] + 0:(dxsrinp[xsr] - 1L),][sort.list(xyzm[xsri[xsr] + 0:(dxsrinp[xsr] - 1L), y], method = "quick", na.last = NA),]
    }
  }else subregions <- FALSE
  
  sr.time <- Sys.time()
  
  if(cdh > 0){
    to <- sort.list(xyzm[, m], decreasing = TRUE, method = "quick", na.last = NA) #order in which to find close neighbours - order of decreasing mass, unstable ordering of ties (quick sort) does not matter
    if(subregions){
      srrgs <- Map(seq, sri, srend - 1L) #list of particle collections for each subregion
      srs <- vapply(seq_len(np), function(p){
        which(p >= sri & p < srend) #to which subregion does this point belong?
      }, integer(1L))
    }
    grlabs <- integer(np) #group to which assigned
    asd <- logical(np) #whether assigned (all FALSE at first)
    lab <- 1L
    for(p in to){
      if(asd[p]) next #already assigned - move to next
      if(subregions){
        sr <- srs[p]
        srrg <- srrgs[[sr]] #particle collection to which this point belongs
        cl <- sri[sr] - 1L + which({
          zsepsq <- if(TwoD) 0 else ((xyzm[p, z] - xyzm[srrg, z])^2)*(cdh/cdv)^2
          if(nan.flag && !TwoD) zsepsq[is.nan(zsepsq)] <- 0 # 0*Inf = NaN, but would like 0 here (when cdv = 0 and equal z values are compared)
          (xyzm[p, x] - xyzm[srrg, x])^2 + (xyzm[p, y] - xyzm[srrg, y])^2 + zsepsq <= cdh^2
        }) #those particle numbers within particle p's subregion to which p is close, including self
        grlabs[cl] <- lab
      }else{
        cl <- which({
          zsepsq <- if(TwoD) 0 else ((xyzm[p, z] - xyzm[, z])^2)*(cdh/cdv)^2
          if(nan.flag && !TwoD) zsepsq[is.nan(zsepsq)] <- 0
          (xyzm[p, x] - xyzm[, x])^2 + (xyzm[p, y] - xyzm[, y])^2 + zsepsq <= cdh^2
        }) #those particle numbers to which particle p is close, including self
        grlabs[cl] <- lab
      }
      lab <- lab + 1L
      asd[cl] <- TRUE #these are now assigned: don't check again (avoids clumping chains)
    }
  }else grlabs <- seq(nrow(xyzm))
  
  cl.time <- Sys.time()

  set(xyzm, NULL, "group", grlabs); setkey(xyzm, group)
  
  #aggregate with data table by group label
  xyzm <- if(TwoD){
    xyzm[, if(identical(length(x), 1L)){
      list(x = x, y = y, m = m)
    }else{
      xy <- llply(list(x = x, y = y), weighted.mean, m)
      c(xy, list(m = sum(m)))
    }, by = group]
  }else{
    xyzm[, if(identical(length(x), 1L)){
      list(x = x, y = y, z = z, m = m)
    }else{
      xyz <- llply(list(x = x, y = y, z = z), weighted.mean, m)
      c(xyz, list(m = sum(m)))
    }, by = group]
  }
  
  gr.time <- Sys.time()
  
  np <- nrow(xyzm)
  
  #which particles are smaller than mm, or which are the smallest particles if there is a set maximum
  if(maxbymm && mm <= 0) mm <- min(xyzm$m)
  smallf <- function() if(forcemmloop <<- maxbymm && np > maxnp){
    pflags <- logical(np)
    msrt <- sort.list(xyzm$m, method = "quick", na.last = NA)
    # the lightest np - maxnp particles are labelled as small
    pflags[msrt[1:(np - maxnp)]] <- TRUE; pflags
  }else xyzm$m < mm
  
  #coalesce particles with mass less than mm with nearest particles until no particles have mass less than mm
  while((mm > 0 && any(small <- smallf())) || forcemmloop){
    if(identical(np, 1L)) break #infinite loop would be possible here
    if(all(small)){
      #all particles are small - return one super-particle
      xyzm <- xyzm[, c(list(weighted.mean(x, m), weighted.mean(y, m)),
                       if(!TwoD) list(weighted.mean(z, m)), list(sum(m)))]
      break
    }
    
    # number of points left in each subregion
    if(subregions){
      npsrs <- mapply(function(start, end) length(unique(grlabs[start:end])), sri, srend - 1L)
      if(any(npsrs == 1L)) subregions <- FALSE else{ #if number of particles is getting too small, then there is danger of infinite loop, so stop using subregions
        sri <- vapply(1:nsr, function(n) sum(npsrs[seq_len(n - 1L)]), integer(1L)) + 1L
        srend <- vapply(1:nsr, function(n) sum(npsrs[seq_len(n)]), integer(1L)) + 1L
        grlabs <- unique(grlabs)
      }
    }

    # find closest neighbour to all small particles
    clst <- integer(np)
    
    # R subtlety: when reassigning by index (e.g. x[1:n] <- ...), the index argument (1:n) is evaluated *after* the expression on the right, so if, say, n changes during the evaluation of the right, the reassignment may not proceed as expected.  This is why smhold is used, to save the previous small vector before it is changed.
    
    smhold <- small
    if(subregions){
      clst[smhold] <- vapply(which(small), function(p){
        if(!small[p]) return(0L) #this one already in a group, so small flag has been reverted
        sr <- which(p >= sri & p < srend) #to which subregion does this point belong?
        srrg <- sri[sr]:(srend[sr] - 1L) #what is the entire set of points in this subregion?
        zsepsq <- if(TwoD) 0 else ((xyzm[p, z] - xyzm[srrg, z])^2)*(cdh/cdv)^2
        if(nan.flag && !TwoD) zsepsq[is.nan(zsepsq)] <- 0
        dsq <- (xyzm[p, x] - xyzm[srrg, x])^2 + (xyzm[p, y] - xyzm[srrg, y])^2 + zsepsq
        dsq[p - sri[sr] + 1L] <- NA # don't want it to find itself!
        cls <- which.min(dsq) + sri[sr] - 1L
        small[cls] <<- FALSE
        cls #returned value is the closest particle's index
      }, integer(1L)) #only calculates for small particles
    }else{
      clst[smhold] <- vapply(which(small), function(p){
        if(!small[p]) return(0L) #this one already in a group, so small flag has been reverted
        zsepsq <- if(TwoD) 0 else ((xyzm[p, z] - xyzm[, z])^2)*(cdh/cdv)^2
        if(nan.flag && !TwoD) zsepsq[is.nan(zsepsq)] <- 0
        dsq <- (xyzm[p, x] - xyzm[, x])^2 + (xyzm[p, y] - xyzm[, y])^2 + zsepsq
        dsq[p] <- NA # don't want it to find itself!
        cls <- which.min(dsq)
        small[cls] <<- FALSE
        cls #returned value is the closest particle's index
      }, integer(1L)) #only calculates for small particles
    }
    
    for(psm in which(small)){
      xyzm[clst[psm],] <- {
        grp <- xyzm[c(clst[psm], psm),] #should always be two rows

        #weighted average position
        xyznew <- if(TwoD){
          llply(grp[, list(x, y)], weighted.mean, grp$m)
        }else{
          llply(grp[, list(x, y, z)], weighted.mean, grp$m)
        }

        as.data.table(c(grp$group[1], xyznew, sum(grp$m)))
      }
      if(subregions) grlabs[psm] <- grlabs[clst[psm]]
    }

    #small particles are now moved to new homes, so to speak, so can delete their original entries
    #mass and centre of mass are conserved
    np <- nrow(xyzm <- xyzm[!small,])
  }

  mm.time <- Sys.time()

  timings <<- diff(c(st.time, subregions = sr.time, close = cl.time, coalescing = gr.time, minmass = mm.time))
  
  xyzm[, group := NULL]
  if(!is.null(cns)) setnames(xyzm, cns)
  if(na.rep && any(nas)) xyzm <- rbind(xyzm, xyzmna)
  rm(forcemmloop, envir = .GlobalEnv)
  return(xyzm)
}
