# Christopher Barry, started on 20/06/2016 at University of Birmingham

# the super-functions used for the propagate and coalesce stages of GDPS

#propagation function----

if(!exists("Ndisppairs")) Ndisppairs <- 2L

#function to calculate the mass distribution from one particle after one timestep
#returns a matrix with six columns: x, y, z, L, z_off, m
prop <- function(state, t.new, Delta.t, newcbf, por, statei = NULL, Rf = 1, sorb = FALSE){
  t.old <- t.new - Delta.t
  
  ### sorb ----
  
  # equilibrium sorption is assumed here
  if(sorb){
    # new sorption
    newstatei <- state[, list(x, y, L, zo, m = m*(1 - 1/Rf))]
    state[, m := m/Rf]
    
    # sorbed mass giving back
    newstate <- statei[, list(x, y, L, zo, m = m/Rf)]
    statei[, m := m*(1 - 1/Rf)]
    
    state <- rbind(state, newstate)
    statei <- rbind(statei, newstatei)
    rm(newstate, newstatei)
  }
  np <- nrow(state)
  
  if(decaysorbed && lambda){
    degraded[tPt] <<- sum(statei$m*(1 - lambda*Delta.t))
    statei[, m := m*exp(-lambda*Delta.t)]
  }
  
  ### advect ----
  
  write(rsptxt(t.new - MFt0, newcbf), paste0(dmrt, ".rsp"))
  write(gd.PTR(state[, list(x = x - MFxy0[1L],
                            y = y - MFxy0[2L],
                            L, zo)], t.old - MFt0),
        paste0(dmrt, ".ptr"))
  if(np > MXP.def) write(dattxt(por, np, dis), paste0(dmrt, ".dat")) #only rewrite when more particles are required (saves memory otherwise)
  
  ## run MODPATH --
  tm <- Sys.time()
  system(paste0(mpexe, " ", dmrt, ".rsp"), intern = TRUE,
         wait = TRUE, show.output.on.console = FALSE)
  # tests whether MODPATH has run successfully by checking whether there is
  #  a new pathline file; there is a five-second grace period in case of
  #  very fast run times and slight timing inconsistencies
  if(!file.exists("pathline") || file.mtime("pathline") < tm - 5){
    stop("MODPATH failed")
  }
  
  ## MassTrack --
  ptldf <- fread("pathline", skip = 1L)
  setnames(ptldf, PTL.headers[1:10])
  
  # which particles released? (i.e. started in active region)
  reld <- seq_len(np) %in% ptldf$ptlno
  
  # time correction
  ptldf[, t := t + MFt0]
  mt <- ptl.masstrack2(ptldf, gwdata, wtop, por, state$m, retain.storage = TRUE,
                       outflux = TRUE, loss = TRUE, tlumpoutflux = TRUE,
                       react.loss = TRUE, outflux.array = FALSE,
                       end.t = t.new, linear.decay = lambda)
  rm(ptldf)
  
  # collect pathline data and extract necessary columns
  # positions converted to absolute grid reference at this stage
  # key value ptlno is retained
  xyzm <- mt$traces[, list(ptlno, x, y, z, z_off, L, t, m)]
  setnames(xyzm, "z_off", "zo")
  
  # convert back to global co-ordinates
  xyzm[, c("x", "y") := list(x + MFxy0[1L], y + MFxy0[2L])]
  
  ### sinks ----
  
  ## update outflux and loss data --
  if(any(mt$traces$ml != 0)){
    fluxout[[tPt]] <<- mt$traces[, list(ts = tPt, J_out = sum(ml)/Delta.t),
                                 by = c("C", "R", "L")]
  }
  
  massloss[["inactive"]][tPt] <<- sum(mt$loss)
  if(lambda) degraded[tPt] <<- degraded[tPt] + mt$traces[, sum(mrl)]
  
  
  ## update masses based on what was lost to sinks --
  state <- state[reld,]
  
  ## average travel speed of each particle --
  # note sum(double(0L)) = 0, so paths with one entry will return v = 0 (in
  #  case of instant capture) - OK
  v <- if(vdepD){
    xyzm[, sum(sqrt(diff(x)^2 + diff(y)^2))/diff(range(t)), by = ptlno]$V1
  }else 1
  
  ## average direction of travel of each particle --
  # atan2(0, 0) = 0, in case of instant capture - OK
  phi <- atan2(xyzm[, y[.N] - y[1L], by = ptlno]$V1,
               xyzm[, x[.N] - x[1L], by = ptlno]$V1)
  
  # get end points of pathlines as new positions (.N is length of subset)
  # rows are duplicated for the dispersion step
  xyzm <- xyzm[, .SD[rep(.N, 2*Ndisppairs)],
               by = ptlno, .SDcols = c("x", "y", "z", "zo", "L", "t", "m")]
  
  if(all(xyzm$m == 0)){
    xyzm <- xyzm[, .SD[.N],
                 by = ptlno, .SDcols = c("x", "y", "zo", "L", "m")]
    xyzm[, ptlno := NULL]
    setcolorder(xyzm, colord)
    return(xyzm) 
    # note will return if nrow(xyzm) = 0, which is good (all(logical(0L)) = TRUE)
  }
  
  # collapse ptlno so that it is a sequence from 1 to Npart; avoids
  #  complications in using ptlno to index if some of the particles have
  #  been lossed or completely drained
  xyzm[, ptlno := .GRP, by = ptlno]
  np <- uniqueN(xyzm$ptlno)
  
  ## disperse ----
  
  #a future development may be to allow sorbed mass to diffuse
  
  #new grid
  dldt <- if(ThreeDD){
    2*(Delta.t*c(DL, DT, DV))^0.5
  }else{
    2*(Delta.t*c(DL, DT))^0.5
  }
  
  #finds the positions relative to the original particle (without rotation yet - in the longitudinal and transverse dimensions) using a random normal scattering in a random orientation
  npd <- np*Ndisppairs
  relpos <- if(ThreeDD){
    #distance
    dispersion <- rnorm(npd)
    dispersion <- dispersion[rep(seq_len(npd), each = 2L)]
    
    #in horizontal plane
    rotation1 <- runif(npd, -pi, pi)
    rotation1 <- rotation1[rep(seq_len(npd), each = 2L)] + c(0, pi)
    
    #from horizontal plane
    rotation2 <- runif(npd, -pi/2, pi/2)
    rotation2 <- rotation2[rep(seq_len(npd), each = 2L)] + c(0, pi)
    
    cbind(dldt[1L]*dispersion*cos(rotation1)*cos(rotation2),
          dldt[2L]*dispersion*sin(rotation1)*cos(rotation2),
          dldt[3L]*dispersion*sin(rotation2))
  }else{
    #distance
    dispersion <- rnorm(npd)
    dispersion <- dispersion[rep(seq_len(npd), each = 2L)]
    
    #horizontal plane
    rotation <- runif(npd, -pi, pi)
    rotation <- rotation[rep(seq_len(npd), each = 2L)] + c(0, pi)
    
    cbind(dldt[1L]*dispersion*cos(rotation),
          dldt[2L]*dispersion*sin(rotation))
  }
  
  ## applies the displacement and rotation according to the particles'
  ##  movements --
  # note that v = 1 if not vdepD; -phi is used because the co-ordinates are
  #  on the LHS of %*% (%*% is not commutative)
  # if 3D dispersion, then the rotation is applied about a vertical axis,
  #  as if the velocity is purely horizontal, which is truer to the way
  #  MODFLOW respresents flow within a layer
  rots <- if(ThreeDD) replicate(np, diag(1, 3L), "array") else array(0, c(2L, 2L, np))
  
  rots[1:2, 1:2,] <- mapply(function(phi, v){
    v*matrix(c(cos(-phi), sin(-phi), -sin(-phi), cos(-phi)), 2L, 2L)
  }, phi, v)
  
  xyzm[, c("x", "y", "z", "m") := {
    rotpos <- relpos[.I,] %*% rots[,, ptlno]
    list(rotpos[, 1L] + x,
         rotpos[, 2L] + y,
         z + if(ThreeDD) rotpos[, 3L] else 0,
         m/(2*Ndisppairs))
  }, by = ptlno]
  
  # MODFLOW timestep number
  mfts <- cellref.loc(t.new, c(0, gwtime) + MFt0)
  if(is.na(mfts)){
    # this covers case for which the last tval is the same as the end time
    #  of the MODFLOW model
    # because cellref.loc uses >= and <, NA is returned in this case
    # MODFLOW timestep number re-assessed
    mfts <- cellref.loc(t.new - Delta.t/100, c(0, gwtime) + MFt0)
  }
  
  if(ThreeDD) xyzm[, c("zo", "L") := {
    C <- cellref.loc(x, gccs + MFxy0[1L], FALSE)
    R <- cellref.loc(y, grcs + MFxy0[2L], TRUE)
    bot <- nc.imtx(gwdata, "elev", cbind(C, R, L + 1L))
    top <- nc.imtx(wtop, "wtop", cbind(C, R, L, mfts))
    
    NLAY <- dis$extent["NLAY"]
    nwet <- mapply(function(ct, rt){
      sum(!is.na(var.get.nc(wtop, "wtop",
                            c(ct, rt, 1L, mfts), c(1L, 1L, NLAY, 1L))))
    }, C, R)
    #depth divides
    # the water table is repeated for the number of dry or partially saturated cells in the column and then the layer bottoms are used
    # the length of ddivs will always be NLAY + 1
    # should ensure that a particle is never labelled with a layer that is dry at its location
    ddivs <- mapply(function(ct, rt){
      c(rep(wtop[ct, rt, NLAY - nwet + 1L, mfts], NLAY - nwet + 1L),
        var.get.nc(gwdata, "elev", c(ct, rt, NLAY - nwet + 1L),
                   c(1L, 1L, nwet)))
    }, C, R)
    Lnew <- apply(ddivs, 2L, function(dd) cellref.loc(z, rev(dd), TRUE))
    
    zonew <- punif(z, dd[cbind(Lnew, 1:2)], dd[cbind(Lnew + 1L, 1:2)])
    
    # retain vloss algorithm to come
  }]
  
  # newpos <- xyzm[, {
  #   #find the rotation matrix
  #   rot <- if(ThreeDD){
  #     # this 3D rotation matrix is always about a vertical axis
  #     # the z dispersion is scaled so that it operates on the z offset rather than z, truer to the way MODFLOW/ MODPATH represent dipping layers (no diagonal connection)
  #     rtmp <- diag(v, 3L); rtmp[1:2, 1:2] <- rots[,, ptlno]
  #     rtmp
  #   }else rots[,, ptlno]
  #   
  #   #rotate
  #   r <- relpos %*% rot
  #   
  #   #shift back to real co-ordinates
  #   xnew <- r[, 1L] + x
  #   ynew <- r[, 2L] + y
  #   
  #   C <- cellref.loc(xnew, gwdata$gccs + MFxy0[1L], FALSE)
  #   R <- cellref.loc(ynew, gwdata$grcs + MFxy0[2L], TRUE)
  #   bot <- gwdata$elev[cbind(C, R, L + 1L)]
  #   top <- wtop[cbind(C, R, L, mfts)]
  #   
  #   # Explanation:
  #   # if 2D dispersion, then simple. z offsets, rather than z values are maintained
  #   # For 3D dispersion
  #   # 1. the horizontally dispersed particles retain z offset values
  #   # 2. the vertically dispersed particles are dispersed by z, not by z offset, and then z offset and layer are subsequently calculated for these two particles
  #   # 3. with vertical dispersion, there is an additional option to retain any mass lost from the top or bottom of the model by mirroring it back in (according to the method of images, valid for linear processes such as dispersion)
  #   # note that with 3D dispersion, new particles 3 and 5 are those that are dispersed vertically
  #   
  #   Lnew <- rep(L, 2L)
  #   zonew <- if(ThreeDD){
  #     c(zo, zo, NA, zo, NA, zo, zo)
  #   } else zo
  #   
  #   if(ThreeDD){
  #     znew3D <- r[c(3L, 5L), 3L] + z
  #     NLAY <- dim(gwdata$elev)[3L] - 1L
  #     nwet <- sum(!is.na(wtop[C[4L], R[4L],, mfts]))
  #     #depth divides
  #     # the water table is repeated for the number of dry or partially saturated cells in the column and then the layer bottoms are used
  #     # the length of ddivs will always be NLAY + 1
  #     # should ensure that a particle is never labelled with a layer that is dry at its location
  #     ddivs <- c(rep(wtop[C[4L], R[4L], NLAY - nwet + 1L, mfts], NLAY - nwet + 1L),
  #                gwdata$elev[C[4L], R[4L], seq(to = NLAY, by = 1L, length.out = nwet) + 1L])
  #     Lnew[c(3L, 5L)] <- cellref.loc(znew3D, rev(ddivs), TRUE)
  #     bot <- ddivs[Lnew[c(3L, 5L)] + 1L]
  #     top <- ddivs[Lnew[c(3L, 5L)]]
  #     zonew[c(3L, 5L)] <- (znew3D - bot)/(top - bot)
  #     
  #     #retain.vloss will reflect back any mass that has been lost out of the model vertically
  #     while(retain.vloss && any(!Lnew %in% 1:(nlay <- dim(gwdata$data)[3L]))){
  #       # determine whether and how far particles are above the water column or below the model domain
  #       # the overshoot is relative to the highest wtop for in non-dry cells at the xy location
  #       wheight <- max(wtop[C[4L], R[4L],, mfts], na.rm = TRUE)
  #       overshoot <- znew3D - wheight #only relevant for particles that are too high
  #       undershoot <- gwdata$elev[C[4L], R[4L], nlay + 1L] - znew3D #ditto for too low
  #       
  #       #new z values by reflection
  #       znew3D[overshoot > 0] <- wheight - overshoot[overshoot > 0]
  #       znew3D[undershoot > 0] <- gwdata$elev[C[4L], R[4L], nlay + 1L] + undershoot[undershoot > 0]
  #       
  #       #re-determine L and z offset
  #       Lnew[c(3L, 5L)] <- cellref.loc(znew3D, rev(ddivs), TRUE)
  #       bot <- ddivs[Lnew[c(3L, 5L)] + 1L]
  #       top <- ddivs[Lnew[c(3L, 5L)]]
  #       zonew[c(3L, 5L)] <- (znew3D - bot)/(top - bot)
  #       
  #       # Very occasionally this loop may need to repeat if particles are reflected out the other side.  But such models are probably not very well built.  Generally speaking, the `while` should be thought of as `if`, as it should only be needed once.
  #     }
  #   }
  #   
  #   list(x = xnew, y = ynew, L = Lnew, zo = zonew, m = {
  #     #mass redistributed to the new particles
  #     m*(if(ThreeDD) c(edge, edge, edge, centre, edge, edge, edge) else c(edge, edge, centre, edge, edge))
  #   })
  # }, by = ptlno]
  # 
  # newpos[, ptlno := NULL]
  
  xyzm[, c("z", "ptlno", "t") := NULL]
  
  setcolorder(xyzm, colord)
  
  if(sorb) return(list(xyzm, statei)) else return(xyzm)
}

#coalescing function----
#finds particles which are close together and lumps them to new particle at centre of mass - essential moderater of particle count
#additionally coalesces particles below a certain mass with nearest non-insignificant particle - this happens after the initial coalescence and avoids spending time tracking particles with negligible mass whilst still ensuring conservation of mass
#expected columns: 1 - x, 2 - y, 3 - z, 4 - zo, 5 - L, 6 - m
coalesce <- function(state, cd, mm = 0, t){
  # in the case that all mass is depleted, return a zero-row data table
  if(all(state$m == 0)) return(data.table(ts = integer(0L),
                                          x = double(0L),
                                          y = double(0L),
                                          L = integer(0L),
                                          zo = double(0L),
                                          m = double(0L)))
  
  # coalescing is done by z offset, rather than z
  # this is more true to the way that both MODFLOW and MODPATH treat
  #  dipping geological units - there is no diagonal cell connection, and
  #  it avoids particles near the bottom or top of the model water column
  #  being moved out of the model domain by coalescence
  
  ## determine lmzo --
  # lmzo means layer minus z offset
  # lthk is the layer thickness of the active particles' cells: the
  #  weighted average is taken and cd["v"] is divided by this to get the
  #  dimensionless vertical coalescing radius.
  lmzo <- state[, L - zo]
  
  # avoids L = 0 anywhere, unlikely to be needed though as MODPATH 5
  #  corrects this already
  if(any(state$zo == 1, na.rm = TRUE)) lmzo[state$zo == 1] <- state[zo == 1, L] - 0.9999
  # column numbers
  C <- cellref.loc(state$x, gccs + MFxy0[1L], F)
  # row numbers
  R <- cellref.loc(state$y, grcs + MFxy0[2L], T)
  # timestep number; need t0 here?
  mfts <- cellref.loc(t, c(0, gwtime) + MFt0, F)
  lthk <- state[, {
    imtop <- cbind(C, R, L, mfts); imbot <- cbind(C, R, L + 1L)
    nc.imtx(wtop, "wtop", imtop) - nc.imtx(gwdata, "elev", imbot)
  }]
  
  # NA may arise if particles are outside model bounds - this instance is
  #  dealt with in the prop function and should not really be considered as
  #  an error as it is a naturally possible occurrence
  cdvdls <- cd["v"]/weighted.mean(lthk, state$m, na.rm = TRUE)
  rm(C, R, lthk)
  
  ## perform coalesce --
  state <- coalesce3Dm(state[, list(x = x, y = y, lmzo = lmzo, m = m)], cd["h"], cdvdls, mm, subregions = T, maxnp = maxp, na.rep = TRUE)
  
  ## work out from the lmzo values what layers and z offsets should be
  ##  given to the particles -- 
  # layer number is integer part of lmzo + 1
  L <- as.integer(ceiling(state$lmzo))
  
  # z offset is mantissa of lmzo
  zo <- 1 - (state$lmzo %% 1)
  
  # no need to calculate z now...
  
  ## delete lmzo column and add the layer and z offset columns --
  set(state, NULL, c("lmzo", "zo", "L"), list(NULL, zo, L))
  setcolorder(state, colord)
  
  return(state)
}

