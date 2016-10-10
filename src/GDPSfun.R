# Christopher Barry, started on 20/06/2016 at University of Birmingham

# the super-functions used for the propagate and coalesce stages of GDPS

#propagation function----
#2D rotation matrix - produces 3D array with a phi given as a vector
rotmat2D <- function(phi, v = 1) structure({
  mapply(function(phi, v) v*c(cos(phi), sin(phi), -sin(phi), cos(phi)), phi, v)
}, dim = c(2L, 2L, length(phi)))

#function to calculate the mass distribution from one particle after one timestep
#returns a matrix with six columns: x, y, z, L, z_off, m
prop <- function(state, t.new, Delta.t, newcbf, por, statei = NULL, Rf = 1, sorb = FALSE){
  t.old <- t.new - Delta.t
  
  #sorb...
  
  #equilibrium sorption is assumed here
  if(sorb){
    #new sorption
    newstatei <- state[, list(x, y, L, zo, m = m*(1 - 1/Rf))]
    state[, m := m/Rf]
    
    #sorbed mass giving back
    newstate <- statei[, list(x, y, L, zo, m = m/Rf)]
    statei[, m := m*(1 - 1/Rf)]
    
    state <- rbind(state, newstate)
    statei <- rbind(statei, newstatei)
    rm(newstate, newstatei)
  }
  np <- nrow(state)
  
  if(decaysorbed) statei[, m := m*exp(-lambda*Delta.t)]
  
  #advect...
  
  write(rsptxt(t.new - MFt0, newcbf), paste0(dmrt, ".rsp"))
  write(dmoc.PTR(cbind(as.matrix(state[, list(x - MFxy0[1L], y - MFxy0[2L], zo, L)]), t.old - MFt0)), paste0(dmrt, ".ptr"))
  if(np > MXP.def) write(dattxt(por, np, dis), paste0(dmrt, ".dat")) #only rewrite when more particles are required (saves memory otherwise)
  
  #run MODPATH
  tm <- Sys.time()
  system(paste0(mpexe, " ", dmrt, ".rsp"), intern = TRUE, wait = TRUE, show.output.on.console = FALSE)
  if(file.mtime("pathline") < tm - 5) stop("MODPATH failed") #small buffer to allow for very fast run times
  
  #MODPATHmass add on
  ptldf <- fread("pathline", skip = 1L)
  reld <- seq_len(np) %in% ptldf[, V1] #which particles released? (i.e. started in active region)
  ptldf[, V6 := V6 + MFt0] # time correction (column not renamed at this point)
  mt <- ptl.masstrack(ptldf, gwdata, wtop, por, state$m, retain.storage = TRUE,
                      outflux = TRUE, loss = TRUE, tlumpoutflux = TRUE, end.t = t.new, linear.decay = lambda)
  rm(ptldf) #relieve memory
  
  #collect pathline data and extract necessary columns
  #positions converted to absolute grid reference at this stage
  #key value ptlno is retained
  xyzm <- mt$traces[, llply(.SD, tail, 1L), by = ptlno, .SDcols = c("x", "y", "z", "z_off", "L", "m")]
  setnames(xyzm, "z_off", "zo")
  xyzm[, c("x", "y") := list(x + MFxy0[1L], y + MFxy0[2L])] # convert back to global co-ordinates
  
  #sinks...
  
  #update outflux and loss data
  if(any(mt$outflux != 0)) fluxout[[tPt]] <<- {
    im <- which(mt$outflux != 0, arr.ind = TRUE)
    vals <- mt$outflux[im]
    dt <- data.table(tPt, im, vals/Delta.t)
    setnames(dt, c("ts", "C", "R", "L", "J_out"))
  }
  massloss[tPt - 1L] <<- sum(mt$loss) + massloss[tPt - 1L]
  
  if(all(xyzm$m == 0)) return(xyzm) #note will return if nrow(xyzm) == 0, which is good (all(logical(0L)) = TRUE)
  
  #update masses based on what was lost to sinks
  #the following line is a TEMPORARY FIX to the problem of particles being dispersed into no flow cells - mass will be lost in this way
  #a proper fix would require a different imprint for transport along boundaries
  #that said, it's possible to get some diffusive loss into aquitard material
  #however, you'd expect to get some back as well over longer timescales (back-diffusion), which this doesn't show
  state <- state[reld,]
  
  xy.old <- state[, list(x, y)] #original positions
  delta.xy <- xyzm[, list(x, y)] - xy.old #not got 3D rotation yet! (not sure that I'd want it anyway - shouldn't dispersion tensor be in line with layers)
  v <- if(vdepD) sqrt(delta.xy$x^2 + delta.xy$y^2)/Delta.t else 1
  
  #disperse...
  
  #a future development may be to allow sorbed mass to diffuse
  
  #new grid
  dldt <- if(ThreeDD){
    2*dispC*(Delta.t*c(DL, DT, DV))^0.5
  }else{
    2*dispC*(Delta.t*c(DL, DT))^0.5
  }
  phi <- atan2(delta.xy$y, delta.xy$x) # orientation of displacement (0 along x axis, increasing anti-clockwise)
  
  #finds the positions relative to the original particle (without rotation yet - in the longitudinal and transverse dimensions) for which the pre-determined FD imprint is true
  relpos <- if(ThreeDD){
    cbind(dldt[1L]*c(-1, 0, 0, 0, 0, 0, 1),
          dldt[2L]*c(0, -1, 0, 0, 0, 1, 0),
          dldt[3L]*c(0, 0, -1, 0, 1, 0, 0))
  }else{
    cbind(dldt[1L]*c(-1, 0, 0, 0, 1),
          dldt[2L]*c(0, -1, 0, 1, 0))
  }
  
  #applies the displacement and rotation according to the particles' movements
  rots <- rotmat2D(-phi, v) # note that v = 1 if not vdepD; -phi is used because the co-ordinates are on the LHS of %*% (%*% is not commutative)
  mfts <- cellref.loc(t.new, c(0, gwdata$time) + MFt0) # MODFLOW timestep number
  if(is.na(mfts)){
    # this covers case for which the last tval is the same as the end time of the MODFLOW model
    # because cellref.loc uses >= and <, NA is returned in this case, so the time
    mfts <- cellref.loc(t.new - Delta.t/100, c(0, gwdata$time) + MFt0) # MODFLOW timestep number re-assessed
  }
  
  newpos <- xyzm[, {
    #find the rotation matrix
    rot <- if(ThreeDD){
      # this 3D rotation matrix is always about a vertical axis
      # the z dispersion is scaled so that it operates on the z offset rather than z, truer to the way MODFLOW/ MODPATH represent dipping layers (no diagonal connection)
      rtmp <- diag(v, 3L); rtmp[1:2, 1:2] <- rots[,, ptlno]
      rtmp
    }else rots[,, ptlno]
    
    #rotate
    r <- relpos %*% rot
    
    #shift back to real co-ordinates
    xnew <- r[, 1L] + x
    ynew <- r[, 2L] + y
    
    C <- cellref.loc(xnew, gwdata$gccs + MFxy0[1L], FALSE)
    R <- cellref.loc(ynew, gwdata$grcs + MFxy0[2L], TRUE)
    bot <- gwdata$elev[cbind(C, R, L + 1L)]
    top <- wtop[cbind(C, R, L, mfts)]
    
    # Explanation:
    # if 2D dispersion, then simple. z offsets, rather than z values are maintained
    # For 3D dispersion
    # 1. the horizontally dispersed particles retain z offset values
    # 2. the vertically dispersed particles are dispersed by z, not by z offset, and then z offset and layer are subsequently calculated for these two particles
    # 3. with vertical dispersion, there is an additional option to retain any mass lost from the top or bottom of the model by mirroring it back in (according to the method of images, valid for linear processes such as dispersion)
    # note that with 3D dispersion, new particles 3 and 5 are those that are dispersed vertically
    
    Lnew <- rep(L, if(ThreeDD) 7L else 5L)
    zonew <- if(ThreeDD){
      c(zo, zo, NA, zo, NA, zo, zo)
    } else zo
    
    if(ThreeDD){
      znew3D <- r[c(3L, 5L), 3L] + z
      NLAY <- dim(gwdata$elev)[3L] - 1L
      nwet <- sum(!is.na(wtop[C[4L], R[4L],, mfts]))
      #depth divides
      # the water table is repeated for the number of dry or partially saturated cells in the column and then the layer bottoms are used
      # the length of ddivs will always be NLAY + 1
      # should ensure that a particle is never labelled with a layer that is dry at its location
      ddivs <- c(rep(wtop[C[4L], R[4L], NLAY - nwet + 1L, mfts], NLAY - nwet + 1L),
                 gwdata$elev[C[4L], R[4L], seq(to = NLAY, by = 1L, length.out = nwet) + 1L])
      Lnew[c(3L, 5L)] <- cellref.loc(znew3D, rev(ddivs), TRUE)
      bot <- ddivs[Lnew[c(3L, 5L)] + 1L]
      top <- ddivs[Lnew[c(3L, 5L)]]
      zonew[c(3L, 5L)] <- (znew3D - bot)/(top - bot)
      
      #retain.vloss will reflect back any mass that has been lost out of the model vertically
      while(retain.vloss && any(!Lnew %in% 1:(nlay <- dim(gwdata$data)[3L]))){
        # determine whether and how far particles are above the water column or below the model domain
        # the overshoot is relative to the highest wtop for in non-dry cells at the xy location
        wheight <- max(wtop[C[4L], R[4L],, mfts], na.rm = TRUE)
        overshoot <- znew3D - wheight #only relevant for particles that are too high
        undershoot <- gwdata$elev[C[4L], R[4L], nlay + 1L] - znew3D #ditto for too low
        
        #new z values by reflection
        znew3D[overshoot > 0] <- wheight - overshoot[overshoot > 0]
        znew3D[undershoot > 0] <- gwdata$elev[C[4L], R[4L], nlay + 1L] + undershoot[undershoot > 0]
        
        #re-determine L and z offset
        Lnew[c(3L, 5L)] <- cellref.loc(znew3D, rev(ddivs), TRUE)
        bot <- ddivs[Lnew[c(3L, 5L)] + 1L]
        top <- ddivs[Lnew[c(3L, 5L)]]
        zonew[c(3L, 5L)] <- (znew3D - bot)/(top - bot)
        
        # Very occasionally this loop may need to repeat if particles are reflected out the other side.  But such models are probably not very well built.  Generally speaking, the `while` should be thought of as `if`, as it should only be needed once.
      }
    }
    
    list(x = xnew, y = ynew, L = Lnew, zo = zonew, m = {
      #mass redistributed to the new particles
      m*(if(ThreeDD) c(edge, edge, edge, centre, edge, edge, edge) else c(edge, edge, centre, edge, edge))
    })
  }, by = ptlno]
  
  newpos[, ptlno := NULL]
  
  setnames(newpos, colord)
  
  if(sorb) return(list(newpos, statei)) else return(newpos)
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
  
  #coalescing is done by z offset, rather than z
  #this is more true to the way that both MODFLOW and MODPATH treat dipping geological units - there is no diagonal cell connection, and it avoids particles near the bottom or top of the model water column being moved out of the model domain by coalescence
  #lmzo means layer minus z offset
  #lthk is the layer thickness of the active particles' cells: the weighted average is taken and cd["v"] is divided by this to get the dimensionless vertical coalescing radius.
  lmzo <- state[, L - zo]
  if(any(state$zo == 1, na.rm = TRUE)) lmzo[state$zo == 1] <- state[zo == 1, L] - 0.9999 #avoids L = 0 anywhere, unlikely to be needed though as MODPATH 5 corrects this already
  C <- cellref.loc(state$x, gwdata$gccs + MFxy0[1L], F) #column numbers
  R <- cellref.loc(state$y, gwdata$grcs + MFxy0[2L], T) #row numbers
  mfts <- cellref.loc(t, c(0, gwdata$time) + MFt0, F) # timestep number; need t0 here?
  lthk <- state[, {
    imtop <- cbind(C, R, L, mfts); imbot <- cbind(C, R, L + 1L)
    wtop[imtop] - gwdata$elev[imbot]
  }]
  cdvdls <- cd["v"]/weighted.mean(lthk, state$m, na.rm = TRUE) #NA may arise if particles are outside model bounds - this instance is dealt with in the prop function and should not really be considered as an error as it is a naturally possible occurrence
  rm(C, R, lthk)
  
  #perform coalesce
  state <- coalesce3Dm(state[, list(x = x, y = y, lmzo = lmzo, m = m)], cd["h"], cdvdls, mm, subregions = T, maxnp = maxp, na.rep = TRUE)
  
  #work out from the lmzo values what layers and z offsets should be given to the particles
  L <- as.integer(ceiling(state$lmzo)) #layer number is integer part of lmzo + 1
  zo <- 1 - (state$lmzo %% 1) #z offset is mantissa of lmzo
  
  #no need to calculate z now...
  
  #delete lmzo column and add the layer and z offset columns
  set(state, NULL, c("lmzo", "zo", "L"), list(NULL, zo, L))
  setcolorder(state, colord)
  
  return(state)
}

