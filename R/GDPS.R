#script file that executes the grid-decoupled plume simulation


# preamble ----------------------------------------------------------------

sim.start <- Sys.time()
library("plyr")
library("rlist")
library("stringr")
library("abind")
library("data.table")
library("ks")
library("sp")
source(paste0(gendir, "MFread.R"))
source(paste0(gendir, "MT3D.R")) #for the array writing function RIARRAY
source(paste0(gendir, "MODPATHmass.R"))

#functions for propagate and coalesce stages
source("C:/Users/cjb309/Documents/GitHub/coalesce/R/coalesce.R")
source("C:/Users/cjb309/Documents/GitHub/GDPS/R/GDPSfun2.R")

od <- getwd(); setwd(mfdir)
if(!exists("sorb")) sorb <- FALSE


# perform checks before starting ------------------------------------------

# check that there is correct number of release rate functions
if(length(rel.fun) != nrow(xy0)) stop("The number of release functions in rel.fun (which should",
                                      "be a list of functions), should equal the number of release",
                                      "points.\n",
                                      "Currently:\nlength(rel.fun) == ", length(rel.fun), "\n",
                                      "nrow(xy0) == ", nrow(xy0))


# load groundwater data ---------------------------------------------------

# new system to keep connection with data files and load in only what is
#  needed for the current time step; different from previous lRAM because
#  no need to start again in file each time

# set up connections
hdcon <- file(paste0(mfrt, ".hds"), "rb")
flcon <- file(paste0(mfrt, ".cbb"), "rb")

# get co-ordinate data
gwdata <- list()
if(reload.DIS || !exists("dis")) dis <- read.DIS(paste0(mfrt, ".dis"))
if(reload.BAS || !exists("bas")) bas <- read.BAS(paste0(mfrt, ".bas"), dis)

gwdata$gccs <- if(identical(length(dis$DELR), 1L)){
  seq(0, by = dis$DELR, length.out = dis$extent["NCOL"] + 1L)
}else cumsum(c(0, dis$DELR))
gwdata$grcs <- if(identical(length(dis$DELC), 1L)){
  seq(0, by = dis$DELC, length.out = dis$extent["NROW"] + 1L)
}else cumsum(c(0, dis$DELC))
gwdata$elev <- if(is.null(dim(dis$elev))){
  sapply(dis$elev, array, dis$extent[c("NCOL", "NROW")],
         simplify = "array")
}else dis$elev

gwtime <- with(dis$sps, {
  data.frame(i = 1:sum(NSTP),
             t = c(Map(function(Dt, N, m, t.end, sp){
               dt1 <- Dt/sum(m^(0:(N - 1L)))
               ts <- cumsum(dt1*m^(0:(N - 1L))) + t.end - Dt
               names(ts) <- paste(sp, 1:N, sep = "_")
               ts
             }, PERLEN, NSTP, TSMULT, cumsum(PERLEN), 1:length(PERLEN)),
             recursive = TRUE))
})

# function to read in new groundwater data for the timestep
# does not require array to be passed in by value, so saving memory
# this function assumes that all the MODFLOW output data is written in
#  order of model time (i.e. the time steps are accessed sequentially)
renew.gwdata <- function(curr.mfts, new.mfts, readFLF, readS){
  # overlap; probable
  if(any(new.mfts %in% curr.mfts)){
    gwdata$data[,,, 1:sum(new.mfts %in% curr.mfts),] <<-
      gwdata$data[,,, which(curr.mfts %in% new.mfts),]
    
    wtop[,,, 1:sum(new.mfts %in% curr.mfts)] <<-
      wtop[,,, which(curr.mfts %in% new.mfts)]
  }
  
  # does array need expanding?
  if(length(new.mfts) > dim(gwdata$data)[4L]){
    # MODFLOW data array
    tmp <- double(prod(dis$extent[c("NCOL", "NROW", "NLAY")])*
                    length(new.mfts)*(3L + readFLF + readS))
    dim(tmp) <- c(dis$extent[c("NCOL", "NROW", "NLAY")],
                  length(new.mfts), unname(3L + readFLF + readS))
    dimnames(tmp)[4:5] <- list(names(gwdata$time)[new.mfts],
                               c("Head",
                                 "FlowRightFace",
                                 "FlowFrontFace",
                                 if(readFLF) "FlowLowerFace",
                                 if(readS) "Storage"))
    tmp[,,, 1:sum(new.mfts %in% curr.mfts),] <-
      gwdata$data[,,, 1:sum(new.mfts %in% curr.mfts),]
    gwdata$data <<- tmp
    rm(tmp)
    
    # water top
    tmp <- double(prod(dis$extent[c("NCOL", "NROW", "NLAY")])*
                    length(new.mfts))
    dim(tmp) <- c(dis$extent[c("NCOL", "NROW", "NLAY")],
                  length(new.mfts))
    dimnames(tmp)[[4L]] <- names(gwdata$time)[new.mfts]
    tmp[,,, 1:sum(new.mfts %in% curr.mfts)] <-
      wtop[,,, 1:sum(new.mfts %in% curr.mfts)]
    wtop <<- tmp
    rm(tmp)
  }
  
  # know when to stop
  last.L <- dis$extent["NLAY"]
  last.mfts <- last(new.mfts)
  
  # read in head data
  mfL <- mfts <- 0L
  cat("reading head data:\n")
  cat(str_dup("/", sum(new.mfts %in% curr.mfts)))
  while(mfts != last.mfts || mfL != last.L){
    # time step, stress period
    tssp <- readBin(hdcon, "integer", 2L, 4L)
    mfts <- gwtime[paste(rev(tssp), collapse = "_"), "i"]
    
    # skip simulation time
    readBin(hdcon, "double", 2L, 4L)
    
    # parameter type, should be Head
    arty <- nicearname(readChar(hdcon, 16L))
    stopifnot(arty == "Head")
    
    # skip CR count, get layer number
    mfL <- readBin(hdcon, "integer", 3L, 4L)[3L]
    
    # which time step index to store the next array (if at all)
    tspos <- which(new.mfts == mfts)
    
    if(length(tspos)){
      # put data into array
      gwdata$data[,, mfL, tspos, arty] <<-
        readBin(hdcon, "double", prod(dis$extent[c("NCOL", "NROW")]), 4L)
      
      # determine wtop
      # wtop (water top) is the top of the water in the cell
      # it is the lower of Head or the cell top
      # in MODPATH, z offset is defined as the proportionate vertical
      #  position between the cell base and the water top
      wtop[,, mfL, tspos] <<-
        ifelse(dis$elev[,, mfL] <= gwdata$data[,, mfL, tspos, "Head"],
               dis$elev[,, mfL], gwdata$data[,, mfL, tspos, "Head"])
      if(mfL == last.L) cat("|")
    }else{
      # this time step is before the requested range (normally only for the
      #  first time round)
      readBin(hdcon, "double", prod(dis$extent[c("NCOL", "NROW")]), 4L)
      if(mfL == last.L) cat(".")
    }
  }; cat("\n")
  
  # read in flow data
  mfts <- 0L
  cat("reading flow data:\n")
  cat(str_dup("/", sum(new.mfts %in% curr.mfts)))
  while(mfts != last.mfts){
    # time step, stress period
    tssp <- readBin(flcon, "integer", 2L, 4L)
    mfts <- gwtime[paste(rev(tssp), collapse = "_"), "i"]
    
    # parameter type
    arty <- nicearname(readChar(flcon, 16L))
    
    # skip CRL count
    readBin(flcon, "integer", 3L, 4L)
    
    # which time step index to store the next array (if at all)
    tspos <- which(new.mfts == mfts)
    
    # if this is a required parameter, put into array, else skip
    if(length(tspos) && arty %chin% dimnames(gwdata$data)[[5L]]){
      gwdata$data[,,, tspos, arty] <<-
        readBin(flcon, "double",
                prod(dis$extent[c("NCOL", "NROW", "NLAY")]), 4L)
      cat(substr(arty, 1L, 1L))
    }else{
      # parameter or time step not requested; skip
      readBin(flcon, "double",
              prod(dis$extent[c("NCOL", "NROW", "NLAY")]), 4L)
      cat(".")
    }
  }; cat("\n")
  
  invisible()
}


# input to MODPATH --------------------------------------------------------

#text is returned - needs to be written separately

#response file
arp <- "@RESPONSE\n"
rsptxt <- function(tlim, newcbf){
  #intro line not needed
  paste0(arp, c(paste0(dmrt, ".nam"), # name file giving model data
                if(tr) "2", #
                if(tr) "1 0",
                "Y",
                paste(tlim, "1"), # terminate simulation at this time
                if(tr) ifelse(newcbf, "1", "2"), # new cbf?
                if(tr) paste0(dmrt, ".cbf"), # name of cbf to be written or read
                "2", # pathline output
                "N", # don't calculate location at specific time points
                "1",
                paste0(dmrt, ".ptr"), # starting locations file
                rep("1", 2L),
                rep("N", 3L),
                "Y"), collapse = "\n")
}

#name file
namtxt <- paste(paste0("DIS    29    \'", mfrt, ".dis\'"),
                paste0("main      10      \'", dmrt , ".dat\'"),
                paste0("budget    17      \'", mfrt, ".cbb\'"),
                paste0("head(binary)  18      \'", mfrt, ".hds\'"), sep = "\n")

# ptr.dat columns x, y, zo, L and a single value for release time
# order is sorted within function, but ptr.dat must be given with the correct column names
# for writing the particle starting locations file
gd.PTR <- function(ptr.dat, rt){
  set(ptr.dat, NULL, c("C", "R", "It", "Jt", "Kt", "rt"), list(0L, 0L, 2L, 2L, 0L, rt))
  setcolorder(ptr.dat, c("C", "R", "L", "x", "y", "zo", "It", "Jt", "Kt", "rt"))
  
  ffmtptr <- mapply(formatC, ptr.dat,
                    width = c(4L, 4L, 3L, 17L, 17L, 17L, 2L, 2L, 2L, 13L),
                    digits = c(0L, 0L, 0L, 8L, 8L, 8L, 0L, 0L, 0L, 3L),
                    format = c("d", "d", "d", "e", "e", "e", "d", "d", "d", "f"))
  
  if(nrow(ptr.dat) == 1L) str_c(ffmtptr, collapse = "") else apply(ffmtptr, 1L, str_c, collapse = "")
}

#write the name file - only one needed because data file can be changed (and even then, particles are called separately)
write(namtxt, paste0(dmrt, ".nam"))

#main data file
dattxt <- function(por, MXP = 1000L, dis = paste0(mfrt, ".dis"), bas = paste0(mfrt, ".bas")){
  if(is.character(dis)) dis <- read.DIS(dis) #dis may be specified as ready-read list or file name to read
  if(is.character(bas)) bas <- read.BAS(bas, dis)
  
  txt <- character(7L)
  
  txt[1] <- paste0(FFf(2^35, 16L, 0L), FFe(999, 16L, 6L, 3L),
                   FFe(10^30, 16L, 6L, 3L), "  ", as.integer(MXP), "  1  1")
  txt[3] <- paste(dis$LAYCBD, collapse = " ")
  txt[4] <- paste(if(is.vector(ib <- bas$IBOUND)) vapply(ib, function(val) RIARRAY(CNSTNT = val, FMTIN_type = "i", FMTIN_w = 3L, flag.no = 10L), character(1)) else
    if(length(dim(ib)) == 2L) RIARRAY(arr = ib, FMTIN_type = "i", FMTIN_w = 3L, flag.no = 10L) else{
      apply(ib, 3L, RIARRAY, FMTIN_type = "i", FMTIN_w = 3L, flag.no = 10L)
    }, collapse = "\n")
  
  #porosity is given either as a single value, a vector (value per layer) or an array (value per cell)
  if(is.vector(por)){
    lpor <- double(dis$extent["NLAY"])
    lpor[] <- por #fills layer values - works for uniform or value per layer
    txt[5] <- paste(vapply(lpor, function(p) RIARRAY(CNSTNT = p), character(1L)), collapse = "\n")
  }else{
    txt[5] <- if(identical(length(dim(por)), 2L)) RIARRAY(arr = por, flag.no = 1) else{
      paste(apply(por, 3L, RIARRAY, flag.no = 1), collapse = "\n")
    }
  }
  
  #TBEGIN - doesn't affect the running of MODPATH; this is the reference time for the start of the MODFLOW model
  txt[6] <- " 0.0"
  
  #starting and ending timesteps to process - it may be that using only the necessary timesteps will give a good saving in MODPATH run time; if you wish to develop this, then a new DAT file will be required for each DMOC step
  txt[7] <- paste(c("", "1", "1", dis$extent["NPER"],
                    dis$sps[dis$extent["NPER"], "NSTP"]), collapse = "  ")
  
  return(paste(txt, collapse = "\n"))
}

MXP.def <- 50000L
if(write.dat || !file.exists(paste0(dmrt, ".dat")))
  write(dattxt(phi_e, MXP.def, dis), paste0(dmrt, ".dat"))



# simulation --------------------------------------------------------------

# time steps: ensured that they do not extend beyond the time period of the
#  MODFLOW model
if(start.t < MFt0) start.t <- MFt0
if(end.t > last(gwtime$t) + MFt0) end.t <- last(gwtime$t) + MFt0
tvals <- seq(start.t, end.t, Delta.t)

# ensure get to end even if duration is not multiple of Delta.t
if(last(tvals) != end.t) tvals <- c(tvals, end.t)

# number of time steps
nts <- length(tvals)

cat("simulation period is from", tvals[1L], "to", tail(tvals, 1L), "\n")

#initialise outflux data
fluxout <- vector("list", nts)

# mobile phase list pre-allocation
mob <- vector("list", nts)

# immobile phase list pre-allocation
if(sorb) immob <- vector("list", nts)

# release particles list pre-allocation
rel <- vector("list", nts)

# initialise lost mass vector and degraded mass vector
degraded <- double(nts)
massloss <- structure(replicate(6L, double(nts), FALSE),
                      names = c("bottom", "left", "top",
                                "right", "other", "inactive"))

# steady state simulations do not need a cbf file
newcbf <- ifelse(tr, newcbf, FALSE)

# if an initial condition is specified
if(exists("load.init") && load.init){
  cat("loading plume from existing dataset\n")
  res.init <- list.load(init.from)
  
  if(identical(rel.fun, "read")) rel.fun <- res.init$release.rates
  if(identical(xy0, "read")) xy0 <- res.init$release.loc[, 1:2]
  
  if(!any(res.init$time < start.t)){
    warning("specified initial state starts after start time for current simulation, so no starting plume is given")
  }else{
    #find best time step to use
    ts.init <- sum(res.init$time < start.t) # the number of time steps before current simulation start time
    
    #adjust start time and time points
    start.t <- res.init$time[ts.init]
    tvals <- unique(c(seq(start.t, tvals[1L], Delta.t), tvals))
    
    #read plume into initial conditions
    mob[[1L]] <- res.init$plume[ts == ts.init]
    mob[[1L]][, c("ts", "z") := list(1L, NULL)]
    if(sorb && "sorbed" %chin% names(res.init)){
      immob[[1L]] <- res.init$sorbed[ts == ts.init]
      immob[[1L]][, c("ts", "z") := list(1L, NULL)]
    }
    
    rm(res.init, ts.init)
  }
}

colord <- c("x", "y", "L", "zo", "m")
relstate0 <- data.table(x = xy0[, 1L], y = xy0[, 2L], L = as.integer(L), zo = zo)

# plot wells on go? currently disabled
# best way to re-enable would be to use read.WEL function
pw <- FALSE

#a rectangle representing the model bound
MFdx <- diff(range(gwdata$gccs))
MFdy <- diff(range(gwdata$grcs))
bbox.poly <- cbind(x = c(0, MFdx, MFdx, 0) + MFxy0[1L],
                   y = c(0, 0, MFdy, MFdy) + MFxy0[2L])

# initialise groundwater data arrays
gwdata$time <- 0
gwdata$data <- array(0, c(dis$extent[c("NCOL", "NROW", "NLAY")], 1L,
                          unname(3L + (dis$extent["NLAY"] > 1L) + tr)),
                     list(NULL, NULL, NULL, "0_0",
                          c("Head",
                            "FlowRightFace",
                            "FlowFrontFace",
                            if(dis$extent["NLAY"] > 1L) "FlowLowerFace",
                            if(tr) "Storage")))
wtop <- array(0, c(dis$extent[c("NCOL", "NROW", "NLAY")], 1L))
new.mfts <- 0L

for(tPt in 2:nts){
  st.time <- Sys.time()
  
  # load required MODFLOW data for this time step
  # - currently loaded time steps
  curr.mfts <- new.mfts
  # - time range required
  t0 <- tvals[tPt - 1L]; t1 <- tvals[tPt]
  # - MODFLOW time steps which are required to cover that time range
  mfts0 <- cellref.loc(t0, c(0, gwtime$t) + MFt0)
  mfts1 <- cellref.loc(t1, c(0, gwtime$t) + MFt0)
  #  -- case in which t1 is the same as the end of the MODFLOW model
  if(is.na(mfts1)) mfts1 <- nrow(gwtime)
  # - update time values describing the ends of these MODFLOW time steps
  gwdata$time <- gwtime[mfts0:mfts1, "t"]
  # - the set of MODFLOW time steps required
  new.mfts <- mfts0:mfts1
  
  # - update the arrays (if necessary)
  if(!identical(curr.mfts, new.mfts)){
    renew.gwdata(curr.mfts, new.mfts, dis$extent["NLAY"] > 1L, tr)
  }
  
  # start where left off
  
  state <- copy(mob[[tPt - 1L]])
  if(is.null(state)){
    state <- data.table(x = double(0L), y = double(0L), L = integer(0L), zo = double(0L), m = double(0L))
  }else state[, ts := NULL]
  
  if(sorb){
    statei <- copy(immob[[tPt - 1L]])
    if(is.null(statei)){
      statei <- data.table(x = double(0L), y = double(0L), L = integer(0L), zo = double(0L), m = double(0L))
    }else statei[, ts := NULL]
  }
  
  t.old <- tvals[tPt - 1L]; t.new <- tvals[tPt]; dt <- t.new - t.old
  
  # mass releases during this time step
  
  #add new release?
  relm <- vapply(rel.fun, function(fun){
    #mass released in time step found by integration of release rate functions
    tpts <- seq(t.old + dt/200, t.new - dt/200, length.out = 100L)
    vs <- vapply(tpts, fun, double(1L)) #fun need not be vectorised
    sum(vs)*dt/100
  }, double(1L))
  if(any(is.na(relm))) warning("source term function returning NA at timestep ", tPt)
  relstate <- if(any(relm > 0, na.rm = TRUE)) cbind(relstate0[relm > 0,], m = relm[relm > 0])
  
  state <- rbind(relstate, state); if(identical(nrow(state), 0L)) state <- NULL
  
  # save releases separately
  if(!is.null(relstate) && nrow(relstate) > 0L){
    rel[[tPt - 1L]] <- copy(relstate)
    rel[[tPt - 1L]][, ts := tPt - 1L]
    setcolorder(rel[[tPt - 1L]], c("ts", "x", "y", "L", "zo", "m"))
  }
  
  if(is.null(state)) next #this time is before the first release or complete clean-up has occurred
  
  #display current timestep details once it is established that something will be happening
  cat("timestep ", tPt, ", up to t = ", t.new, ", with ", nrow(state), " particles\n", sep = "")
  
  rls.time <- Sys.time()
  
  # propagation: execute the transport algorithm for the current time step
  # advection, sinks and reactions, dispersion
  
  OUTts <- prop(state, t.new, dt, newcbf, phi_e,
                if(sorb) statei, if(sorb) Rf, sorb,
                new.mfts)
  if(sorb){
    state <- OUTts[[1L]]; statei <- OUTts[[2L]]
  }else state <- OUTts
  rm(OUTts)
  newcbf <- F #only needed the first time
  state[, c("L", "zo") := list(ifelse(is.na(L), 0L, L),
                               ifelse(is.na(zo), NA, zo))]
  if(sorb) statei[, c("L", "zo") := list(ifelse(is.na(L), 0L, L),
                                         ifelse(is.na(zo), NA, zo))]
  
  prop.time <- Sys.time()
  
  # coalescence: clump nearby particles together for efficient representation of concentration field
  
  #it may be that all particles are abstracted - in which case there is no active mass at this timestep
  mob[[tPt]] <- if(nrow(state) > 50L){
    coalesce(state, cd, mm, tvals[tPt])
  }else state
  
  if(sorb){
    immob[[tPt]] <- if(nrow(statei) > 50L){
      coalesce(statei, cd, mm, tvals[tPt])
    }else statei
  }
  
  mob[[tPt]][, ts := tPt]
  setcolorder(mob[[tPt]], c("ts", colord)) #reorder
  if(sorb){
    immob[[tPt]][, ts := tPt]
    setcolorder(immob[[tPt]], c("ts", colord))
  }
  
  co.time <- Sys.time()
  
  # clean up of escaped particles
  
  #check which particles are in model bound
  #any which are not are deleted and their mass is saved in massloss (vector with one value per timestep)
  inmodxy <- point.in.polygon(mob[[tPt]]$x, mob[[tPt]]$y, bbox.poly[, "x"], bbox.poly[, "y"]) == 1L
  inmodz <- mob[[tPt]]$L >= 1L & mob[[tPt]]$L <= dis$extent["NLAY"] & !is.na(mob[[tPt]]$L)
  if(any(!(inmod <- inmodxy & inmodz))){
    # the following algorithm determines which dimension (2D) the mass has escaped out of
    outmod <- mob[[tPt]][!inmod]
    outmod[, c("bottom", "left", "top", "right", "other") := {
      modmid <- colMeans(bbox.poly)
      yovx <- MFdy/MFdx
      yintp <- modmid[2L] - yovx*modmid[1L]
      yintn <- modmid[2L] + yovx*modmid[1L]
      list(y <= MFxy0[2L] & y <= yintp + yovx*x & y < yintn - yovx*x,
           x <= MFxy0[1L] & y > yintp + yovx*x & y <= yintn - yovx*x,
           y >= MFxy0[2L] + MFdy & y >= yintp + yovx*x & y > yintn - yovx*x,
           x >= MFxy0[1L] + MFdx & y > yintp + yovx*x & y <= yintn - yovx*x,
           is.na(x) | is.na(y))
    }]
    outmod[, {
      massloss[["bottom"]][tPt] <<- sum(m[bottom], na.rm = TRUE)
      massloss[["left"]][tPt] <<- sum(m[left], na.rm = TRUE)
      massloss[["top"]][tPt] <<- sum(m[top], na.rm = TRUE)
      massloss[["right"]][tPt] <<- sum(m[right], na.rm = TRUE)
      massloss[["other"]][tPt] <<- sum(m[other])
    }]
    mob[[tPt]] <- mob[[tPt]][inmod,]
  }
  
  # live plotting if requested
  
  if(plot.on.go && tPt != 1L && !is.null(mob[[tPt - 1L]])){
    maxm <- max(mob[[tPt - 1L]]$m, rel[[tPt - 1L]]$m)
    mfts <- cellref.loc(tvals[tPt - 1L], c(0, gwdata$time) + MFt0)
    for(lay in sort(unique(c(mob[[tPt - 1L]]$L, rel[[tPt - 1L]]$L)))){
      #plot model active region, with constant heads shown in blue
      with(gwdata,
           MFimage(bas$IBOUND[,, lay],
                   gccs + MFxy0[1L], grcs + MFxy0[2L],
                   col = c("blue", "grey", "white"), zlim = c(-1, 1),
                   xlab = "easting", ylab = "northing"))
      # plot particles, with grey shade indicating mass
      setkey(mob[[tPt - 1L]], m)
      mob[[tPt - 1L]][L == lay,
                      points(x, y,
                             col = grey((1 - m/maxm)*.7),
                             cex = .3, pch = 16L)]
      if(!is.null(rel[[tPt - 1L]]))
        rel[[tPt - 1L]][L == lay, points(x, y, col = "darkred",
                                         cex = .3, pch = 16L)]
      
      #add title and indication of mass magnitude
      title(main = paste0("t = ", tvals[tPt - 1L], ", layer ", lay),
            sub = paste0("maximum mass = ", signif(maxm, 4L)))
      if(pw){
        #plot active wells in this layer
        with(gwdata, MFimage(data[,, lay, mfts, "Wells"] != 0,
                             gccs + MFxy0[1L], grcs + MFxy0[2L], 0:1, c("transparent", "red"),
                             add = TRUE))
        #plot active wells in other layers, semi-transparent
        if(any(lay != unique(L))){
          with(gwdata, MFimage(rowSums(data[,, -lay, mfts, "Wells", drop = FALSE], dims = 2L) != 0,
                               gccs + MFxy0[1L], grcs + MFxy0[2L], 0:1, c("transparent", "#FF000080"),
                               add = TRUE))
        }
      }
    }
  }
  
  # print execution time summary for the time step
  
  plot.time <- Sys.time()
  print(diff(c(st.time, release = rls.time, propagate = prop.time, coalesce = co.time, plot = plot.time)))
}

close(hdcon)
close(flcon)

# post-process ------------------------------------------------------------

rm(list = c("state", if(sorb) "statei"))
cat("Simulation complete. Now organising and post-processing results.\n")

cat("binding particle and outflux data into single data tables...\n")
mob <- rbindlist(mob)
setkey(mob, ts)

if(sorb) immob <- rbindlist(immob)
if(sorb) setkey(immob, ts)

rel <- rbindlist(rel)
setkey(rel, ts)

fluxout <- rbindlist(fluxout[sapply(fluxout, is.data.table)], use.names = TRUE)
setcolorder(fluxout, c("ts", "C", "R", "L", "J_out"))
setkey(fluxout, ts)

massloss <- do.call(cbind, massloss)

#determine z values
cat("z calculation...\n")
zexpr <- expression({
  C <- cellref.loc(x, gccs + MFxy0[1L], FALSE)
  R <- cellref.loc(y, grcs + MFxy0[2L], TRUE)
  mfts <- cellref.loc(ts, c(0, time) + MFt0, FALSE)
  bot <- elev[cbind(C, R, L + 1L)]
  thk <- wtop[cbind(C, R, L, mfts)] - bot
  bot + zo*thk
})

mob[, z := with(gwdata, eval(zexpr))]
if(sorb) immob[, z := with(gwdata, eval(zexpr))]
rel[, z := with(gwdata, eval(zexpr))]

setcolorder(mob, c("ts", letters[24:26], "L", "zo", "m"))
if(sorb) setcolorder(immob, c("ts", letters[24:26], "L", "zo", "m"))
setcolorder(rel, c("ts", letters[24:26], "L", "zo", "m"))

#kernel density estimate and plot
if(!ThreeDK) nkcell <- nkcell[1:2] # for safety
ksOUT <- array(0, dim = c(nkcell, ts = length(tvals)))
ksDATA <- NULL
Vkcell <- MFdx*MFdy*(if(ThreeDK) diff(Kzlim) else 1)/prod(nkcell[1:ifelse(ThreeDK, 3L, 2L)]) # ks cell volume
cat("Kernel Smooth: timestep      ")
mob[, {
  #perform the kernel smooth in 2 or 3 dimensions
  k <- kde(cbind(x, y, if(ThreeDK) z),
           H = diag(c(rep(smd[1L]^2, 2L), if(ThreeDK) smd[2L]^2), ifelse(ThreeDK, 3L, 2L)), # should be ^3 if ThreeDK?
           bgridsize = nkcell, binned = TRUE,
           xmin = c(MFxy0, if(ThreeDK) Kzlim[1L]), # minima in all dimensions
           xmax = c(MFxy0 + c(MFdx, MFdy), if(ThreeDK) Kzlim[2L]), # maxima in all dimensions
           w = m/mean(m)) # weights (kde insists that weights sum to the number of points, corrected later)
  
  # scale the results to represent concentration
  if(ThreeDK){
    ksOUT[,,, ts] <<- k$estimate*sum(m)
  }else{
    ksOUT[,, ts] <<- k$estimate*sum(m)
  }
  
  # save the evaluation points and the H matrix
  if(is.null(ksDATA)) ksDATA <<- k[c("eval.points", "H")]
  
  cat("\b\b\b\b\b", FFI(ts, 5L), sep = "")
}, by = ts]; cat("\n")

# remove unneeded time steps
ksDATA$time <- tvals
ksDATA$nkcell <- nkcell
ksDATA$Vkcell <- Vkcell

sim.end <- Sys.time()

cat("saving...\n")
list.save(list(plume = mob,
               sorbed = if(sorb) immob,
               release = rel,
               KSplume = list(k = ksOUT, info = ksDATA, "smooth" = smd[if(ThreeDK) 1:2 else 1L],
                              "number of divisions" = nkcell, "kcell volume or area" = Vkcell),
               fluxout = fluxout,
               degradedmass = degraded,
               lostmass = massloss,
               time = tvals,
               D = list("3Ddisp" = ThreeDD,
                        "D" = c(DL = DL, DT = DT, DV = if(ThreeDD) DV),
                        "vdepD" = vdepD,
                        "retain vertical loss" = if(ThreeDD) retain.vloss else NA),
               react = mget(c("sorb", if(sorb) "Rf", "lambda", "decaysorbed")),
               porosity = phi_e,
               release.loc = data.frame(xy0, L = L, zo = zo),
               release.rates = rel.fun,
               MFbounds = list(origin = c(MFxy0, t = MFt0),
                               bounds = bbox.poly),
               coalesce = mget(c("cd", "mm", "maxp")),
               description = description,
               timings = c(start = sim.start, end = sim.end)),
          file = paste0(dmrt, "_", info, ".rds"))

setwd(od)
cat("Execution complete.  Results saved to\n", mfdir, dmrt, "_", info, ".rds\n", sep = "")
