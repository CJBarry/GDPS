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
source("//COLLES-11893/Users/cjb309/Dropbox/Scripts/R/coalesce.R")
source("//COLLES-11893/Users/cjb309/Documents/GitHub/GDPS/src/GDPSfun2.R")

od <- getwd(); setwd(mfdir)
dis <- read.DIS(paste0(mfrt, ".dis"))
if(plot.on.go || write.dat) bas <- read.BAS(paste0(mfrt, ".bas"), dis)
if(!exists("sorb")) sorb <- FALSE

#reorganise groundwater data if not done so (similar idea to CBF and FTL files for MODPATH and MT3D respectively)
if(!exists("lRAM")) lRAM <- FALSE
if(!exists("lRAM2")) lRAM2 <- FALSE


# perform checks before starting ------------------------------------------

# check that there is correct number of release rate functions
if(length(rel.fun) != nrow(xy0)) stop("The number of release functions in rel.fun (which should",
                                      "be a list of functions), should equal the number of release",
                                      "points.\n",
                                      "Currently:\nlength(rel.fun) == ", length(rel.fun), "\n",
                                      "nrow(xy0) == ", nrow(xy0))


# load groundwater data ---------------------------------------------------

if(!lRAM && !lRAM2 && (reload.gwdata || !exists("gwdata") || !exists("wtop"))){
  gwdata <- if(file.exists(fnm <- paste0(mfrt, ".rds")) && !fresh.mfdata) list.load(fnm) else GWdata.save(mfdir, mfrt, fnm)
  gwdata$data[,,,, -1L][is.na(gwdata$data[,,,, -1L])] <- 0 #MODPATHmass uses flows in adjacent cells, so there can't be NA values
  
  #transient water surface top - either cell top or Head (water table) at any point
  wtop <- with(gwdata, {
    lt <- structure(rep(c(elev[,, -dim(elev)[3]]), times = length(time)), dim = dim(data)[1:4])
    wt <- adrop(data[,,,, "Head", drop = F], c(F, F, F, F, T))
    ifelse(lt > wt, wt, lt)
  }); wtop[is.na(wtop)] <- 999 #can't be having NA values here
  
  #save on memory - only keep what is needed
  gwdata$data <- with(gwdata, data[,,,, dimnames(data)[[5]] %in% c("Head",
                                                                   "FlowRightFace",
                                                                   "FlowFrontFace",
                                                                   "FlowLowerFace",
                                                                   "Storage",
                                                                   "Wells"), drop = F])
}

# lRAM means there is limited RAM available, so the whole transient dataset cannot simply be loaded in.  Therefore only the necessary time steps for each GDPS time step are loaded
if(lRAM){
  # get the co-ordinate information
  dis <- read.DIS(paste0(mfrt, ".dis"))
  gwdata <- list(data = NULL,
                 grcs = {
                   rsp <- dis$DELC
                   if(identical(names(rsp), "CNSTNT")) rsp <- rep(rsp, dis$extent["NROW"])
                   cumsum(c(0, rsp))
                 },
                 gccs = {
                   csp <- dis$DELR
                   if(identical(names(csp), "CNSTNT")) csp <- rep(csp, dis$extent["NCOL"])
                   cumsum(c(0, csp))
                 },
                 elev = dis$elev,
                 time = readHDS.arr(paste0(mfrt, ".hds"), time.only = TRUE, show.help = FALSE))
  mftime <- gwdata$time
  sp_ts.ts <- t(sapply(str_split(names(mftime), "_"), identity))
  mode(sp_ts.ts) <- "integer"
  
  # expression for loading the required timesteps each GDPS timestep
  datatimeexpr <- expression({
    hds <- readHDS.arr(paste0(mfrt, ".hds"), sp_ts = nects)
    cbb <- readCBB.arr(paste0(mfrt, ".cbb"), sp_ts = nects)
    ds <- dim(cbb) + c(rep(0L, 4L), 1L)
    dnms <- c(dimnames(cbb)[1:4],
              list(c("Head", dimnames(cbb)[[5L]])))
    list(data = array(c(hds$Head, cbb), ds, dnms),
         time = hds$time)
  })
  
  wtopexpr <- expression(with(gwdata, {
    wtop <- array(dim = dim(data)[1:4])
    for(tdim in 1L:dim(data)[4L]){
      wtab <- adrop(data[,,, tdim, "Head", drop = FALSE], c(F, F, F, T, T))
      wtop[,,, tdim] <- ifelse(wtab > elev, wtab, elev)
    }
    wtop
  }))
}

if(lRAM2 && reload.gwdata || !exists("gwdata") || !exists("wtop")){
  COLs <- Crange[1L]:Crange[2L]
  ROWs <- Rrange[1L]:Rrange[2L]
  
  dis <- read.DIS(paste0(mfrt, ".dis"))
  
  lr2gwdexpr <- expression({
    gwdata <- list(data = NULL,
                   grcs = {
                     rsp <- dis$DELC
                     if(identical(names(rsp), "CNSTNT")) rsp <- unname(rep(rsp, dis$extent["NROW"]))
                     cumsum(c(0, rsp))[dis$extent["NROW"] - rev(c(ROWs, last(ROWs) + 1L))]
                   },
                   gccs = {
                     csp <- dis$DELR
                     if(identical(names(csp), "CNSTNT")) csp <- unname(rep(csp, dis$extent["NCOL"]))
                     cumsum(c(0, csp))[c(COLs, last(COLs) + 1L)]
                   },
                   elev = dis$elev[COLs, ROWs,, drop = FALSE],
                   time = NULL)
    
    hds <- readHDS.arr(paste0(mfrt, ".hds"), CRs = list(COLs, ROWs))
    
    # doesn't matter if not all of artys are actually found in the file
    cbb <- readCBB.arr(paste0(mfrt, ".cbb"), hds = hds, CRLs = list(COLs, ROWs, "all"),
                       artys = c("Storage", "FlowRightFace", "FlowFrontFace", "FlowLowerFace", "Wells"))
    
    gwdata$time <- hds$time
    gwdata$data <- array(c(hds$Head, cbb), dim = dim(cbb) + c(0L, 0L, 0L, 0L, 1L),
                         dimnames = c(dimnames(cbb)[1:4], list(c("Head", dimnames(cbb)[[5L]]))))
    
    rm(cbb, hds)
    
    gwdata
  })
  
  if(fresh.mfdata || !file.exists(paste0(mfrt, lRAM2suff, ".rds"))){
    gwdata <- eval(lr2gwdexpr)
    list.save(gwdata, paste0(mfrt, lRAM2suff, ".rds"))
  }else gwdata <- list.load(paste0(mfrt, lRAM2suff, ".rds"))
  
  #check that the correct columns and rows have been read - stored in dimnames
  #if not, then re-read gwdata from model results and save as new list
  # note that "2" == 2L returns TRUE
  if(`||`(length(COLs) != dim(gwdata$data)[1L] || length(ROWs) != dim(gwdata$data)[2L],
     any(COLs != dimnames(gwdata$data)[[1L]]) || any(ROWs != dimnames(gwdata$data)[[2L]]))){
    gwdata <- eval(lr2gwdexpr)
    list.save(gwdata, paste0(mfrt, lRAM2suff, ".rds"))
  }

  wtop <- with(gwdata, {
    lt <- structure(rep(c(elev[,, -dim(elev)[3]]), times = length(time)), dim = dim(data)[1:4])
    wt <- adrop(data[,,,, "Head", drop = F], c(F, F, F, F, T))
    ifelse(lt > wt, wt, lt)
  }); wtop[is.na(wtop)] <- 999 #can't be having NA values here
}

#a rectangle representing the model bound
MFdx <- diff(range(gwdata$gccs)); MFdy <- diff(range(gwdata$grcs))
bbox.poly <- cbind(x = c(0, MFdx, MFdx, 0) + MFxy0[1L], y = c(0, 0, MFdy, MFdy) + MFxy0[2L])

#find HDRY
if(file.exists(paste0(mfdir, mfrt, ".lpf"))){
  HDRY <- scan(paste0(mfdir, mfrt, ".lpf"), list(integer(), double(), integer()), 1L, comment.char = "#")[[2L]]
}else if(file.exists(paste0(mfdir, mfrt, ".bcf"))){
  HDRY <- scan(paste0(mfdir, mfrt, ".bcf"), list(integer(), double(), integer()), 1L, comment.char = "#")[[2L]]
}else{
  if(interactive()){
    HDRY <- eval(parse(text = (readline("no value for HDRY found from LPF or BCF files; give the expected value for dry cells.  Note that unidentified dry cells can lead to an infinite loop.  Put value: "))))
  }else warning("No value for HDRY found.  This value is normally found as the second item of the BCF or LPF package files.  Make the appropriate file available or else write a line in the input script: \"HDRY <- ...\" to define.")
}
if(exists("HDRY") && !lRAM) wtop[abs(wtop) > abs(HDRY)*.99 & abs(wtop) < abs(HDRY)*1.01] <- NA
#dry cells might as well be no-flow for this algorithm
#approximate matching because very big HDRY values cause problems with double precision for exact matching
#assumed that HDRY is well outside range of proper head values


# analytical dispersion parameters ----------------------------------------

#number of dispersion steps per advection step and ready-solved imprint
#this has been solved using the analytical diffusion equation and integrating over the areas (or volumes) occupied by each new particle's region; each new particle is placed at the centre of mass of its region
#the algorithm finds the grid for which these numbers are true, which seems a bit backward, but that is the novelty of this code!
#refer to a2.p340
if(ThreeDD){
  #determined by numerical integration
  # rc <- 0.96645
  centre <- 0.4
  edge <- 0.1
  dispC <- 1.187558
}else{
  #analytically determined
  rc <- sqrt(log(2))
  centre <- 1 - exp(-rc^2) # 0.5
  edge <- exp(-rc^2)/4 # (1 - centre)/4 = 0.125
  dispC <- sqrt(8)/(pi*exp(-rc^2))*(rc*exp(-rc^2) + (sqrt(pi)/2)*(pracma::erfc(rc)))
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
  
  txt[1] <- paste0(FFf(2^35, 16, 0), FFe(999, 16, 6, 3), FFe(10^30, 16, 6, 3), "  ", as.integer(MXP), "  1  1")
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
  txt[7] <- paste(c("", "1", "1", dis$extent["NPER"], dis$sps[dis$extent["NPER"], "NSTP"]), collapse = "  ")
  
  return(paste(txt, collapse = "\n"))
}

MXP.def <- 50000L
if(write.dat || !file.exists(paste0(dmrt, ".dat"))) write(dattxt(phi_e, MXP.def, dis), paste0(dmrt, ".dat"))



# simulation --------------------------------------------------------------

#time steps: ensured that they do not extend beyond the time period of the MODFLOW model
tvals <- seq(ifelse(start.t >= MFt0, start.t, MFt0),
             ifelse(end.t < tail(gwdata$time, 1L) + MFt0, end.t, tail(gwdata$time, 1L) + MFt0), Delta.t)

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
degraded <- massloss <- double(length(tvals))

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

#plot wells on go? not an option with lRAM = TRUE
if(plot.on.go && !lRAM) pw <- "Wells" %in% dimnames(gwdata$data)[[5]] else pw <- FALSE
if(plot.on.go && lRAM2){
  truegw <- list(gccs = {
    csp <- dis$DELR
    if(identical(names(csp), "CNSTNT")) csp <- unname(rep(csp, dis$extent["NCOL"]))
    cumsum(c(0, csp))
  }, grcs = {
    rsp <- dis$DELC
    if(identical(names(rsp), "CNSTNT")) rsp <- unname(rep(rsp, dis$extent["NROW"]))
    cumsum(c(0, rsp))
  })
}

for(tPt in 2:nts){
  st.time <- Sys.time()
  
  # if using low RAM option one, load the required MODFLOW time steps into gwdata
  
  if(lRAM){
    t0 <- tvals[tPt]; t1 <- tvals[tPt + 1L]
    mfts0 <- cellref.loc(t0, c(0, mftime) + MFt0)
    mfts1 <- cellref.loc(t1, c(0, mftime) + MFt0)
    nects <- sp_ts.ts[mfts0:mfts1,]
    gwdata[c("data", "time")] <- eval(datatimeexpr)
    
    wtop <- eval(wtopexpr)
    if(exists("HDRY")) wtop[abs(wtop) > abs(HDRY)*.99 & abs(wtop) < abs(HDRY)*1.01] <- NA
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
  
  OUTts <- prop(state, t.new, dt, newcbf, phi_e, if(sorb) statei, if(sorb) Rf, sorb)
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
    massloss[tPt] <- sum(mob[[tPt]]$m[!inmod])
    mob[[tPt]] <- mob[[tPt]][inmod,]
  }
  
  # live plotting if requested
  
  if(plot.on.go && tPt != 1L && !is.null(mob[[tPt - 1L]])){
    maxm <- max(mob[[tPt - 1L]]$m, rel[[tPt - 1L]]$m)
    mfts <- cellref.loc(tvals[tPt - 1L], c(0, gwdata$time) + MFt0)
    for(lay in sort(unique(c(mob[[tPt - 1L]]$L, rel[[tPt - 1L]]$L)))){
      #plot model active region, with constant heads shown in blue
      with(if(lRAM2) truegw else gwdata,
           MFimage(bas$IBOUND[,, lay],
                   gccs + MFxy0[1L], grcs + MFxy0[2L],
                   col = c("blue", "grey", "white"), zlim = c(-1, 1),
                   xlab = "easting", ylab = "northing"))
      #plot particles, with opacity indicating mass
      mob[[tPt - 1L]][L == lay, points(x, y, col = rgb(.63, .13, .94, m[L == lay]/maxm), pch = 16L)]
      if(!is.null(rel[[tpt - 1L]]))
        rel[[tPt - 1L]][L == lay, points(x, y, col = rgb(.63, .13, .94, m[L == lay]/maxm), pch = 16L)]
      
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
           gridsize = nkcell[1:ifelse(ThreeDK, 3L, 2L)],
           xmin = with(gwdata, c(c(min(gccs), min(grcs)) + MFxy0, if(ThreeDK) Kzlim[1L])), # minima in all dimensions
           xmax = with(gwdata, c(c(max(gccs), max(grcs)) + MFxy0, if(ThreeDK) Kzlim[2L])), # maxima in all dimensions
           w = m/mean(m)) # weights (kde insists that weights sum to the number of points, so normalised concentration is returned; the result is corrected later to convert to true concentration)
  
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

if(save.res) cat("saving...\n")
if(save.res) list.save(list(plume = mob,
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
                                     "retain vertical loss" = if(ThreeDD) retain.vloss else NA),
                            react = mget(c("sorb", "Rf", "lambda", "decaysorbed")),
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
