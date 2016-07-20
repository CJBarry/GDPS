#script file that executes the grid-decoupled plume simulation

#----

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
source("C:/Users/cjb309/Documents/GitHub/GDPS/src/coalesce.R")
source("C:/Users/cjb309/Documents/GitHub/GDPS/src/GDPSfun.R")

od <- getwd(); setwd(mfdir)
dis <- read.DIS(paste0(mfrt, ".dis"))
if(plot.on.go || write.dat) bas <- read.BAS(paste0(mfrt, ".bas"), dis)

#reorganise groundwater data if not done so (similar idea to CBF and FTL files for MODPATH and MT3D respectively)
if(reload.gwdata || !exists("gwdata") || !exists("wtop")){
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
    HDRY <- as.double(readline("no value for HDRY found from LPF or BCF files; give the expected value for dry cells.  Note that unidentified dry cells can lead to an infinite loop.  Put value: "))
  }else{warning("No value for HDRY found.  This value is normally found as the second item of the BCF or LPF package files.  Make the appropriate file available or else write a line in the input script: \"HDRY <- ...\" to define.")}
}
wtop[abs(wtop) > abs(HDRY)*.99 & abs(wtop) < abs(HDRY)*1.01] <- NA
#dry cells might as well be no-flow for this algorithm
#approximate matching because very big HDRY values cause problems with double precision for exact matching
#assumed that HDRY is well outside range of proper head values

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

#functions for writing MODPATH input files----
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

#PTR starting locations file
ptr.line <- function(ptr.dat1) paste(mapply(FFgen, c(0L, 0L, ptr.dat1[c(4, 1:3)], 2L, 2L, 0L, ptr.dat1[5L]),
                                            c("i", "i", "i", "e", "e", "e", "i", "i", "i", "f"),
                                            c(4L, 4L, 3L, 17L, 17L, 17L, 2L, 2L, 2L, 9L),
                                            c(NA, NA, NA, 8L, 8L, 8L, NA, NA, NA, 2L)), collapse = "")

# ptr.dat in order x, y, zo, L, rt
dmoc.PTR <- function(ptr.dat){
  ptr.dat <- as.matrix(ptr.dat)
  
  lns <- apply(ptr.dat, 1L, ptr.line)
  
  return(str_c(lns, collapse = "\n"))
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

#simulation----

#time steps: ensured that they do not extend beyond the time period of the MODFLOW model
tvals <- seq(ifelse(start.t >= MFt0, start.t, MFt0),
             ifelse(end.t < tail(gwdata$time, 1L) + MFt0, end.t, tail(gwdata$time, 1L) + MFt0), Delta.t)
if(last(tvals) == end.t) tvals <- c(tvals, end.t) # ensure get to end even if duration is not multiple of Delta.t
cat("simulation period is from", tvals[1L], "to", tail(tvals, 1L), "\n")

#initialise outflux data
# massout <- array(0, with(dis, c(extent[c("NCOL", "NROW", "NLAY")], nts = length(tvals))))
fluxout <- rep(data.table(ts = integer(0L),
                          C = integer(0L),
                          R = integer(0L),
                          L = integer(0L),
                          J_out = double(0L)), length(tvals))

#initialise lost mass vector
massloss <- double(length(tvals))

colord <- c("x", "y", "L", "zo", "m")
relstate0 <- data.table(x = xy0[, 1L], y = xy0[, 2L], L = as.integer(L), zo = zo)

newcbf <- ifelse(tr, newcbf, FALSE) #ss simulations do not need a cbf file

# mobile phase list pre-allocation
mob <- rep(list(data.table(ts = integer(0L),
                           x = double(0L),
                           y = double(0L),
                           L = integer(0L),
                           zo = double(0L),
                           m = double(0L))), length(tvals))

# immobile phase list pre-allocation
if(sorb) immob <- rep(list(data.table(ts = integer(0L),
                                      x = double(0L),
                                      y = double(0L),
                                      L = integer(0L),
                                      zo = double(0L),
                                      m = double(0L))), length(tvals))

# if an initial condition is specified
if(load.init){
  res.init <- list.load(init.from)
  if(!any(res.init$time < start.t)){
    warning("specified initial state starts after start time for current simulation, so no starting plume is given")
  }else{
    #find best time step to use
    ts.init <- sum(res.init$time < start.t) # the number of time steps before current simulation start time
    
    #adjust start time and time points
    start.t <- res.init$time[ts.init]
    tvals <- unique(c(seq(start.t, tvals[1L], Delta.t), tvals))
    
    #read plume into initial conditions
    mob[[1L]] <- res.init$plume[ts == ts.init]; mob[[1L]][, z := NULL]
    if(sorb && "sorbed" %chin% names(res.init)){
      immob[[1L]] <- res.init$sorbed[ts == ts.init]; immob[[1L]][, z := NULL]
    }
    
    rm(res.init, ts.init)
  }
}

for(tPt in 2:length(tvals)){
  st.time <- Sys.time()
  
  state <- copy(mob[[tPt - 1L]]); state[, ts := NULL]
  if(sorb){statei <- copy(immob[[tPt - 1L]]); statei[, ts := NULL]}
  t.old <- tvals[tPt - 1L]; t.new <- tvals[tPt]; dt <- t.new - t.old
  
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
  if(nrow(relstate) > 0L) mob[[tPt - 1L]] <- rbind(cbind(ts = tPt - 1L, relstate), mob[[tPt - 1L]])
  
  if(is.null(state)) next #this time is before the first release
  
  rls.time <- Sys.time()
  
  cat("timestep ", tPt, ", up to t = ", t.new, ", with ", nrow(state), " particles\n", sep = "")
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
  
  #check which particles are in model bound
  #any which are not are deleted and their mass is saved in massloss (vector with one value per timestep)
  inmodxy <- point.in.polygon(mob[[tPt]]$x, mob[[tPt]]$y, bbox.poly[, "x"], bbox.poly[, "y"]) == 1L
  inmodz <- mob[[tPt]]$L >= 1L & mob[[tPt]]$L <= dis$extent["NLAY"] & !is.na(mob[[tPt]]$L)
  if(any(!(inmod <- inmodxy & inmodz))){
    massloss[tPt] <- sum(mob[[tPt]]$m[!inmod])
    mob[[tPt]] <- mob[[tPt]][inmod,]
  }
  
  if(plot.on.go && tPt != 1L && !is.null(mob[[tPt - 1L]])){
    mob[[tPt - 1L]][, {
      maxm <- max(m)
      for(lay in sort(unique(L))){
        with(gwdata, MFimage(bas$IBOUND[,, lay], gccs + MFxy0[1L], grcs + MFxy0[2L], col = c("blue", "grey", "white"), zlim = c(-1, 1), show.range = F, xlab = "easting", ylab = "northing"))
        points(x[L == lay], y[L == lay], col = rgb(1, 0, 0, m[L == lay]/maxm), pch = 16L)
        title(main = paste0("t = ", tvals[tPt - 1L], ", layer ", lay), sub = paste0("maximum mass = ", signif(maxm, 4L)))
      }
    }]
  }
  
  plot.time <- Sys.time()
  print(diff(c(st.time, release = rls.time, propagate = prop.time, coalesce = co.time, plot = plot.time)))
}

#post-process----
rm(state, statei)
cat("Simulation complete. Now organising and post-processing results.\n")

cat("binding particle and outflux data into single data tables...\n")
mob <- rbindlist(mob)
setkey(mob, ts)
if(sorb){
  immob <- rbindlist(immob)
  setkey(immob, ts)
}
fluxout <- rbindlist(fluxout[sapply(fluxout, is.data.table)])
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
setcolorder(mob, c("ts", letters[24:26], "L", "zo", "m"))
if(sorb) setcolorder(immob, c("ts", letters[24:26], "L", "zo", "m"))

#kernel density estimate and plot
if(!ThreeDK) nkcell <- nkcell[1:2] # for safety
ksOUT <- array(NA_real_, dim = c(nkcell, ts = length(tvals)))
ksDATA <- NULL
Vkcell <- MFdx*MFdy*(if(ThreeDK) diff(Kzlim) else 1)/prod(nkcell) # ks cell volume
cat("Kernel Smooth: timestep      ")
mob[, {
  #perform the kernel smooth in 2 or 3 dimensions
  k <- kde(cbind(x, y, if(ThreeDK) z),
           H = diag(c(rep(smd[1L]^2, 2L), if(ThreeDK) smd[2L]^2), ifelse(ThreeDK, 3L, 2L)), # should be ^3 if ThreeDK?
           gridsize = nkcell,
           xmin = c(MFxy0, if(ThreeDK) Kzlim[1L]), # minima in all dimensions
           xmax = c(MFxy0 + c(MFdx, MFdy), if(ThreeDK) Kzlim[2L]), # maxima in all dimensions
           w = m/mean(m)) # weights (kde insists that weights sum to the number of points, corrected later)
  
  # scale the results to represent concentration
  if(ThreeDK){
    ksOUT[,,, ts] <<- k$estimate*sum(m)/Vkcell
  }else{
    ksOUT[,, ts] <<- k$estimate*sum(m)/Vkcell
  }
  
  # save the evaluation points and the H matrix
  if(is.null(ksDATA)) ksDATA <<- k[c("eval.points", "H")]
  
  cat("\b\b\b\b\b", FFI(ts, 5L), sep = "")
}, by = ts]; cat("\n")

# remove unneeded time steps
natss <- apply(ksOUT, ifelse(ThreeDK, 4L, 3L), function(ts) all(is.na(ts)))

ksOUT <- do.call(`[`, c(list(ksOUT), rep(list(bquote()), ifelse(ThreeDK, 3L, 2L)), list(!natss, drop = F)))
ksDATA$time <- tvals[!natss]

if(save.res) cat("saving...\n")
if(save.res) list.save(list(plume = mob,
                            sorbed = if(sorb) immob,
                            KSplume = list(ksOUT, ksDATA, "smooth" = smd,
                                           "number of divisions" = nkcell),
                            lostmass = massloss,
                            fluxout = fluxout,
                            time = tvals,
                            D = list("3Ddisp" = ThreeDD,
                                     "D" = c(DL = DL, DT = DT, DV = if(ThreeDD) DV),
                                     "retain vertical loss" = if(ThreeDD) retain.vloss else NA),
                            release.loc = xy0,
                            release.rates = rel.fun,
                            MFbounds = list(origin = c(MFxy0, t = MFt0),
                                            bounds = bbox.poly),
                            description = description),
                       file = paste0(dmrt, "_", info, ".rds"))

cat("plotting kernel-smoothed plume...\n")
# contour levels to plot: essentially ... 0.01, 0.0316, 0.1, 0.316, 1 ... (whole and half powers of 10)
plevs <- signif(10^(ceiling(log10(max(mob$m)/Vkcell)))*10^((-12:0)/2), 3L)
xlm <- unname({xrg <- quantile(mob$x, c(.01, .99)); dxrg <- diff(xrg); xrg + dxrg*c(-.2, .2)})
ylm <- unname({yrg <- quantile(mob$y, c(.01, .99)); dyrg <- diff(yrg); yrg + dyrg*c(-.2, .2)})

# limit to 20 contour plots
to.plot <- seq(1L, sum(!natss), ceiling(sum(!natss)/20))

l_ply(which(!natss)[to.plot], function(ts){
  mtx <- if(ThreeDK) rowMeans(ksOUT[,,, ts, drop = F], dims = 2L) else ksOUT[,, ts, drop = T]
  contour(ksDATA$eval.points[[1L]], ksDATA$eval.points[[2L]], mtx, levels = plevs,
          main = paste("t =", tvals[ts]), xlab = "x", ylab = "y", xlim = xlm, ylim = ylm)
})

setwd(od)
cat("Execution complete.  Results saved to\n", mfdir, dmrt, "_", info, ".rds\n", sep = "")
