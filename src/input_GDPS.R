#Christopher Barry, started on 23/02/2016 at University of Birmingham

# uses a simple MODFLOW model (dmoc2f) to develop features of working with groundwater model data
# this development stage incorporates MODPATHmass to account for mass loss to sinks and to save mass flux to receptors
# the input mass has also been adjusted so as to ensure correct total mass input when the timesteps do not match up with the release period well

#input----

#check:
# if you have changed the MODFLOW outputs, you must put fresh.mfdata <- TRUE, or an old version will be used
# if you have changed the MODFLOW model or the chosen porosity, you must put write.dat <- TRUE

#directories: R scripts and groundwater model working directories
gendir <- "C:/Users/cjb309/Documents/GitHub/tools/"
mfdir <- "" #the directory in which the MODFLOW input and output files are stored
mpexe <- "Mpathr5_0.exe" #full path and directory for the MODPATH version 5 executable
source(paste0(gendir, "td.R")) # do not change
# use the td(day, month, year) function for a consistent date-based time system in which days are whole numbers
# td(31, 12, 1899) = 1, so that it is consistent with the MS Excel date numbering system after 01/03/1900
# td will produce consistent negative results for any date before 30/12/1899 as well

#root name for MODFLOW files, and is the simulation transient?  freshly reformat MODFLOW data?
mfrt <- "" #root names for MODFLOW files (without extension)
tr <- TRUE #is the model transient
fresh.mfdata <- FALSE
reload.gwdata <- FALSE
newcbf <- FALSE

#root name for dmoc files and any extra info to append to saved file name, as well as a description string to attach to the list
dmrt <- "" # root name for the GDPS outputs
save.res <- TRUE # save the results or merely keep in RAM (almost always want TRUE)
info <- "" # a short string to append to file name to avoid overwriting (e.g. "run1")
description <- paste0("") # as much text as you want to describe the simulation - this information will be saved with the result

#source description: x, y, layer, z offset and flux as function of time
#xy0, L and zo must be vectors whose lengths are all equal, the number of distinct release points in 3D space
#rel.fun must be a list of functions of one variable (time) whose length is as above
xy0 <- cbind(x = c(),
             y = c())
L <- as.integer(c())
zo <- c()
rel.fun <- list(function() ...)

#dispersivity; if vdepD the DL and DT represent aL and aT (dispersion coefficients rather that dispersivities)
#DV is vertical dispersion coefficient or dispersivity if ThreeDD is TRUE
#retain.vloss states whether loss of mass out of top or bottom of model should be disallowed (it will be mirrored back)
DL <- 10
DT <- 1
vdepD <- TRUE
if(ThreeDD <- TRUE){
  DV <- .1
  retain.vloss <- TRUE
}

#reactions
if(sorb <- TRUE) Rf <- 2
lambda <- 0
decaysorbed <- FALSE

#timestep size
Delta.t <- 100

#simulation length
#use the td(d, m, y) function for consistent relation between day values and dates
start.t <- td(1, 1, 1925)
end.t <- td(1, 1, 2100)

#coalescing distance (horiz and vert); minimum mass to keep particle (otherwise merged with nearest non-small particle)
#how many time steps for each coalescence? (avoids needing to make small cd for small time step)
#and what is the maximum number of particles permissible at the end of a time step?
cd <- c(h = 10, v = 5)
mm <- 0.00001
maxp <- 25000L #R tip: use the L qualifier to specify that a number should stored as an integer

#MODPATH model set-up; write MODPATH dat input file afresh?
phi_e <- 0.1 #effective porosity
write.dat <- TRUE

#at what time value does the MODFLOW model start, and where is its origin?
MFt0 <- 0 # often a good idea to use the td function here
MFxy0 <- c(x = 0, y = 0)

#plot while evaluating?
plot.on.go <- TRUE

#kernel smooth paramters: smoothing distance and number of cells in each dimension
ThreeDK <- FALSE
smd <- c(h = 15, v = 2)
nkcell <- c(100L, 100L)
Kzlim <- c(0, 10)

source("GDPS.R")
