# analysis of GDPS results

source("C:/Users/cjb309/Dropbox/Scripts/R/plottools.R")
library("plyr"); library("data.table"); library("RNetCDF")

#' "Which" for NetCDF arrays
#' 
#' @description 
#' A version of "which" for NetCDF datasets that does not require prior
#' loading of arrays and avoids loading the whole of large arrays at once.
#'
#' @param ncfile NetCDF connection object
#' @param variable variable name within NetCDF data set
#' @param FUN function to operate on array, which should produce a logical
#' vector or array
#' @param arr.ind,useNames passed to which
#' @param size.threshold how large can the array be before the function
#' performs on chunks?
#' @param ... additional arguments to FUN
#'
#' @return
#' integer vector (arr.ind = FALSE) or matrix (arr.ind = TRUE)
#' @export
#'
#' @examples
which.nc <- function(ncfile, variable, FUN = Negate(is.na), ..., arr.ind = FALSE,
                     useNames = TRUE, size.threshold = 1e7){
  # get dimensions of array
  ards <- sapply(lapply(var.inq.nc(ncfile, variable)$dimids,
                        dim.inq.nc, ncfile = ncfile),
                 `[[`, "length")
  nards <- length(ards)
  arsize <- prod(ards)
  FUN <- match.fun(FUN)
  
  if(arsize > size.threshold){
    # large array, perform in chunks
    nchunk <- prod(ds) %/% size.threshold + 1L
    chsize <- ceiling(last(ds)/nchunk)
    
    # define chunk indices - roughly equal sections along highest index
    chunks <- lapply(1:nchunk, function(i){
      (1:last(ds))[(1:last(ds) %/% chsize + 1L) == i]
    })
    
    sections <- lapply(chunks, function(ch){
      x <- var.get.nc(ncfile, variable, c(rep(NA, nards - 1L), ch[1L]),
                      c(rep(NA, nards - 1L), last(ch) - ch[1L] + 1L),
                      collapse = FALSE)
      wh <- which(FUN(x, ...), arr.ind, useNames)
      if(arr.ind){
        wh[, nards] <- wh[, nards] + ch[1L] - 1L
      }else wh <- wh + prod(ards[-nards])*(ch[1L] - 1L)
      wh
    })
    
    if(arr.ind){
      do.call(rbind, sections)
    }else c(sections, recursive = TRUE)
  }else{
    # small array, perform all together
    x <- var.get.nc(ncfile, variable, collapse = FALSE)
    which(FUN(x, ...), arr.ind, useNames)
  }
}

#plots colour-flooded kernel-smoothed plume with optional particle overlay
#also plots the model boundaries using the bas file and the active wells at the time represented by the plot
#the plots are saved as png files
plot.plume <- function(res, gwdata, bas, folder, prefix, p = TRUE, sc = TRUE, tss = NULL, to.png = TRUE,
                       width = 7.5, height = 6, unit = "in", resolution = 144, mai = c(1, 1, .5, .5),
                       pcol = "#00000020", plot.pars = list(), breaks = 10^seq(-6, 0, .5),
                       time.to.date = TRUE, all.tss = FALSE, ...){
  if(is.null(tss) && !to.png && !all.tss) stop("really plot all time steps in R graphics device?\n",
                                               "use all.tss = TRUE if so")
  if(is.null(tss)) tss <- seq_along(res$time)
  
  with(as.list(res$MFbounds$origin), {
    for(l in ls()) assign(paste0("MF", l, "0"), eval(as.name(l)), parent.env(environment()))
  })
  
  l_ply(tss, function(tpt){
    if(to.png){
      png(paste0(folder, prefix, "_", gsub(" ", "0", FFI(tpt, 3L)), if(p) "p", ".png"),
          width, height, unit, res = resolution, ...)
      on.exit(dev.off())
      par(mai = mai)
    }
    
    ThreeDK <- length(res$KSplume$info$eval.points) == 3L
    with(res$KSplume, do.call(image, c(list(info$eval.points[[1L]], info$eval.points[[2L]],
                                            if(ThreeDK) rowMeans(k[,,, tpt], dims = 2L) else k[,, tpt],
                                            col = colfl(12L), breaks = breaks,
                                            xlab = "easting", ylab = "northing", asp = 1,
                                            main = if(time.to.date){
                                              paste("year", floor(res$time[tpt]/365.25 + 1900))
                                            }else paste("time =", signif(res$time[tpt], 3L))),
                                       plot.pars)))
    
    gccs <- var.get.nc(gwdata, "gccs")
    grcs <- var.get.nc(gwdata, "grcs")
    MFimage(bas$IBOUND[,, 1], gccs + MFx0, grcs + MFy0, c(-1, 1),
            c("blue", "grey", "transparent"), add = TRUE)
  
    mfts <- cellref.loc(res$time[tpt],
                        c(0, var.get.nc(gwdata, "time")) + MFt0)
    
    if("Wells" %chin% var.get.nc(gwdata, "parameters")){
      MFimage(var.get.nc(gwdata, "Wells",
                         c(1, 1, 1, mfts), c(NA, NA, NA, 1)),
              gccs, grcs, 0:1, c("transparent", "darkred"), add = TRUE)
    }
    
    if(p) res$plume[ts == tpt,
                    points(x, y, pch = 16L, cex = .4, col = pcol)]
    
    if(sc) points(unique(res$release.loc[, 1:2], MARGIN = 1L),
                  col = "purple", cex = 1.5, lwd = 2)
  })
}

#the following assumes that obsC is a list of data frames with columns for time and concentration
#obsC must be named with each element referring to a column and row: "C<column>R<row>" e.g. "C28R21"
#in this way the data is matched with the simulation
#obsC.unitmatch gives the scaling to get from the units of res to the units of obs
#e.g. if res is in kg/m^3 and obsC is in ug/l, then obsC.unitmatch should be 1e6, which is the default
#if multiple.res = TRUE, res should be a list of results lists which will be plotted on top of each other
plot.abstracted <- function(res, gwdata, obsC, year.range = NULL, lcol = "red",
                            obsC.unitmatch = 1e6, unit.label = "\u00b5g/l", multiple.res  = FALSE,
                            mr.labels = paste("simulation", 1:length(res)), legpos = "topright",
                            time.to.date = TRUE){
  WellCR <- which.nc(gwdata, "Wells", `!=`, 0, arr.ind = TRUE)
  rownames(wellCR) <- apply(wellCR, 1L, function(cr) paste0("C", cr[1L], "R", cr[2L]))
  
  MFt0 <- if(!multiple.res) res$MFbounds$origin["t"] else sapply(res, with, MFbounds$origin["t"])
  
  tadapt <- if(time.to.date) function(t) t/365.25 + 1900 else identity
  
  l_ply(names(obsC), function(w){
    obs <- obsC[[w]]
    cr <- wellCR[w,]
    
    Jexpr <- expression({
      Jtmp <- double(length(time))
      fluxout[C == cr[1] & R == cr[2], Jtmp[ts] <<- sum(J_out), by = ts]
      Jtmp
    })
    
    J <- if(!multiple.res) with(res, eval(Jexpr)) else lapply(res, function(r) with(r, eval(Jexpr)))
    
    if(any(c(J, recursive = TRUE) > 0) || any(obs$conc > 0)){
      Q <- if(!multiple.res){
        mfwell <- cbind(time = var.get.nc(gwdata, "time"),
                        Q = -colSums(var.get.nc(gwdata, "Wells",
                                                c(cr, NA, NA),
                                                c(1L, 1L, NA, NA))))
        with(mfwell, approx(time + MFt0, colSums(-data[cr[1], cr[2],,, "Wells"]),
                            res$time, "constant", f = 1)$y)
      }else lapply(1:length(res), function(n){
        mfwell <- cbind(time = var.get.nc(gwdata, "time"),
                        Q = -colSums(var.get.nc(gwdata, "Wells",
                                                c(cr, NA, NA),
                                                c(1L, 1L, NA, NA))))
        with(mfwell, approx(time + MFt0[n], colSums(-data[cr[1], cr[2],,, "Wells"]),
                            res$time, "constant", f = 1)$y)
      })
      
      conc <- if(!multiple.res) J/Q*obsC.unitmatch else Map(`/`, lapply(J, `*`, obsC.unitmatch), Q)
      
      ctmp <- c(conc, recursive = TRUE)
      ylm <- c(0, max(ctmp[is.finite(ctmp)], obs$conc))
      
      if(multiple.res) lcols <- lighten(rep(lcol, length(res)),
                                        2^seq(1, -1, length.out = length(res)))
      
      # with long time steps there may only be one concentration result within an abstraction period, which would result in a line not being plotted, because the concentration result would be surrounded by NA values
      # this is a solution, a slight cheat really - it replicates the isolated reading halfway between it and the previous time value
      # note that this does not work if the isolated reading is the first or last time step
      if(!multiple.res){
        if(any(isolated <- {
          dwnf <- diff(which(!is.finite(conc)))
          sapply(which(dwnf == 2L), function(i) sum(dwnf[1:(i - 1L)] - 1L) + i + 1L)
        })){
          t.to.plot <- {
            ttmp <- double(length(res$time) + length(isolated))
            itofill1 <- 1:length(res$time) + sapply(1:length(res$time), function(i) sum(i >= isolated))
            itofill2 <- which(!1:length(ttmp) %in% itofill1)
            ttmp[itofill1] <- res$time
            ttmp[itofill2] <- sapply(isolated, function(i) median(ttmp[itofill1[i] - c(2L, 0L)]))
            ttmp[itofill2 + 1L] <-
              sapply(isolated, function(i) quantile(ttmp[itofill1[i] - c(1L, -1L)], 2/3))
            ttmp
          }
          
          C.to.plot <- {
            ttmp <- double(length(conc) + length(isolated))
            itofill1 <- 1:length(conc) + sapply(1:length(conc), function(i) sum(i >= isolated))
            itofill2 <- which(!1:length(ttmp) %in% itofill1)
            ttmp[itofill1] <- conc
            ttmp[itofill2] <- ttmp[itofill1[isolated]]
            ttmp
          }
        }else{
          t.to.plot <- res$time
          C.to.plot <- conc
        }
      }else{
        if(any(c(isolated <- lapply(conc, function(cnc) {
          dwnf <- diff(which(!is.finite(cnc)))
          sapply(which(dwnf == 2L), function(i) sum(dwnf[1:(i - 1L)] - 1L) + i + 1L)
        }),
                 recursive = TRUE))){
          t.to.plot <- Map(function(res2, iso){
            ttmp <- double(length(res2$time) + length(iso))
            itofill1 <- 1:length(res2$time) + sapply(1:length(res2$time), function(i) sum(i >= iso))
            itofill2 <- which(!1:length(ttmp) %in% itofill1)
            ttmp[itofill1] <- res2$time
            ttmp[itofill2] <- sapply(iso, function(i) median(ttmp[itofill1[i] - c(2L, 0L)]))
            ttmp[itofill2 + 1L] <-
              sapply(iso, function(i) quantile(ttmp[itofill1[i] - c(1L, -1L)], 2/3))
            ttmp
          }, res, isolated)
          
          C.to.plot <- Map(function(C, iso){
            ttmp <- double(length(C) + length(iso))
            itofill1 <- 1:length(C) + sapply(1:length(C), function(i) sum(i >= iso))
            itofill2 <- which(!1:length(ttmp) %in% itofill1)
            ttmp[itofill1] <- C
            ttmp[itofill2] <- ttmp[itofill1[iso]]
            ttmp
          }, conc, isolated)
        }else{
          t.to.plot <- lapply(res, `[[`, "time")
          C.to.plot <- conc
        }
      }
      
      plot(tadapt(if(!multiple.res) t.to.plot else t.to.plot[[1L]]),
           if(!multiple.res) C.to.plot else C.to.plot[[1L]],
           type = "l", lwd = 2, col = lcol,
           xlim = year.range, ylim = ylm, xlab = if(time.to.date) "date" else "time",
           ylab = paste0("abstracted concentration (", unit.label, ")"),
           main = w)
      if(multiple.res){
        Map(function(n, lc){
          lines(tadapt(t.to.plot[[n]]), C.to.plot[[n]], lwd = 2, col = lc)
        }, 2:length(res), lcols[-1L])
        
        legend(legpos, legend = mr.labels, lwd = 2, col = lcols)
      }
      points(tadapt(obs$time), obs$conc, cex = 1.5, lwd = 2)
    }
  })
}

# plot the mass balance through time
plot.mass.balance <- function(res, time.to.date = TRUE, main = ""){
  # input from sources
  in.sc <- with(res, {
    mtmp <- double(length(time))
    release[, mtmp[ts + 1L] <<- sum(m), by = ts]
    cumsum(mtmp)
  })
  
  # output to sinks
  out.abs <- with(res, {
    mtmp <- double(length(time))
    dt <- diff(time)
    fluxout[, mtmp[ts] <<- sum(J_out)*dt[ts - 1L], by = ts]
    cumsum(mtmp)
  })
  
  # active mass, mobile
  act.pl <- with(res, {
    mtmp <- double(length(time))
    plume[, mtmp[ts] <<- sum(m), by = ts]
    mtmp
  })
  
  # active mass, immobile
  act.i <- with(res, if(res$react$sorb){
    mtmp <- double(length(time))
    if("sorbed" %in% ls()) sorbed[, mtmp[ts] <<- sum(m), by = ts]
    mtmp
  }else double(length(time)))
  
  # mass lost from model
  # - backward compatible for single vector
  out.lost <- cumsum(rowSums(as.matrix(res$lostmass)))
  
  # degraded mass
  out.dec <- cumsum(res$degradedmass)
  
  plmx <- max(in.sc, act.pl, act.i, out.abs, out.lost, out.dec)
  date <- if(time.to.date) res$time/365.25 + 1900 else res$time
  
  ylm <- c(0, max(in.sc, act.pl, act.i, out.abs, out.lost, out.dec, na.rm = TRUE))
  
  plot(date, in.sc, type = "l", col = "red", ylim = ylm,
       xlab = if(time.to.date) "date" else "time", ylab = "mass (kg)", main = main)
  lines(date, act.pl)
  lines(date, act.i, lty = 2)
  lines(date, out.abs, col = "green")
  lines(date, out.lost, col = "green", lty = 3)
  lines(date, out.dec, col = "green", lty = 2)
  lines(date, in.sc - out.abs - out.dec - out.lost, col = "blue")
  lines(date, act.pl + act.i, col = "blue", lwd = 2, lty = 2)
  
  legend("topleft",
         legend = c("sources", "sinks", "degraded", "lost", "mobile", "immobile",
                    "sources - sinks", "total active"),
         lty = c(1, 1, 2, 3, 1, 2, 1, 2), lwd = c(rep(1, 7), 2),
         col = c("red", rep("green", 3), rep("black", 2), rep("blue", 2)))
}
