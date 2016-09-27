# analysis of GDPS results

source("C:/Users/cjb309/Dropbox/Scripts/R/plottools.R")
library("plyr")

#plots colour-flooded kernel-smoothed plume with optional particle overlay
#also plots the model boundaries using the bas file and the active wells at the time represented by the plot
#the plots are saved as png files
plot.plume <- function(res, gwdata, bas, folder, prefix, p = TRUE, sc = TRUE, tss = NULL,
                       width = 7.5, height = 6, unit = "in", resolution = 144, mai = c(1, 1, .5, .5),
                       pcol = "#00000020", plot.pars = list(), ...){
  if(is.null(tss)) tss <- seq_along(res$time)
  
  with(as.list(res$MFbounds$origin), {
    for(l in ls()) assign(paste0("MF", l, "0"), eval(as.name(l)), parent.env(environment()))
  })
  
  l_ply(tss, function(tpt){
    png(paste0(folder, prefix, "_", gsub(" ", "0", FFI(tpt, 3L)), if(p) "p", ".png"),
        width, height, unit, res = resolution, ...)
    on.exit(dev.off())
    par(mai = mai)
    
    with(res$KSplume, do.call(image, c(list(info$eval.points[[1L]], info$eval.points[[2L]],
                                            rowMeans(k[,,, tpt], dims = 2L),
                                            col = colfl(12L), breaks = 10^seq(-8, -2, .5),
                                            xlab = "easting", ylab = "northing", asp = 1,
                                            main = paste("year", floor(res$time[tpt]/365.25 + 1900))),
                                       plot.pars)))
    
    with(gwdata, MFimage(bas$IBOUND[,, 1], gccs + MFx0, grcs + MFy0, c(-1, 1),
                         c("blue", "grey", "transparent"), add = TRUE))
    
    mfts <- cellref.loc(res$time[tpt], c(0, gwdata$time) + MFt0)
    
    with(gwdata, MFimage(rowSums(data[,,, mfts, "Wells"], dims = 2L) != 0,
                         gccs + MFx0, grcs + MFy0, 0:1, c("transparent", "darkred"), add = TRUE))
    
    if(p) res$plume[ts == tpt, points(x, y, pch = 16L, cex = .4, col = pcol)]
    
    if(sc) points(unique(res$release.loc[, c("x", "y")], MARGIN = 1L), col = "purple", cex = 1.5, lwd = 2)
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
                            mr.labels = paste("simulation", 1:length(res)), legpos = "topright"){
  wellCR <- which(rowSums(gwdata$data[,,,, "Wells"], dims = 2L) != 0, TRUE)
  rownames(wellCR) <- apply(wellCR, 1L, function(cr) paste0("C", cr[1L], "R", cr[2L]))
  
  MFt0 <- if(!multiple.res) res$MFbounds$origin["t"] else sapply(res, with, MFbounds$origin["t"])
  
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
        with(gwdata, approx(time + MFt0, colSums(-data[cr[1], cr[2],,, "Wells"]),
                            res$time, "constant", f = 1)$y)
      }else lapply(1:length(res), function(n){
        with(gwdata, approx(time + MFt0[n], colSums(-data[cr[1], cr[2],,, "Wells"]),
                            res[[n]]$time, "constant", f = 1)$y)
      })
      
      conc <- if(!multiple.res) J/Q*obsC.unitmatch else Map(`/`, lapply(J, `*`, obsC.unitmatch), Q)
      
      ctmp <- c(conc, recursive = TRUE)
      ylm <- c(0, max(ctmp[is.finite(ctmp)], obs$conc))
      
      if(multiple.res) lcols <- lighten(rep(lcol, length(res)),
                                        2^seq(1, -1, length.out = length(res)))
      
      plot((if(!multiple.res) res$time else res[[1]]$time)/365.25 + 1900,
           if(!multiple.res) conc else conc[[1L]],
           type = "l", lwd = 2, col = lcol,
           xlim = year.range, ylim = ylm, xlab = "date",
           ylab = paste0("abstracted concentration (", unit.label, ")"),
           main = w)
      if(multiple.res){
        Map(function(n, lc){
          lines(res[[n]]$time/365.25 + 1900, conc[[n]], lwd = 2, col = lc)
        }, 2:length(res), lcols[-1L])
        
        legend(legpos, legend = mr.labels, lwd = 2, col = lcols)
      }
      points(obs$time/365.25 + 1900, obs$conc, cex = 1.5, lwd = 2)
    }
  })
}
