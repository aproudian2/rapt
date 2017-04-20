#
# This file contains methods for visualizing APT data.
#

#### polyCurve ####
#' Shade regions under a \code{\link[MALDIquant]{MassSpectrum}}.
#'
#' @seealso
#' \url{www.fromthebottomoftheheap.net/2013/01/11/shading-regions-under-a-curve}
#' for more details.
# Add ability to shade regions
polyCurve <- function(x, y, from, to, n = 50, miny,
                      col = "red", border = NA) {
  drawPoly <- function(fun, from, to, n = 50, miny, col, border) {
    Sq <- seq(from = from, to = to, length = n);
    polygon(x = c(Sq[1], Sq, Sq[n]),
            y = c(miny, fun(Sq), miny),
            col = col, border = border);
  }
  lf <- length(from);
  stopifnot(identical(lf, length(to)))
  if(length(col) != lf)
    col <- rep(col, length.out = lf);
  if(length(border) != lf)
    border <- rep(border, length.out = lf);
  if(missing(miny))
    miny <- min(y);
  interp <- approxfun(x = x, y = y);
  mapply(drawPoly, from = from, to = to, col = col, border = border,
         MoreArgs = list(fun = interp, n = n, miny = miny));
  invisible();
}
#### prettyPlot ####
#' Plot a pretty \code{\link[MALDIquant]{MassSpectrum}}.
#'
#' \code{prettyPlot}
# Lift the assumption of a log10 transformed spectrum.
prettyPlot <- function(ms, rng = NA,  xlim = NULL, ylim = NULL, lwd = 1,
                       main = "Mass Spectrum", hold.par = F) {
  if(log10(max(ms@intensity)) > 1) {
    stop("prettyPlot should only be used with log10 transformed spectra.")
  }
  ms.old <- par(no.readonly = T);
  par(mgp = c(2.2, 0.5, 0), cex = 1.5, las =  1, mar = c(3.5, 3.5, 1.5, 1));
  plot(ms, yaxt = "n", tck = -0.02, xlim = xlim, ylim = ylim, lwd = lwd,
       xlab = "Mass-to-Charge-State Ratio (Da)", ylab = "Counts (arb.)",
       main = main);
  ms.lims <- par("usr");
  ms.exp <- sapply(1:ms.lims[4], function(x){
    paste0(c("10^", x), collapse = "")
  });
  axis(side = 2, at = 1:ms.lims[4], labels = NA, tck = -0.02);
  axis(side = 2, at = 1:ms.lims[4], labels = parse(text = ms.exp),
       las = 1, lwd = 0);
  ms.tck <- lapply(1:2, axTicks);
  ms.diff <- ms.tck[[1]][2] - ms.tck[[1]][1];
  ms.tck[[1]] <- seq(
    min(c(ms.tck[[1]], 0)),
    max(ms.tck[[1]] + ms.diff),
    ms.diff / 5
  );
  ms.tck[[2]] <- log10(rep(1:9, ms.lims[4]+1) * rep(10^(0:ms.lims[4]),
                                                    each = 9));
  invisible(lapply(1:2, function(x){
    axis(x, at = ms.tck[[x]], labels = NA, tck = -0.01)
  }));
  if (!is.na(rng)) {
    polyCurve(ms@mass, ms@intensity, rng$start, rng$end, col = rng$color);
    lines(ms);
  }
  if (hold.par){
    par(ms.old);
  }
}
#### detectorGIF ####
#' Create a GIF that steps through detector hits.
#'
#' \code{detectorGIF} generates a GIF that shows the positions of detector hits.
#'
detectorGIF <- function(ato, len = 100000, range = NULL,
                        size = 7, name = 'detector', path = './', save.pdf = F,
                        delay = 10, ...) {
  gif.wd <- getwd()
  setwd(path)
  gif.pdf <- paste0(name, '.pdf')
  gif.gif <- paste0(name, '.gif')
  pdf(file = gif.pdf, width = size, height = size, bg = 'white')
  gif.n <- floor(dim(ato)[1]/len)
  for(i in 0:gif.n) {
    with(ato[i*len+seq_len(len),],
         plot(dy, dx, asp = 1, ... = ...))
  }
  dev.off()
  # convert pngs to one gif using ImageMagick
  gif.sys <- paste('convert - delay', delay, '-background white',
                   '-alpha background -density 100 +antialias',
                   gif.pdf, gif.gif)
  system(gif.sys)
  if(!save.pdf)
    file.remove(list.files(pattern = gif.pdf))
  setwd(gif.wd)
}

