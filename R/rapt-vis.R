### Mass Spectrum ###
# Shade regions under a MassSpectrum
# Based on Gavin Simpson:
# www.fromthebottomoftheheap.net/2013/01/11/shading-regions-under-a-curve
polyCurve <- function(x, y, from, to, n = 50, miny,
                      col = "red", border = col) {
  drawPoly <- function(fun, from, to, n = 50, miny, col, border) {
    Sq <- seq(from = from, to = to, length = n)
    polygon(x = c(Sq[1], Sq, Sq[n]),
            y = c(miny, fun(Sq), miny),
            col = col, border = border)
  }
  lf <- length(from)
  stopifnot(identical(lf, length(to)))
  if(length(col) != lf)
    col <- rep(col, length.out = lf)
  if(length(border) != lf)
    border <- rep(border, length.out = lf)
  if(missing(miny))
    miny <- min(y)
  interp <- approxfun(x = x, y = y)
  mapply(drawPoly, from = from, to = to, col = col, border = border,
         MoreArgs = list(fun = interp, n = n, miny = miny))
  invisible()
}
# Plot a pretty MassSpectrum
prettyPlot <- function(ms, rng = NA,  xlim = NULL, ylim = NULL) {
  ms.old <- par(no.readonly = T)
  par(mgp = c(1.9, 0.5, 0), cex = 1.5, las =  1, mar = c(3, 3, 1.5, 0.5))
  plot(ms, yaxt = "n", tck = -0.02, xlim = xlim, ylim = ylim,
       xlab = "Mass-to-Charge-State Ratio (Da)", ylab = "Counts (arb.)",
       main = "Mass Spectrum")
  ms.lims <- par("usr")
  ms.exp <- sapply(1:ms.lims[4], function(x){paste0(c("10^", x), collapse = "")})
  axis(side = 2, at = 1:ms.lims[4], labels = NA, tck = -0.02)
  axis(side = 2, at = 1:ms.lims[4], labels = parse(text = ms.exp), las = 1, lwd = 0)
  ms.tck <- lapply(1:2, axTicks)
  ms.tck[[1]] <- seq(min(ms.tck[[1]]), max(ms.tck[[1]]), (ms.tck[[1]][2] - ms.tck[[1]][1]) / 5)
  ms.tck[[2]] <- log10(rep(1:9, ms.lims[4]+1) * rep(10^(0:ms.lims[4]), each = 9))
  invisible(lapply(1:2, function(x){axis(x, at = ms.tck[[x]], labels = NA, tck = -0.01)}))
  if (!is.na(rng)) {
    rng <- do.call(rbind,unlist(rng, recursive = F))
    polyCurve(ms@mass,ms@intensity,rng[,1],rng[,2])
    lines(ms)
  }
  par(ms.old)
}
