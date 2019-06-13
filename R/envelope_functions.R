# Functions having to do with the calculation of and display of envelopes ran
# on random re-samplings of pp3 patterns

#### subSquare ####
#' Select a subsection from the center of a \code{\link[spatstat]{pp3}} data
#' file
#'
#' Given an original \code{\link[spatstat]{pp3}} object, \code{subSquare} will
#' select a rectangular prism centered at the center of the original point
#' pattern, and return a \code{\link[spatstat]{pp3}} object of the subsection.
#'
#' @param orig The original \code{\link[spatstat]{pp3}} object
#' @param win Numerical vector containing the dimensions for the box that you
#'   would like to select: c(xdim, ydim, zdim) (e.g. c(10,10,10)).
#' @return Returns a \code{\link[spatstat]{pp3}} object of the selected box,
#'   shifted so that the origin is still at (0,0,0).

subSquare <- function(orig, win) {

  orig.domain <- domain(orig)
  orig.center <- c(mean(orig.domain$xrange),
                   mean(orig.domain$yrange),
                   mean(orig.domain$zrange))
  xs <-orig.center[1] - (win[1]/2)
  ys <-orig.center[2] - (win[2]/2)
  zs <-orig.center[3] - (win[3]/2)
  xb <-orig.center[1] + (win[1]/2)
  yb <-orig.center[2] + (win[2]/2)
  zb <-orig.center[3] + (win[3]/2)

  xr <- c(xs, xb)
  yr <- c(ys, yb)
  zr <- c(zs, zb)
  sub.box <- box3(xrange = xr, yrange = yr, zrange = zr)

  tflist <- inside.boxx(orig, w = sub.box)

  sub <- orig[tflist]

  xrn <- c(0, xb-xs)
  yrn <- c(0, yb-ys)
  zrn <- c(0, zb-zs)
  sub.box.new <- box3(xrange = xrn, yrange = yrn, zrange = zrn)

  coo <- coords(sub)

  coo$x <- coo$x - xs
  coo$y <- coo$y - ys
  coo$z <- coo$z - zs

  sub.new <- createSpat(coo,win=sub.box.new)

  return(sub.new)
}

#### percentSelect ####
#' Randomly select a percent of the points in a \code{\link[spatstat]{pp3}}
#' object.
#'
#' Function randomly selects a certain percent of points within an original
#' \code{\link[spatstat]{pp3}} object. This function was created to be used in
#' random relabeling of point patterns.
#'
#' @param perc The fraction of points from the original pattern that are to be
#'   selected. A value between 0 and 1.
#' @param pattern The original \code{\link[spatstat]{pp3}} object to be selected
#'   from.
#' @param s Seed for the random selection
#' @return A \code{\link[spatstat]{pp3}} object containing only the selected
#'   points.
percentSelect <- function(perc, pattern) {

  reLabel <- rlabel(pattern,
                    labels = c(rep("A", round(npoints(pattern) * perc)),
    rep("B", npoints(pattern) - round(npoints(pattern) * perc))))
  inds <- which(marks(reLabel) == "A")

  newPattern <- reLabel[inds]
  return(newPattern)
}

#### envPlot ####
#' Plot envelopes of K3est test
#'
#' Plot the results of envelope calculations from the \code{\link{pK3est}},
#' \code{\link{pG3est}}, or \code{\link{pF3est}} functions, with the ability to
#' choose the percentiles for plotting.
#'
#' @param tests The return file from \code{p(K/G/F)3est} or the first, [[1]],
#'   entry in the list returned by \code{p(K/G/F)3est} with \code{anom = TRUE}.
#' @param percentiles Numerical vector of percentiles that you want to see the
#'   envelopes for. Each between 0 and 1.
#' @param ylim Numerical vector containing the min and max values for the y axis
#'   on the plot.
#' @param xlim Numerical vector containing the min and max values for the x axis
#'   on the plot.
#' @param leg True or falsel whether to show the automatically generated legend.
#' @param colors List of color names to make the envelopes.
#' @param ... Arguments to be passed into \code{plot()}.
#' @return Nothing.

envPlot <- function(tests, percentiles = c(0.999, 0.99, 0.97),
                    ylim = c(-3,3), xlim = c(0, ceiling(max(tests[,1]))),
                    leg = TRUE, ...) {
  color <- c("lightskyblue", "mediumpurple", "lightpink")

  color <- colors
  # break up data into r values and test results
  rvals <- tests[,1]
  tvals <- tests[,2:ncol(tests)]

  nTests <- ncol(tvals) # number of tests done
  prange <- percentiles * nTests # get the range of indeces for which each percentile spans

  sortedtVals <- t(apply(tvals, 1, sort)) # sort the results at each r value from lowest to highest
  percentileIndicesBig <- round(nTests/2) + floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
  percentileIndicesSmall <- round(nTests/2) - floor(prange/2) # do the same for the low end

  # grab out the columns from the sorted test results that we will plot
  toPlotBigs <- matrix(0, nrow = nrow(tvals), ncol = length(percentiles))
  toPlotSmalls <- matrix(0, nrow = nrow(tvals), ncol = length(percentiles))
  for(i in 1:length(percentiles)) {
    toPlotBigs[,i] <- sortedtVals[,percentileIndicesBig[i]]
    toPlotSmalls[,i] <- sortedtVals[,percentileIndicesSmall[i]]
  }

  # plot the envelopes from the percentile data
  #par(oma = c(0, 2, 0, 0))
  toplt <- data.frame(rvals,tvals[,1])

  plot(toplt, type = "n", xlab = "r",
       ylab = expression(sqrt('K'[3]*'(r)')*'  Anomaly'),
       ylim = ylim, xlim = xlim, xaxt = "n", ...)
  axis(1, at = 0:xlim[2], labels=FALSE)
  axis(1, at = seq(0,xlim[2],by=2))
  a <- c(rvals$V1, rev(rvals$V1))
  for(i in 1:length(percentiles)) {
    polygon(a, c(toPlotBigs[,i], rev(toPlotSmalls[,i])), col = color[i])#,border=color[i],lwd=2)
  }
  abline(h = 0, lty = 2, lwd = 1, col="black")
  if(leg == TRUE) {
    legend(0, ylim[2], legend = c(paste(toString(percentiles[1]*100), "% AI"),
                                  paste(toString(percentiles[2]*100), "% AI"),
                                  paste(toString(percentiles[3]*100), "% AI")),
           col=c(color[1], color[2], color[3]),
           lty = c(1,1,1), lwd = c(10,10,10))

  }
}

#### finite_deriv ####
#' Find the numerical derivative of a finite set of points.
#'
#' Uses the central difference method to calculate the numerical derivative for
#' a set of x,y data. Vectors x and y must be the same length.
#'
#' @param x The x values for your set of points.
#' @param y The y vales for your set of points.
#'
#' @return A vector with the same length as x containing estimated derivative at
#'   each x value.

finite_deriv <- function(x,y) {
  if(length(x) != length(y)) {
    print("x and y must be same length.")
    return()
  }
  d <- vector("numeric", length(x))
  for(i in 1:length(x)) {
    if(i == 1){
      d[i] <- (y[i+1] - y[i])/(x[i+1] - x[i])
    }else if(i == length(x)) {
      d[i] <- (y[i] - y[i-1])/(x[i] -x[i-1])
    }else {
      d[i] <- (y[i+1] - y[i-1])/(x[i+1] - x[i-1])
    }
  }
  return(d)
}

#### argmax ####
#' Find the peaks of a finite data set using smoothing.
#'
#' This function fits a rolling polynomial interpolation to a set of data and
#' finds maximums and minimums in data based on these interpolating functoins.
#' See
#' \url{https://stats.stackexchange.com/questions/36309/how-do-i-find-peaks-in-a-dataset}

argmax <- function(x, y, w = 1, ...) {
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x = x[i.max], i = i.max, y.hat = y.smooth)
}
