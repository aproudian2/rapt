#
# This file contains functions relating to the calculation of metrics from the 3D K function
#

#### k3metrics ####
#' Extract metrics from 3D K function output
#'
#' Takes as inputs the results from \code{\link{anomK3est}} function performed
#' on a clustered data set. Returns five different metrics of this result.
#'
#' @param rvals.new The radius values from the \code{\link{anomK3est}}. These
#'   need to be trimmed so that they start past the inhibition radius of the
#'   underlying RCP pattern.
#' @param tvals.new The returned K values from the \code{\link{anomK3est}}.
#'   These need to be trimmed to the same indices as rvals.new.
#'
#' @return A list of [[1]] Km, [[2]] Rm, [[3]] Rdm, [[4]] Rddm, [[5]] Kdm for
#'   the input K function. If no first K peak is found, returns a list of 5
#'   NaNs.
#' @seealso \code{\link{anomK3est}}
#' @export

k3metrics <- function(rvals.new, tvals.new, toplot, ...) {
  if(any(is.infinite(tvals.new))){
    return(list(NA, NA, NA, NA, NA))
  }

  peak.info <- argmax(rvals.new, tvals.new, w = 3, span = 0.08)

  if(is.na(peak.info$x[1])){
    return(list(NA, NA, NA, NA, NA))
  }
  span <- (peak.info$x[1]/7)*(0.3)
  peak.info <- argmax(rvals.new, tvals.new, w = 3, span = span)
  peak.info$neg <- argmax(rvals.new, -1*tvals.new, w = 3, span = span)

  peak.info$deriv <- finite_deriv(rvals.new, peak.info$y.hat)

  peak.info$derivsm <- argmax(rvals.new, -1*peak.info$deriv, w = 3,
                              span = span)
  peak.info$derivsm_neg <- argmax(rvals.new, peak.info$deriv, w = 3,
                                  span = span)

  peak.info$dderiv <- finite_deriv(rvals.new, -1*peak.info$derivsm$y.hat)
  peak.info$dderivsm <- argmax(rvals.new, peak.info$dderiv, w = 3,
                               span = span)


  peak.info$ddderiv <- finite_deriv(rvals.new, peak.info$dderivsm$y.hat)
  peak.info$ddderivsm <- argmax(rvals.new, peak.info$ddderiv, w = 3,
                                span = span)

  lb <- peak.info$i[1]
  ub <- (peak.info$derivsm$i[1] + 2*peak.info$neg$i[1])/3
  if(is.na(ub)) {
    ub <- length(rvals.new)
  }
  Rddm_ind <- peak.info$ddderivsm$i[peak.info$ddderivsm$i >
                                      lb & peak.info$ddderivsm$i < ub][1]

  # stuff if you want to plot
  #browser()
  if(toplot == TRUE) {
    plot(rvals.new, tvals.new, type = "n", ...)

    lines(rvals.new, peak.info$y.hat, lwd = 2)
    points(peak.info$x[1], tvals.new[peak.info$i[1]], pch = 17, cex = 2, col= "gray60")
    #points(peak.info$neg$x, tvals.new[peak.info$neg$i], pch = 17, cex = 2, col="black")

    lines(rvals.new, -peak.info$derivsm$y.hat, lwd = 2, col = "red")
    points(peak.info$derivsm$x, peak.info$deriv[peak.info$derivsm$i],
           col="hotpink", cex = 2, pch = 17)

    lines(rvals.new, peak.info$dderivsm$y.hat, lwd = 2, col = "purple")
    #points(peak.info$dderivsm$x, peak.info$dderiv[peak.info$dderivsm$i],
           #col="green", cex = 2, pch = 17)

    lines(rvals.new, peak.info$ddderivsm$y.hat,lwd = 2, col = "blue")
    points(rvals.new[Rddm_ind], peak.info$ddderiv[Rddm_ind],
           col="deepskyblue", cex = 2, pch = 17)

    legend(max(rvals.new)*0.75, max(tvals.new)*0.9, legend = c('Tmax/Rmax', 'Tdmin/Rdmin','Rd3max'),
           col= c('gray60', 'hotpink', 'deepskyblue'), pch = 17)
  }

  return(list(peak.info$y.hat[peak.info$i[1]], #Km
              peak.info$x[1], #Rm
              peak.info$derivsm$x[1], #Rdm
              rvals.new[Rddm_ind[1]], #Rddm
              -peak.info$derivsm$y.hat[peak.info$derivsm$i[1]])) #Kdm
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
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- zoo::rollapply(zoo::zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x = x[i.max], i = i.max, y.hat = y.smooth)
}
