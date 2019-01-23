#Functions to analze and produce results from data

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

k3metrics <- function(rvals.new, tvals.new){

  #plot(rvals.new,tvals.new, type = "n", ...)
  peak.info <- argmax(rvals.new, tvals.new, w = 3, span = 0.08)
  #lines(rvals.new, peak.info$y.hat)
  #points(peak.info$x,tvals.new[peak.info$i],pch = 6, cex = 2, col="red")
  if(is.na(peak.info$x[1])){
    return(list(NaN, NaN, NaN, NaN, NaN))
  }
  span <- (peak.info$x[1]/7)*(0.3)
  peak.info <- argmax(rvals.new, tvals.new, w = 3, span = span)
  #lines(rvals.new, peak.info$y.hat)
  #points(peak.info$x,tvals.new[peak.info$i],pch = 6, cex = 2, col="red")

  peak.info$deriv <- finite_deriv(rvals.new, peak.info$y.hat)
  #points(rvals.new,peak.info$deriv)
  peak.info$derivsm <- argmax(rvals.new, -1*peak.info$deriv, w = 3, span = span)
  peak.info$derivsm_neg <- argmax(rvals.new, peak.info$deriv, w = 3, span = span)
  #lines(rvals.new, -peak.info$derivsm$y.hat, lwd = 2, col = "red")
  #points(peak.info$derivsm$x,peak.info$deriv[peak.info$derivsm$i], col="blue", cex = 2)

  peak.info$dderiv <- finite_deriv(rvals.new, -1*peak.info$derivsm$y.hat)
  peak.info$dderivsm <- argmax(rvals.new, peak.info$dderiv, w = 3, span = span)
  #points(rvals.new,peak.info$dderiv)
  #lines(rvals.new, peak.info$dderivsm$y.hat, lwd = 2, col = "purple")
  #points(peak.info$dderivsm$x,peak.info$dderiv[peak.info$dderivsm$i], col="green", cex = 2, pch = 19)

  peak.info$ddderiv <- finite_deriv(rvals.new, peak.info$dderivsm$y.hat)
  peak.info$ddderivsm <- argmax(rvals.new, peak.info$ddderiv, w = 3, span = span)
  #points(rvals.new, peak.info$ddderiv)
  lb <- peak.info$i[1]
  ub <- (peak.info$derivsm$i[1] + peak.info$derivsm_neg$i[1])/2
  if(is.na(ub)){
    ub <- length(rvals.new)
  }
  Rddm_ind <- peak.info$ddderivsm$i[peak.info$ddderivsm$i > lb & peak.info$ddderivsm$i < ub][1]
  #lines(rvals.new, peak.info$ddderivsm$y.hat,lwd = 2, col = "blue")
  #points(rvals.new[Rddm_ind],peak.info$ddderiv[Rddm_ind], col="green", cex = 1, pch = 19)

  return(list(peak.info$y.hat[peak.info$i[1]], #Km
              peak.info$x[1], #Rm
              peak.info$derivsm$x[1], #Rdm
              rvals.new[Rddm_ind[1]], #Rddm
              -peak.info$derivsm$y.hat[peak.info$derivsm$i[1]])) #Kdm
}
