

#' G features
#'
#' @description
#' Features calculated from G function (G3est).
#'
#' @details
#' Takes one set of radius values and corresponding G values of two functions as input and returns
#' features for machine learning
#'
#' @param rvals vector of radius values
#' @param gvals_new vector of new G3est values from which to subtract expected values
#' @param gvals_old vector of expected G3est values
#'
#' @return features to be used for machine learning prediction of clustering behavior.  The features are
#' the max difference between the functions (\emph{max_diff}), the r value at that difference (\emph{max_diff_r}),
#' the minimum difference (\emph{min_diff}), and the r value when the difference is 0 (\emph{zero_diff_r}).  Note that
#' (\emph{zero_diff_r}) must occur between the maximum and minimum difference
#'
#' @export
g3features <- function(rvals, gvals_new, gvals_old){
  diff <- gvals_new - gvals_old
  max_diff <- max(diff)
  max_diff_r <- rvals[which.max(diff)]
  min_diff <- min(diff)
  min_diff_r <- rvals[which.min(diff)]

  if (min_diff_r < max_diff_r) {
    zero_diff <- absmin(diff[(which.min(diff)):(which.max(diff))])
  }

  if (min_diff_r > max_diff_r) {
    zero_diff <- absmin(diff[(which.max(diff)):(which.min(diff))])
  }
  zero_diff_r <- rvals[which(diff == zero_diff)]
  out <- c(first(max_diff),first(max_diff_r), first(min_diff), first(zero_diff_r))
  return(out)
}

#' F features
#'
#' #' @description
#' Features calculated from F function (F3est).
#'
#' @details
#' Takes one set of radius values and corresponding F values of two functions as input and returns
#' features for machine learning
#'
#' @param rvals vector of radius values
#' @param fvals_new vector of new F3est values from which to subtract expected values
#' @param fvals_old vector of expected F3est values
#'
#' @return features to be used for machine learning prediction of clustering behavior.  The features are
#' the minimum difference between the functions (\emph{min_diff}), and the value of F (fvals_new) at
#' the minimum difference (\emph{min_diff_F})
#' @export
f3features = function(rvals, fvals_new, fvals_old) {
  diff = fvals_new -fvals_old

  min_diff = min(diff)
  min_diff_F = fvals_new[which.min(diff)]

  out = c(first(min_diff), first(min_diff_F))
  return(out)
}


#' GX features
#'
#' @description
#' Features calculated from G function (G3est).
#'
#' @details
#' Takes one set of radius values and corresponding G values of two functions as input and returns
#' features for machine learning
#'
#' @param rvals vector of radius values
#' @param gvals_new vector of new G3cross values from which to subtract expected values
#' @param gvals_old vector of expected G3cross values
#'
#' @return features to be used for machine learning prediction of clustering behavior.  The features are
#' the minimum difference between the functions (\emph{min_diff}),
#' the radius at which gvals_new = 0.95 \emph{ninty_fifth_percentile},
#' and the full width at half max of the peak formed by the difference in the functions \emph{FWHM}
#'
#' @export
g3Xfeatures<- function(rvals, gvals_new, gvals_old) {
  diff = gvals_new - gvals_old
  min_diff = min(diff)
  max_diff = max(diff)
  abs_max_ind = which.max(abs(diff))
  abs_max = diff[abs_max_ind]
  which.min(abs(0.95-gvals_new))
  ninty_fifth_percentile <- rvals[which.min(abs(0.95-gvals_new))]
  ind = which(diff ==absmax(diff))
  first = which(diff[1:ind]-absmax(diff)/2  ==
                  absmin(diff[1:ind] - absmax(diff)/2))
  second =ind +which(diff[ind:length(diff)]-absmax(diff)/2  ==
                       absmin(diff[ind:length(diff)] - absmax(diff)/2))
  FWHM = rvals[second] - rvals[first]
  out = c(first(min_diff), first(ninty_fifth_percentile), first(FWHM), first(abs_max), first(max_diff))
  return(out)
}

#### k3features ####
#' Extract features from 3D K function output
#' This is the same as k3metrics.  This function has update terminology ("metric" to "feature") to be more accurate
#'
#' Takes as inputs the results from \code{\link{anomK3est}} function performed
#' on a clustered data set. Returns five different features of this result.
#'
#' @param rvals_new The radius values from the \code{\link{anomK3est}}. These
#'   need to be trimmed so that they start past the inhibition radius of the
#'   underlying RCP pattern.
#' @param tvals_new The returned K values from the \code{\link{anomK3est}}.
#'   These need to be trimmed to the same indices as rvals_new.
#'
#' @return A list of [[1]] Km, [[2]] Rm, [[3]] Rdm, [[4]] Rddm, [[5]] Kdm for
#'   the input K function. If no first K peak is found, returns a list of 5
#'   NaNs.
#' @seealso \code{\link{anomK3est}}
#' @export
k3features <- function(rvals_new, tvals_new, toplot, ...) {
  if(any(is.infinite(tvals_new))){
    return(list(NA, NA, NA, NA, NA))
  }

  peak.info <- argmax(rvals_new, tvals_new, w = 3, span = 0.08)

  if(is.na(peak.info$x[1])){
    return(list(NA, NA, NA, NA, NA))
  }
  span <- (peak.info$x[1]/7)*(0.3)
  peak.info <- argmax(rvals_new, tvals_new, w = 3, span = span)
  peak.info$neg <- argmax(rvals_new, -1*tvals_new, w = 3, span = span)

  peak.info$deriv <- finite_deriv(rvals_new, peak.info$y.hat)

  peak.info$derivsm <- argmax(rvals_new, -1*peak.info$deriv, w = 3,
                              span = span)
  peak.info$derivsm_neg <- argmax(rvals_new, peak.info$deriv, w = 3,
                                  span = span)

  peak.info$dderiv <- finite_deriv(rvals_new, -1*peak.info$derivsm$y.hat)
  peak.info$dderivsm <- argmax(rvals_new, peak.info$dderiv, w = 3,
                               span = span)


  peak.info$ddderiv <- finite_deriv(rvals_new, peak.info$dderivsm$y.hat)
  peak.info$ddderivsm <- argmax(rvals_new, peak.info$ddderiv, w = 3,
                                span = span)

  lb <- peak.info$i[1]
  ub <- (peak.info$derivsm$i[1] + 2*peak.info$neg$i[1])/3
  if(is.na(ub)) {
    ub <- length(rvals_new)
  }
  Rddm_ind <- peak.info$ddderivsm$i[peak.info$ddderivsm$i >
                                      lb & peak.info$ddderivsm$i < ub][1]

  # stuff if you want to plot
  #browser()
  if(toplot == TRUE) {
    plot(rvals_new, tvals_new, type = "n", ...)

    lines(rvals_new, peak.info$y.hat, lwd = 2)
    points(peak.info$x[1], tvals_new[peak.info$i[1]], pch = 17, cex = 2, col= "gray60")
    #points(peak.info$neg$x, tvals_new[peak.info$neg$i], pch = 17, cex = 2, col="black")

    lines(rvals_new, -peak.info$derivsm$y.hat, lwd = 2, col = "red")
    points(peak.info$derivsm$x, peak.info$deriv[peak.info$derivsm$i],
           col="hotpink", cex = 2, pch = 17)

    lines(rvals_new, peak.info$dderivsm$y.hat, lwd = 2, col = "purple")
    #points(peak.info$dderivsm$x, peak.info$dderiv[peak.info$dderivsm$i],
    #col="green", cex = 2, pch = 17)

    lines(rvals_new, peak.info$ddderivsm$y.hat,lwd = 2, col = "blue")
    points(rvals_new[Rddm_ind], peak.info$ddderiv[Rddm_ind],
           col="deepskyblue", cex = 2, pch = 17)

    legend(max(rvals_new)*0.75, max(tvals_new)*0.9, legend = c('Tmax/Rmax', 'Tdmin/Rdmin','Rd3max'),
           col= c('gray60', 'hotpink', 'deepskyblue'), pch = 17)
  }

  return(list(peak.info$y.hat[peak.info$i[1]], #Km
              peak.info$x[1], #Rm
              peak.info$derivsm$x[1], #Rdm
              rvals_new[Rddm_ind[1]], #Rddm
              -peak.info$derivsm$y.hat[peak.info$derivsm$i[1]])) #Kdm
}

#' absmax and absmin
#' @description absolute maximum and minimum
#' @param x value to return absolute maximum or minimum of
#' @export
absmax = function(x){
  x[which.max(abs(x))]
}

#' absmax and absmin
#' @description absolute maximum and minimum
#' @param x value to return absolute maximum or minimum of
#' @export
#'
absmin = function(x){
  x[which.min(abs(x))]
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
#' @export

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
#' finds maximums and minimums in data based on these interpolating functions.
#' See
#' \url{https://stats.stackexchange.com/questions/36309/how-do-i-find-peaks-in-a-dataset}
#' @export

argmax <- function(x, y, w = 1, ...) {
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- zoo::rollapply(zoo::zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x = x[i.max], i = i.max, y.hat = y.smooth)
}

