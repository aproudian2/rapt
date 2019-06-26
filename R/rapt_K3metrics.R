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

k3metrics <- function(rvals.new, tvals.new, toplot) {
  peak.info <- argmax(rvals.new, tvals.new, w = 3, span = 0.08)

  if(any(is.infinite(tvals.new))){
    return(list(NA, NA, NA, NA, NA))
  }

  #lines(rvals.new, peak.info$y.hat)
  #points(peak.info$x,tvals.new[peak.info$i],pch = 6, cex = 2, col="red")
  if(is.na(peak.info$x[1])){
    return(list(NA, NA, NA, NA, NA))
  }
  span <- (peak.info$x[1]/7)*(0.3)
  peak.info <- argmax(rvals.new, tvals.new, w = 3, span = span)

  peak.info$deriv <- finite_deriv(rvals.new, peak.info$y.hat)
  #points(rvals.new,peak.info$deriv)
  peak.info$derivsm <- argmax(rvals.new, -1*peak.info$deriv, w = 3,
                              span = span)
  peak.info$derivsm_neg <- argmax(rvals.new, peak.info$deriv, w = 3,
                                  span = span)

  peak.info$dderiv <- finite_deriv(rvals.new, -1*peak.info$derivsm$y.hat)
  peak.info$dderivsm <- argmax(rvals.new, peak.info$dderiv, w = 3,
                               span = span)
  #points(rvals.new,peak.info$dderiv)

  peak.info$ddderiv <- finite_deriv(rvals.new, peak.info$dderivsm$y.hat)
  peak.info$ddderivsm <- argmax(rvals.new, peak.info$ddderiv, w = 3,
                                span = span)
  #points(rvals.new, peak.info$ddderiv)
  lb <- peak.info$i[1]
  ub <- (peak.info$derivsm$i[1] + peak.info$derivsm_neg$i[1])/2
  if(is.na(ub)) {
    ub <- length(rvals.new)
  }
  Rddm_ind <- peak.info$ddderivsm$i[peak.info$ddderivsm$i >
                                      lb & peak.info$ddderivsm$i < ub][1]

  # stuff if you want to plot
  #browser()
  if(toplot == TRUE) {
    plot(rvals.new, tvals.new, type = "n")

    lines(rvals.new, peak.info$y.hat)
    points(peak.info$x, tvals.new[peak.info$i], pch = 6, cex = 2, col="red")

    lines(rvals.new, -peak.info$derivsm$y.hat, lwd = 2, col = "red")
    points(peak.info$derivsm$x, peak.info$deriv[peak.info$derivsm$i],
           col="blue", cex = 2)

    lines(rvals.new, peak.info$dderivsm$y.hat, lwd = 2, col = "purple")
    points(peak.info$dderivsm$x, peak.info$dderiv[peak.info$dderivsm$i],
           col="green", cex = 2, pch = 19)

    lines(rvals.new, peak.info$ddderivsm$y.hat,lwd = 2, col = "blue")
    points(rvals.new[Rddm_ind], peak.info$ddderiv[Rddm_ind],
           col="green", cex = 1, pch = 19)
  }

  return(list(peak.info$y.hat[peak.info$i[1]], #Km
              peak.info$x[1], #Rm
              peak.info$derivsm$x[1], #Rdm
              rvals.new[Rddm_ind[1]], #Rddm
              -peak.info$derivsm$y.hat[peak.info$derivsm$i[1]])) #Kdm
}

#### localk3metrics ####
#' Extract variance of R_max over local K functions from a point pattern.
#'
#' This function takes the output of \code{\link{localK3est}}, calculates R_max
#' for the individual local K measures around each point, and returns the
#' variance of these measures.
#'
#' @param kl The output object from a \code{\link{localK3est}} function call.
#' @param start The x value index to start searching for mertics at (to avoid
#'   strange behavior near r = 0).
#' @param nsamp How many samples of R_max to take from the sample. If
#'   \code{NULL}, defaults to the number of points in the pattern (this can be
#'   very computationally expensive for large patterns).
#' @return A signle numerical value which is the variance of the R_max metric
#'   over \code{nsamp} individual points in the pattern.
#' @seealso \code{\link{localK3est}}, \code{\link[spatstat]{K3est}},
#'   \code{\link{kseries2}}

localk3metrics <- function(kl, start, nsamp = NULL){
  span <- NA
  i <- 1
  while(is.na(span)){
    r.n <- kl$r[start:nrow(kl)]
    t.n <- kl[start:nrow(kl),i+3]

    peak.info <- argmax(r.n, t.n, w = 3, span = 0.1)
    span <- (peak.info$x[1]/7)*(0.3)

    i <- i + 1
  }

  if(is.null(nsamp)){
    nsamp <- np
  }

  pks2 <- rep(0,nsamp)
  pts <- sample(1:np, nsamp, replace = FALSE)
  for(i in 1:nsamp){
    t.n <- kl[15:nrow(kl),pts[i]]
    a <- argmax(r.n, t.n, w = 3, span = span)
    pks2[i] <- a$x[1]
  }

  b <- var(pks2, na.rm = TRUE)
  return(b)
}

#### kseries2 ####
#' Helper function to allow for parallelization of metric extraction from
#' cluster simulations
#'
#' This function takes a single RCP pattern and creates and tests the K function
#' metrics on different cluster properties passed into the function. Used to
#' parallelize large sets of cluster property runs.
#'
#' @param j The index of the RCP pattern to run (note that this might only be
#'   able to be applied on the mines HPC under Galen's account right now as
#'   written)
#' @param p The total number of RCP patterns available.
#' @param tot A list of cluster properties to test. This should be a list of
#'   vectors containing c(r, den, rb, gbp).
#' @param maxr Maximum r to calculate K function to.
#' @param nr Number of r values to evaluate k function at.
#' @param toSub A vector of values to substitute for the \code{\link{anomK3est}}
#'   function.
#' @param rcp_path String holding the file path to the directory holding the RCP
#'   'FinalConfig' and 'system' files.
#' @param verbose \code{TRUE} or \code{FALSE}. Whether to output update files to
#'   a junk folder.
#' @param junk_path String holding the file path to the directory where the
#'   update files should go if \code{verbose = TRUE}.
#' @param s Random seed.
#'
#' @return Matrix containing 5 metrics and the seed value for each parameter
#'   combination given in the \code{tot} list.

kseries2 <- function(j, p ,tot, maxr, nr, toSub,
                     rcp_path = '~/Research/point_patterns/Final',
                     verbose = FALSE,
                     junk_path = '~/Research/junk/',
                     s){
  #t1 <- Sys.time()
  under.nums <- seq(2,(p+1),1)
  under.nums[length(under.nums)] <- 1
  over.nums <- seq(1,p,1)

  #upload
  under <- read.rcp(paste(rcp_path, '/FinalConfig', toString(under.nums[j]), sep=''),
                    paste(rcp_path, '/system', toString(under.nums[j]), sep=''),
                    scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste(rcp_path, '/FinalConfig', toString(over.nums[j]), sep=''),
                   paste(rcp_path, '/system', toString(over.nums[j]), sep=''),
                   scaleUp = TRUE,newRadius = 0.5)

  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))

  set.seed(s)
  cnt <- j*length(tot)*round(runif(1, 1, 10000))

  outtemp <- matrix(NA, nrow = length(tot), ncol = 5)

  for(i in 1:length(tot)){
    #print(i)
    if(length(tot[[1]]) == 4){
      cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",
                             cr=tot[[i]][1],
                             den = tot[[i]][2],
                             rb = TRUE, rbp = tot[[i]][1]*tot[[i]][3],
                             gb = TRUE, gbp = c(0, tot[[i]][4]),
                             s = cnt)
    }else{
      cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",
                             cr=tot[[i]][1],
                             den = tot[[i]][2],
                             gb = TRUE, gbp = c(0, tot[[i]][3]),
                             s = cnt)
    }
    if(is.numeric(cluster)){
      outtemp[i,] <- c(NA, NA, NA, NA, NA)
      a <- as.data.frame(1)
      if(verbose == TRUE){
        if(i > 1){
          file.remove(paste(junk_path, '/', toString(j), '_', toString(i-1), '.csv', sep = ''))
        }
        write.csv(a, file = paste(junk_path, '/', toString(j), '_', toString(i), '.csv', sep = ''))
      }
      next
    }
    result <- anomK3est(cluster[[1]],toSub,maxr,nr)
    rvals <- result$r
    tvals <- result$trans

    # get out that peak info son
    rvals.new <- rvals[15:length(rvals)]
    tvals.new <- tvals[15:length(rvals)]

    #get those metrics out
    metrics <- k3metrics(rvals.new, tvals.new, FALSE)

    outtemp[i,] <- c(metrics[[1]], metrics[[2]], metrics[[3]], metrics[[4]], metrics[[5]])

    rm(cluster, result, rvals, tvals, rvals.new, tvals.new)
    gc()

    cnt <- cnt + 1

    a <- as.data.frame(1)

    if(verbose == TRUE){
      if(i > 1){
        file.remove(paste(junk_path, '/', toString(j), '_', toString(i-1), '.csv', sep = ''))
      }
      write.csv(a, file = paste(junk_path, '/', toString(j), '_', toString(i), '.csv', sep = ''))
    }
  }

  rm(over, under, over.big, under.big)
  gc()

  return(outtemp)
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
  y.max <- zoo::rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x = x[i.max], i = i.max, y.hat = y.smooth)
}
