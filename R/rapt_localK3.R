#
# This file contains functions realting to the 3D local K function.
#

#### localK3est ####
#' Extends \code{\link[spatstat.core]{localK}} to 3D.
#'
#' Outputs local K function information for a \code{\link[spatstat.geom]{pp3}}
#' pattern. Currently only implemented for the translation edge correction.
#'
#' @param X A \code{\link[spatstat.geom]{pp3}} object to test.
#' @param rmax See \code{\link[spatstat.core]{K3est}}.
#' @param nrval See \code{\link[spatstat.core]{K3est}}.
#' @param correction Currently only "translation" implemented.
#'
#' @return Data frame with columns for the local K function around each point in
#'   the pattern.
#' @export

localK3est <- function(X, rmax = NULL, nrval = 128, correction = "translation") {
  stopifnot(inherits(X, "pp3"))
  correction <- pickoption("correction", correction,
    c(
      translation = "translation",
      trans = "translation",
      isotropic = "isotropic",
      iso = "isotropic",
      best = "isotropic"
    ),
    multi = TRUE
  )

  if (correction != "translation") {
    msg <- paste(
      "Local K3est function is currently only implemented for the",
      "translation edge correction."
    )
    print(msg)
    return()
  }

  B <- X$domain
  if (is.null(rmax)) {
    rmax <- diameter(B) / 2
  }
  r <- seq(from = 0, to = rmax, length.out = nrval)
  np <- npoints(X)

  # extract the x,y,z ranges as a vector of length 6
  flatbox <- unlist(B[1:3])

  # extract coordinates
  coo <- coords(X)

  u <- localk3engine(coo$x, coo$y, coo$z, flatbox,
    rmax = rmax, nrval = nrval
  )
  um <- matrix(u, nrow = nrval, ncol = np)
  kavg <- apply(um, 1, mean)

  rseq <- seq(0, rmax, len = nrval)
  rtheo <- (4 / 3) * pi * rseq^3
  K <- matrix(c(rseq, rtheo, kavg, u), nrow = nrval, ncol = (np + 3), byrow = FALSE)
  K <- as.data.frame(K)

  names <- vector("character", np + 3)
  names[1] <- "r"
  names[2] <- "theo"
  names[3] <- "kavg"
  for (i in 1:np) {
    names[i + 3] <- paste("k", toString(i), sep = "")
  }

  names(K) <- names

  return(K)
}

#### anomlocalK3est ####
#' Perform localK3est with 50th percentile of a RRL subtracted off.
#'
#' Similar to \code{\link{anomK3est}}, but returns the anomaly K3est for each
#' point in the pattern.
#'
#' @param X The \code{\link[spatstat.geom]{pp3}} object to be tested.
#' @param toSub The vector of values to subtract from the square root of the
#'   results of the K function applied to X.
#' @param rmax See \code{\link[spatstat.core]{K3est}}.
#' @param nrval See \code{\link[spatstat.core]{K3est}}.
#'
#' @return Date frame with columns of the anomaly K test for each point in the
#'   pattern.
#' @export

anomlocalK3est <- function(X, toSub, rmax, nrval) {
  a <- localK3est(X, rmax = rmax, nrval = nrval, correction = "translation")
  for (i in 2:ncol(a)) {
    a[, i] <- sqrt(a[, i]) - toSub
  }

  return(a)
}

#### localk3engine ####
#' Wrapper function for the C code under \code{\link{localK3est}}.
localk3engine <- function(x, y, z, box = c(0, 1, 0, 1, 0, 1), rmax = 1, nrval = 100) {
  res <- .C("RcallK3local",
    as.double(x), as.double(y), as.double(z),
    as.integer(length(x)),
    as.double(box[1L]), as.double(box[2L]),
    as.double(box[3L]), as.double(box[4L]),
    as.double(box[5L]), as.double(box[6L]),
    as.double(0), as.double(rmax),
    as.integer(nrval),
    f = as.double(numeric(nrval)),
    num = as.double(numeric(nrval)),
    denom = as.double(numeric(nrval)),
    full = as.double(numeric(nrval * length(x))),
    PACKAGE = "rapt"
  )
  return(res$full)
}
