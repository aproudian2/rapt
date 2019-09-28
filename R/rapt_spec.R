#
# This file contains methods for working with the mass spectrum of APT data
#

#### rangeCount ####
#' Count the hits within a mass range
#'
#' To get the hit information within a mass range, use \code{\link{rangePOS}}
#' @param pos A data.frame. The pos to be ranged
#' @param start The start of the mass range
#' @param end The end of the mass range
#' @return The number of hits falling within the range.
#' @seealso \code{\link{rngCount}}
#' @export
rangeCount <- function(pos, start, end) {
  n <- with(pos, sum(mass > start & mass < end))
  return(n)
}

#### rangePOS ####
#' Extract hits within a given mass range.
#'
#' rangePOS extracts the rows of a \code{POS} or \code{ATO} object whose mass
#' is within the provided range.
#'
#' @param pos A data.frame. The pos to be ranged
#' @param start The start of the mass range
#' @param end The end of the mass range
#' @return A data.frame of the same structure as \code{pos} containing only hits
#' in the provided range
#' @seealso \code{\link{rngPOS}}
#' @export
rangePOS <- function(pos, start, end) {
  with(pos, pos[mass > start & mass < end,])
}

#### rngCount ####
#' Count the number of hits for each entry within a \code{RNG} object
#'
#' @param pos The pos to range
#' @param rng The ranges
#' @return A data.frame containing the name of each range, the number of counts
#' and fraction of the total ranged counts.
#' @seealso \code{\link{readRRNG}}, \code{\link{rangeCount}}
#' @export
rngCount <- function(pos, rng) {
  cts <- apply(rng, 1, function (x) {
    rangeCount(pos, as.numeric(x['start']), as.numeric(x['end']))
  })
  tot <- sum(cts)
  dat <- data.frame(name = rng$name, counts = cts, fraction = cts/tot)
  return(dat)
}

#### rngPOS ####
#' Extract hits according to a RNG and create a new pos
#'
#' rngPOS extracts the rows of a \code{POS} or \code{ATO} object whose mass
#' is within the ranges of the provided \code{RRNG}.
#'
#' @param pos The pos or ato to extract hits
#' @param rng The ranges to extract
#' @return A data.frame of the same structure as pos, with an appended column
#'   called "mark" that carries the name of the ion
#' @seealso \code{\link{readRRNG}}, \code{\link{rangePOS}}
#' @export
rngPOS <- function(pos, rng) {
  hits <- apply(rng, 1, function (x) {
    rows <- rangePOS(pos, as.numeric(x['start']), as.numeric(x['end']))
    rows$mark <- x['name']
    return(rows)
  })
  dat <- do.call(rbind, hits)
  return(dat)
}

#### rangeMassSpectrum ####
# Range peaks at a specified level
rangeMassSpectrum <- function(ms, start, end, threshold = 0.2) {
  m <- ms@mass
  m.in <- m >= start & m <= end
  m.clip <- m[m.in]
  int <- ms@intensity[m.in]
  int.max <- max(int)
  wh.all <- int >= threshold * int.max
  wh <- with(rle(wh.all), # this method is not right: range must include max!
             rep(lengths == max(lengths[values]) & values, lengths))
  range(m.clip[wh])
}

#### rangeMassPeaks ####

#### fitIonInit ####
# Create initialization for nls fitting of MassSpectrum peaks based on ions
fitIonInit <- function(ms, ions, sd = 0.5, rel = NULL,
                       mass.shift = 0, noise = 0, peak.intensity = NULL,
                       charge = NULL, threshold = 10) {
  if (is.null(rel)) {
    rel <- rep(1, length(ions))
  } else {
    stopifnot(length(rel) == length(ions))
  }
  if (is.null(charge)) {
    charge <- rep(1, length(ions))
  } else {
    stopifnot(length(charge) == length(ions))
  }
  if (is.null(peak.intensity)) {
    peak.intensity <- max(ms@intensity)
  }
  ion.form <- ionFormula(ions, charge = charge, threshold = threshold)
  names(sd) <- "s0"
  names(rel) <- paste0("a", seq_along(ions))
  names(mass.shift) <- "m0"
  names(noise) <- "n0"
  names(peak.intensity) <- "a0"
  ion.start <- c(as.list(sd), as.list(rel),
                 as.list(mass.shift), as.list(noise), as.list(peak.intensity))
  ion.init <- list(formula = ion.form, start = ion.start)
  return(ion.init)
}

#### ionFormula ####
# Helper function for fitIons
ionFormula <- function(ions, charge = 1, threshold = 10) {
  data('isotopes', package = 'enviPat', envir = environment())
  iso <- enviPat::isopattern(isotopes, ions,
                             charge = charge, threshold = threshold,
                             verbose = FALSE)
  ion.init <- mapply(function(X,s) {
    data.frame(
      lambda = c((X[,2]/100)),
      mu = c(X[,1])
    )
  }, iso, SIMPLIFY = FALSE)
  n.ions <- seq_along(ions)
  fs <- mapply(function(n, I, chg) {
    paste(paste0("a", n, " * ", I$lambda,
           " * exp(-(mass - ", I$mu, " - m0)**2 /",
           " (2 * (s0/", chg, ")**2))",
           collapse = " + "))
  }, n.ions, ion.init, charge, SIMPLIFY = FALSE)
  f.noise <- "n0"
  fs <- paste(paste(fs, collapse = ' + '), f.noise, sep = " + ")
  f.full <- paste("intensity ~ a0 * (", fs, ")")
  return(f.full)
}

#### fitGMM ####
# Fit a GMM to a region of a MassSpectrum
fitGMM <- function(ms, start, stop) {

}
