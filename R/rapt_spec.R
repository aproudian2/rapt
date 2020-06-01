#
# This file contains methods for working with the mass spectrum of APT data
#

#### Spectral Ranging ####

### rangeCount ###
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
  n <- sum(pos$mass > start & pos$mass < end)
  return(n)
}

### rangePOS ###
#' Extract hits within a given mass range.
#'
#' rangePOS extracts the rows of a \code{POS} or \code{ATO} object whose mass
#' is within the provided range.
#'
#' @param pos A data.frame. The \code{POS} or \code{ATO} to be ranged
#' @param start The start of the mass range
#' @param end The end of the mass range
#' @return A data.frame of the same structure as \code{pos} containing only hits
#' in the provided range
#' @seealso \code{\link{rngPOS}}
#' @export
rangePOS <- function(pos, start, end) {
  pos[pos$mass > start & pos$mass < end,]
}

### rngCount ###
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

### rngPOS ###
#' Extract hits according to a RNG and create a new pos
#'
#' rngPOS extracts the rows of a \code{POS} or \code{ATO} object whose mass
#' is within the ranges of the provided \code{RRNG}.
#'
#' @param pos The pos or ato from which to extract hits
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

### rangeMassSpectrum ###
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

### rangeMassPeaks ###

#### Spectral Fitting ####

### fitIonInit ###
#' Create initialization for nls fitting of MassSpectrum peaks based on ions
#'
#' fitIonInit creates the formula and starting values for fitting mass spectra
#' using \code{\link[stats]{nls}} based on molecular formulae.
#'
#' @param ions Character. A character vector of molecular formulae sutable to be
#'   passed to \code{\link[enviPat]{isopattern}}.
#' @param sd Numeric. The estimated standard deviation (or equivalent for
#' peak = "emg") of the fitted peak.
#' @param rel Numeric. A numeric vector of relative heights for each ion to
#'   initialize the starting values of the peak magnitudes.
#' @param peak Character. A string specifying what sort of peak shape to fit.
#'   Acceptable values are "gaussian" (the default) and "emg". See details.
#' @param mass.shift Numeric. The estimated mass shift of the measured mass
#'   spectrum relative to the absolute mass position of the ions.
#' @param noise Numeric. The estimated value of the constant noise floor of the
#'   mass spectrum.
#' @param tau Numeric. The estimated value of the skewness parameter for
#'   peak = "emg". Ignored for peak = "gaussian".
#' @param charge Numeric.
#' @param threshold Numeric.
#'
#' @return A list containing elements of the fitting formula, the starting fit
#'   values, and the lower limits for each parameter.
#'
#' @details Fitting peak shape is either a gaussian (peak = "gaussian") or
#' exponentially modified gaussian (peak = "emg").
#'
#' @seealso \code{\link[stats]{nls}}, \code{\link[enviPat]{isopattern}}
#'
# Add upper limits for port algorithm fitting
# Add ion specific sd / tau option
fitIonInit <- function(ions, sd = 0.5, rel = NULL, peak = 'gaussian',
                       mass.shift = 0, noise = 0, tau = 0.1,
                       charge = NULL, threshold = 10) {
  if (is.null(rel)) {
    rel <- rep(1, length(ions))
  } else {
    stopifnot(length(rel) == length(ions))
  }
  if (is.null(charge)) {
    charge <- rep(1, length(ions))
  } else if (length(charge) == 1) {
    charge <- rep(charge, length(ions))
  } else {
    stopifnot(length(charge) == length(ions))
  }
  stopifnot(tau > 0)
  ion.form <- ionFormula(ions, charge = charge, threshold = threshold,
                         peak = peak)
  names(sd) <- "s0"
  names(rel) <- paste0("a", seq_along(ions))
  names(mass.shift) <- "m0"
  names(noise) <- "n0"
  names(tau) <- "t0"
  ion.start <- c(as.list(sd), as.list(rel),
                 as.list(mass.shift), as.list(noise))
  if (peak == 'emg') {
    ion.start <- c(ion.start, as.list(tau))
  }
  ion.lower <- rep_len(0, length(ion.start))
  ion.lower[names(ion.start) == "m0"] <- -1
  ion.lower[names(ion.start) == "t0"] <- 1e-6
  ion.init <- list(formula = ion.form, start = ion.start, lower = ion.lower)
  attr(ion.init, 'isotopes') <- attr(ion.form, 'isotopes')
  return(ion.init)
}

### ionFormula ###
# Helper function for fitIons
ionFormula <- function(ions, charge = 1, threshold = 10, peak = 'gaussian') {
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
  if (is.null(peak) | peak == 'gaussian') {
    fs <- mapply(function(n, I) {
      paste(paste0("a", n, " * ", I$lambda,
                   " * exp(-(mass - ", I$mu, " - m0)**2 /",
                   " (2 * s0**2))",
                   collapse = " + "))
    }, n.ions, ion.init, SIMPLIFY = FALSE)
  } else if (peak == 'emg') {
    fs <- mapply(function(n, I) {
      paste(paste0("emg(mass, a", n, "*", I$lambda,", ",
                   I$mu, "-m0, s0, t0)", collapse = " + "))
    }, n.ions, ion.init, SIMPLIFY = FALSE)
  } else {
    stop('Peak shape must be one of "gaussian" or "emg."')
  }
  f.noise <- "n0"
  fs <- paste(paste(fs, collapse = ' + '), f.noise, sep = " + ")
  f.full <- paste("intensity ~", fs)
  attr(f.full, 'isotopes') <- iso
  return(f.full)
}

### erfc ###
# Complementary error function
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)

### emg ###
# Exponentially modified gaussian function
emg <- function(x,h,m,s,t) {
  (h*s/t)*sqrt(pi/2)*exp(1/2*(s/t)^2-(x-m)/t)*erfc(sqrt(1/2)*(s/t-(x-m)/s))
}

#### fitGMM ####
# Fit a GMM to a region of a MassSpectrum
fitGMM <- function(ms, start, stop) {

}
