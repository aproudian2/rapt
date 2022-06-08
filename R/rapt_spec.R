#
# This file contains methods for working with APT mass spectra
#

#### as.data.frame.MassSpectrum ####
as.data.frame.MassSpectrum <- function(M) {
  data.frame(mass = M@mass, intensity = M@intensity)
}

#### Spectral Ranging ####

### rangeCount ###
#' Count Hits in Range
#'
#' `rangeCount` counts the number of hits in the `POS` or `ATO` that fall within
#' the provided mass range. To get the hit information within a mass range, use
#' \code{\link{rangePOS}}.
#'
#' @param pos A `data.frame`. The `POS` to be ranged
#' @param start Numeric. The start of the mass range
#' @param end Numeric. The end of the mass range
#'
#' @return The number of hits falling within the range.
#'
#' @family ranging functions
#'
#' @export
rangeCount <- function(pos, start, end) {
  stopifnot(length(start) == 1)
  stopifnot(length(stop) == 1)
  n <- sum(pos$mass > start & pos$mass < end)
  return(n)
}

### rangePOS ###
#' Extract Hits in Range
#'
#' \code{rangePOS} extracts the rows of a `POS` or `ATO` object whose
#' mass is within the provided range. To count the number of hits within a mass
#' range, use \code{\link{rangeCount}}.
#'
#' @param pos The `POS` or `ATO` to be ranged. A `data.frame` created by
#'   \code{\link{readPOS}} or \code{\link{readATO}}
#' @param start Numeric. The start of the mass range
#' @param end Numeric. The end of the mass range
#'
#' @return A `data.frame` of the same structure as `pos` containing only hits
#' in the provided range
#'
#' @family ranging functions
#'
#' @export
rangePOS <- function(pos, start, end) {
  stopifnot(length(start) == 1)
  stopifnot(length(stop) == 1)
  pos[pos$mass > start & pos$mass < end,]
}

### rngCount ###
#' Count Hits in RRNG
#'
#' `rngCount` counts the number of hits for each entry within a `RRNG`
#' object. To get hit information within a `RRNG`, use
#' \code{\link{rngPOS}}.
#'
#' @param pos The `POS` or `ATO` to be ranged. A `data.frame` created by
#'   \code{\link{readPOS}} or \code{\link{readATO}}
#' @param rng The `RRNG` ranges to apply. A `data.frame` created by
#'   \code{\link{readRRNG}}
#' @param simplify logical. Whether to simplify the counts by ion name;
#'   default is `FALSE`
#'
#' @return A `data.frame` containing the name of each range, the number of
#'   counts and fraction of the total ranged counts.
#'
#' @family ranging functions
#'
#' @export
rngCount <- function(pos, rng, simplify = FALSE) {
  cts <- apply(rng, 1, function (x) {
    rangeCount(pos, as.numeric(x["start"]), as.numeric(x["end"]))
  })
  tot <- sum(cts)
  dat <- data.frame(name = rng$name, counts = cts, fraction = cts/tot)
  if(simplify) {
    dat <- aggregate(dat[,-1], by = list(name = dat$name), FUN = sum)
    dat <- dat[order(dat$counts, decreasing = TRUE),]
  }
  return(dat)
}

### rngPOS ###
#' Extract Hits in RRNG
#'
#' \code{rngPOS} extracts the rows of a \code{POS} or \code{ATO} object whose
#' mass is within the ranges of the provided \code{RRNG}. To count the hits
#' within a \code{RRNG}, use \code{\link{rngCount}}
#'
#' @param pos The `POS` or `ATO` from which to extract hits. A `data.frame`
#'   created by \code{\link{readPOS}} or \code{\link{readATO}}
#' @param rng The `RRNG` ranges to apply. A `data.frame` created by
#'   \code{\link{readRRNG}}
#'
#' @return A data.frame of the same structure as \code{pos}, with an appended
#'   column called "mark" that carries the name of the ion
#'
#' @family ranging functions
#'
#' @export
rngPOS <- function(pos, rng) {
  hits <- apply(rng, 1, function (x) {
    rows <- rangePOS(pos, as.numeric(x["start"]), as.numeric(x["end"]))
    rows$mark <- x["name"]
    return(rows)
  })
  dat <- do.call(rbind, hits)
  return(dat)
}

### rangeMassSpectrum ###
# Range peaks at a specified level
rangeMassSpectrum <- function(ms, start, end, threshold = 0.2) {
  spatstat.geom::verifyclass(ms, "MassSpectrum")
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
#' Initialize `nls` Fitting Using Ions
#'
#' \code{fitIonInit} creates the formula and starting values for fitting mass
#' spectra using \code{\link[stats]{nls}} based on molecular formulae. The
#' isotopic signatures are generated using \code{\link[enviPat]{isopattern}}.
#'
#' @param ions Character. A character vector of molecular formulae sutable to be
#'   passed to \code{\link[enviPat]{isopattern}}.
#' @param sd Numeric. The estimated standard deviation (or equivalent for
#'   `peak = "emg"`) of the fitted peak.
#' @param rel Numeric. A numeric vector of relative heights for each ion to
#'   initialize the starting values of the peak magnitudes.
#' @param peak Character. A string specifying what sort of peak shape to fit.
#'   Acceptable values are `"gaussian"` (the default) and `"emg"`. See Details.
#' @param mass.shift Numeric. The estimated mass shift of the measured mass
#'   spectrum relative to the absolute mass position of the ions.
#' @param noise Numeric. The estimated value of the constant noise floor of the
#'   mass spectrum.
#' @param tau Numeric. The estimated value of the skewness parameter for
#'   `peak = "emg"`. Ignored for `peak = "gaussian"`.
#' @param charge Integer. The charge(s) of the ion(s). Either a single value for
#'   all ions or a vector of values of the same length as `ions`. Defaults to 1
#' @param threshold Numeric. The lowest relative
#'
#' @return A list containing elements of the fitting formula, the starting fit
#'   values, and the lower limits for each parameter.
#'
#' @details
#' The fitting peak shape is either a gaussian (`peak = "gaussian"`) or
#' exponentially modified gaussian (`peak = "emg"`).  The `"gaussian"` peak
#' shape is defined as
#'   \deqn{rel[i] * \lambda[ij] *
#'     exp(-(mass - m[ij] - mass.shift)^2 / (2 * sd ^2))}
#'  where \eqn{\lambda[ij]} is the intensity of the *j*th isotope of `ions[i]`
#'  normalized to the most intense isotope of `ions[i]` and \eqn{m[ij]} is the
#'  mass of the *j*th isotope of `ions[i]`. The `"emg"` peak shape
#'  is defined as
#'    \deqn{(rel[i] * \lambda[ij] * sd / tau) * sqrt(pi/2) *
#'      exp(1/2*(sd/tau)^2-(mass-m[j]-mass.shift)/tau) *
#'      erfc(sqrt(1/2)*(sd/tau-(mass-m[j]-mass.shift)/sd))}
#'
#' @family spectral fitting functions
#'
#' @seealso \code{\link[stats]{nls}}, \code{\link[enviPat]{isopattern}}
#'
#' @export
# Add upper limits for port algorithm fitting
# Add ion specific sd / tau option
fitIonInit <- function(ions, sd = 0.5, rel = NULL, peak = "gaussian",
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
  if (peak == "emg") {
    ion.start <- c(ion.start, as.list(tau))
  }
  ion.lower <- rep_len(0, length(ion.start))
  ion.lower[names(ion.start) == "m0"] <- -1
  ion.lower[names(ion.start) == "t0"] <- 1e-6
  ion.init <- list(formula = ion.form, start = ion.start, lower = ion.lower)
  attr(ion.init, "isotopes") <- attr(ion.form, "isotopes")
  return(ion.init)
}

### ionFormula ###
#' Helper Function for \code{\link{fitIonInit}}
#'
#' Creates the formulae that get passed to \code{\link[stats]{nls}}
#'
#' @param ions Character. Molecular formula formatted in a way that is
#'   acceptable to \code{\link[enviPat]{isopattern}}
#' @param charge Numeric. Charge of ion.
#' @param threshold Numeric. Percentage threshold of isotope relative to highest
#'   intensity isotope. Defaults to 10
#' @param peak Character. One of `"gaussian"` or `"emg"`. Default is
#'   `peak = "gaussian"`. See Details.
#'
#' @return A character string that concatinates all the formulae for the ions
#'   and their respetive isotopes, along with a constant noise term.
#'
#' @details
#' The fitting peak shape is either a gaussian (`peak = "gaussian"`) or
#' exponentially modified gaussian (`peak = "emg"`).  The `"gaussian"` peak
#' shape is defined as
#'   \deqn{rel[i] * \lambda[ij] *
#'     exp(-(mass - m[ij] - mass.shift)^2 / (2 * sd ^2))}
#'  where \eqn{\lambda[ij]} is the intensity of the *j*th isotope of `ions[i]`
#'  normalized to the most intense isotope of `ions[i]` and \eqn{m[ij]} is the
#'  mass of the *j*th isotope of `ions[i]`. The `"emg"` peak shape
#'  is defined as
#'    \deqn{(rel[i] * \lambda[ij] * sd / tau) * sqrt(pi/2) *
#'      exp(1/2*(sd/tau)^2-(mass-m[j]-mass.shift)/tau) *
#'      erfc(sqrt(1/2)*(sd/tau-(mass-m[j]-mass.shift)/sd))}
#'
#' @family spectral fitting functions
#'
#' @seealso \code{\link[enviPat]{fitIonInit}}, \code{\link[enviPat]{isopattern}}
ionFormula <- function(ions, charge = 1, threshold = 10, peak = "gaussian") {
  data("isotopes", package = "enviPat", envir = environment())
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
      paste(paste0("rapt:::emg(mass, a", n, "*", I$lambda,", ",
                   I$mu, "-m0, s0, t0)", collapse = " + "))
    }, n.ions, ion.init, SIMPLIFY = FALSE)
  } else {
    stop('Peak shape must be one of "gaussian" or "emg."')
  }
  f.noise <- "n0"
  fs <- paste(paste(fs, collapse = ' + '), f.noise, sep = " + ")
  f.full <- paste("intensity ~", fs)
  attr(f.full, "isotopes") <- iso
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
# @family spectral fitting functions
fitGMM <- function(ms, start, stop) {

}
