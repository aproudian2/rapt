#
# This file contains methods for importing and conditioning APT data
#

#### Load Data ####
# Create generic for .ato and .pos
#' Read a POS file.
#'
#' \code{readPOS} reads a POS file (from IVAS) into a data frame.
#'
#' @param filepath A string. The file path to the POS file.
#' @return A dataframe with columns corresponding to the x, y, and z positions
#'   of the reconstruction and the mass-to-charge ratio.
#' @references Local Electrode Atom Probe Tomography: A User's Guide
#' @seealso \code{\link{readATO}}
#' @example X <- readPOS("R45_0001-v01.pos")
#' @export
readPOS <- function(filepath) {
  pos.len <- file.info(filepath)['size'] / 4
  pos.len <- as.numeric(pos.len)
  pos.raw <- readBin(
    filepath, what = 'numeric',
    size = 4, endian = 'big', n = pos.len
  )
  pos.mat <- matrix(pos.raw, ncol = 4, byrow = T)
  pos.dat <- as.data.frame(pos.mat)
  names(pos.dat) <- c("x", "y", "z", "mass")
  pos.name <- sub(".pos", "", basename(filepath))
  attr(pos.dat, "metaData") <- list(
    name = pos.name
  )
  return(pos.dat);
}

#' Read an ATO file.
#'
#' \code{readATO} reads an ATO file (from IVAS) into a data frame.
#'
#' @param filepath A string. The file path to the ATO file.
#' @return A dataframe.
#' @references Local Electrode Atom Probe Tomography: A User's Guide
#' @seealso \code{\link{readPOS}}
#' @example X <- readATO("R45_0001-v01.ato")
#' @export
readATO <- function(filepath) {
  ato.len <- file.info(filepath)['size'] / 4
  ato.len <- as.numeric(ato.len)
  ato.file <- file(filepath, open = 'rb')
  seek(ato.file, where = 8)
  ato.raw <- readBin(
    ato.file, what = 'numeric',
    size = 4, endian = 'little', n = ato.len
  )
  ato.mat <- matrix(ato.raw, ncol = 14, byrow = T)
  ato.dat <- as.data.frame(ato.mat)
  names(ato.dat) <- c("x", "y", "z", "mass", "clusID", "pIndex", "Vdc",
                      "TOF", "dx", "dy", "Vp", "shank", "FouR", "FouI");
  ato.name <- sub(".ato", "", basename(filepath))
  attr(ato.dat, "metaData") <- list(
    name = ato.name
  )
  close(ato.file)
  return(ato.dat)
}

#' Create a list of mass spec formulae from a table (class \code{cform}).
#'
#' \code{readForm} reads a CSV file of formulae into a data frame of class
#' \code{cform}. The CSV should have headers specifying the name of the
#' compound, its chemical formula, and optionally the color used for plotting.
readForm <- function(filepath) {
  form.dat <- read.csv(filepath, stringsAsFactors = F)
  form.chk <- which(
    enviPat::check_chemform(rapt::isotopes, form.dat[, 'formula'])[, 'warning'])
  if(length(form.chk) == 0) {
    form.dat <- form.dat
    form.dat <- as.data.frame(form.dat)
    class(form.dat) <- c("data.frame", "cform")
    return(form.dat)
  }else {
    print(paste("The formulae:", form.dat[form.chk, 'formula'],
                "are not valid."))
  }
}

#' Read an RNG file (from IVAS) into a RNG object
readRNG <- function(rng) {

}

#### Condition Data ####
## Add ability to specify marks
#' Create a \code{\link[spatstat]{pp3}} object from a POS or ATO data frame.
#'
#' @param pos A POS or ATO data frame.
#' @param win The domain of the data.
#' @return A \code{\link[spatstat]{pp3}} with the x,y,z positions of the hits in
#'   the supplied POS or ATO.
#' @seealso \code{\link{readPOS}}, \code{\link{readATO}}
#' @export
createSpat <- function(pos, win = NULL) {
  pp3.box <- win;
  if(is.null(win)) {
    pp3.box <- sapply(pos[1:3], range);
  }
  pp3.dat <- pp3(pos$x, pos$y, pos$z, pp3.box);
  attr(pp3.dat, "metaData") <- attr(pos, "metaData");
  return(pp3.dat);
}

#' Create a \code{\link[spatstat]{ppp}} from an ATO.
#'
#' \code{createDet} generates a \code{\link[spatstat]{ppp}} of detector hits
#' from an ATO.
#'
#' @param ato An ATO data frame.
#' @param window An object of class \code{\link[spatstat]{owin}}. If NULL,
#'   a window will be calculated from the data using
#'   \code{\link[spatstat]{ripras}}.
#' @return A \code{\link[spatstat]{ppp}} with the positions of the detector hits
#'   from the ATO.
#' @seealso \code{\link{readATO}}, \code{\link[spatstat]{ppp}},
#'   \code{\link[spatstat]{ripras}}
#' @export
createDet <- function(ato, window = NULL) {
  if(is.null(window))
    window <- ripras(ato$dx, ato$dy)
  det.dat <- ppp(ato$dx, ato$dy, window = window)
  return(det.dat)
}

#' Create a \code{\link[MALDIquant]{MassSpectrum}} from a POS or ATO
#' data frame.
#'
#' \code{createSpec} generates a \code{\link[MALDIquant]{MassSpectrum}} object
#' with a specified resolution from an ATO or POS data frame (like that created
#' by \code{\link(readPOS)}).
#'
#' @param pos A POS or ATO data frame.
#' @param res The desired \code{\link[MALDIquant]{MassSpectrum}} resolution.
#' @return A \code{\link[MALDIquant]{MassSpectrum}} from the \code{mass} field
#'   of the POS or ATO, with the resolution set by \code{res}.
#'
#' @details
#' The input POS or ATO is binned by mass values; the resolution parameter sets
#' the width of the mass bins used in \code{\link[base]hist}} to create the
#' input to the \code{\link[MALDIquant]{createMassSpectrum}} call, and also acts
#' as a tolerance around the spectrum minimum and maximum mass. The minimum of
#' the mass value is zero.
#' @seealso \code{\link{readPOS}}, \code{\link{readATO}},
#'   \code{\link[MALDIquant]{MassSpectrum}}
#' @export
createSpec <- function(pos, res = 0.001, snip = NULL) {
  ms.max <- max(pos[,"mass"])
  ms.max <- ms.max + res
  ms.min <- min(pos[,"mass"])
  ms.min <- ms.min - res
  if(ms.min < 0) {
    ms.min <- 0
  }
  ms.breaks <- seq(ms.min, ms.max, res)
  ms.hist <- hist(pos[,"mass"], ms.breaks, plot = F)
  ms.dat <- createMassSpectrum(
    ms.hist$mids[-1], ms.hist$counts[-1],
    metaData = attr(pos, "metaData")
  )
  return(ms.dat)
}

#' Create a TOF spectrum from an ATO data frame
createTOF <- function(pos, res = 0.001) {
  ms.max <- max(pos[,"TOF"])
  ms.max <- ms.max + res
  ms.min <- min(pos[,"TOF"])
  ms.min <- ms.min - res
  ms.breaks <- seq(ms.min, ms.max, res)
  ms.hist <- hist(pos[,"TOF"], ms.breaks, plot = F)
  ms.dat <- createMassSpectrum(
    ms.hist$mids, ms.hist$counts,
    metaData = attr(pos, "metaData")
  )
  return(ms.dat)
}

createForm <- function(df) {
  form.chk <- which(
    check_chemform(isotopes, df[, "formula"])[, "warning"]);
  if(length(form.chk) == 0) {
    class(df) <- c("data.frame", "cform")
    return(df);
  }else {
    print(paste("The formulae:", form.dat[form.chk], "are not valid."));
  }
}

#### Write Data ####
# Export methods?
