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
#' @seealso \code{\link{rngCount}}, \code{\link{rangePOS}}
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
#' @return A data.frame of the same structure as \code{pos} containing only hits
#' in the provided range
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
#' @param pos The pos or ato to extract hits
#' @param rng The ranges to extract
#' @return A data.frame of the same structure as pos
#' @seealso \code{\link{rangePOS}}
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

#### autoRange ####
# Range peaks at a specified level
autoRange <- function(pos, rng, width = 10) {

}
