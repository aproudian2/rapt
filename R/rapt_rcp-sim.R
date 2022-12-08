#
# This file contains functions relating to dealing with RCP patterns.
#

### read.rcp ###
#' Read an RCP Simulation
#'
#' Reads in a RCP Generator output file as a \code{\link[spatstat.pp3]{pp3}} object.
#'
#' @param fpath_config The file path to the RCP FinalConfig file.
#' @param fpath_sys The file path to the associated RCP system file. Must only
#'   be specified if `scaleUp == TRUE`.
#' @param scaleUP Boolean. If `TRUE`, scales RCP so particles have new
#'   radius. If `FALSE` (the default), RCP stays as generated.
#' @param newRadius Numeric. If `scaleUP = TRUE`, this is the new radius
#' that the RCP particles will be scaled to. Default is 0.5.
#'
#' @details
#' The RCP generation is described in Desmond & Weeks,
#' "Random Close Packing of Disks and Spheres in Confined Geometries",
#' *Physical Review E*, **80** (5), 051305 (2009):
#' \url{https://doi.org/10.1103/PhysRevE.80.051305}; the algorithm can be
#' downloaded at
#' \url{http://www.physics.emory.edu/faculty/weeks/ken/RCPAlgorithm.html}.
#'
#' @return A \code{\link[spatstat.pp3]{pp3}} object of the RCP pattern.
#'
#' @family simulation functions
#'
#' @references Desmond & Weeks,
#'   "Random Close Packing of Disks and Spheres in Confined Geometries",
#'   \emph{Physical Review E}, \strong{80} (5), 051305 (2009):
#'   \url{https://doi.org/10.1103/PhysRevE.80.051305}
#' @seealso
#' \href{http://www.physics.emory.edu/faculty/weeks/ken/RCPAlgorithm.html}{RCP
#' Generator Algorithm}
#'
#' @export
read.rcp <- function(fpath_config, fpath_sys = NULL,
                     scaleUp = FALSE, newRadius = 0.5) {
  if (is.null(fpath_sys & scaleUp)) {
    stop("fpath_sys must be specified for scaling")
  }
  rcp <- read.table(fpath_config,
    sep = " ", col.names = c("x", "y", "z", "type")
  )
  pp <- createSpat(rcp[, c("x", "y", "z")])
  if (scaleUp == TRUE) {
    r <- scan(fpath_sys, skip = 7, n = 1)
    pp <- scaleRCP(pp, newRadius = newRadius, oldRadius = r)
  }
  return(pp)
}

#### scaleRCP ####
#' Scale RCP Data
#'
#' Given an original raw RCP output file, scale the 3D point cloud as if the
#' spheres have a specific radius.
#'
#' Input the original radius, taken from the "system" file output from the RCP
#' generation, and the new desired radius, `scaleRCP` will scale the RCP
#' point cloud as if the original generation hapened with spheres of the
#' specified radius.
#'
#' @param pp3file A \code{\link[spatstat.pp3]{pp3}} object containing the RCP
#'   generated 3D point positions.
#' @param newRadius The radius that you would like to scale to. Default is 0.5
#'   units
#' @param oldRadius The radius of the original RCP generation, taken from the
#'   "system" file output from the generation. Required input.
#' @param win A \code{\link[spatstat.geom]{box3}} object indicating the size of the
#'   original rcp generation. Can leave blank if these are integer value.
#' @return Will return a \code{\link[spatstat.pp3]{pp3}} object with the scaled
#'   point positions.
#'
#' @family simulation functions
#'
#' @export
scaleRCP <- function(pp3file, newRadius = 0.5, oldRadius = NULL, win = NULL) {
  if (is.null(oldRadius)) {
    print("Old radius is required.")
    return(NULL)
  }

  s <- (newRadius / oldRadius)

  if (is.null(win)) {
    pp3.domain <- domain(pp3file)
    pp3.box <- box3(
      xrange = s * round(pp3.domain$xrange),
      yrange = s * round(pp3.domain$yrange),
      zrange = s * round(pp3.domain$zrange)
    )
  } else {
    pp3.box <- box3(
      xrange = s * win$xrange,
      yrange = s * win$yrange,
      zrange = s * win$zrange
    )
  }
  rcp_xyz <- coords(pp3file) * s
  toReturn <- createSpat(rcp_xyz, win = pp3.box)

  return(toReturn)
}

#### stitch.size ###
# keep marks; add ability to specify by repetitions
#' Stitches together a RCP point pattern with periodic boundary conditions
#'
#' Similar to \code{\link{stitch}}. Instead of only inputting the number of
#' repetitions in each dimension, `stitch.size` allows you to specify the
#' domain size that you want to return, even if it is not an integer multiple of
#' the original dimensions.
#'
#' @param pp3file A \code{\link[spatstat.pp3]{pp3}} object containing the RCP
#'   generated 3D point positions.
#' @param win A \code{\link[spatstat.geom]{box3}} object indicating the size of the
#'   original rcp generation. Can leave blank.
#' @param boxSize A numeric vector of the dimensions of the final
#'   \code{\link[spatstat.pp3]{pp3}} object: c(xmax,ymax,zmax). Assumes that
#'   (xmin,ymin,zmin) = (0,0,0).
#' @return A \code{\link[spatstat.pp3]{pp3}} object.
#'
#' @family simulation functions
#'
#' @export
stitch.size <- function(pp3file, win = NULL, boxSize) {
  if (is.null(win)) {
    pp3.domain <- domain(pp3file)
    pp3.domain <- box3(
      xrange = pp3.domain$xrange,
      yrange = pp3.domain$yrange,
      zrange = pp3.domain$zrange
    )
  } else {
    pp3.domain <- win
  }

  reps <- ceiling(c(
    boxSize[1] / pp3.domain$xrange[2],
    boxSize[2] / pp3.domain$yrange[2],
    boxSize[3] / pp3.domain$zrange[2]
  ))
  pp3.box <- box3(
    xrange = reps[1] * pp3.domain$xrange,
    yrange = reps[2] * pp3.domain$yrange,
    zrange = reps[3] * pp3.domain$zrange
  )

  ogCoords <- as.matrix(coords(pp3file))
  n <- nrow(ogCoords)
  newpp3 <- matrix(NaN, n * reps[1] * reps[2] * reps[3], 3)
  ind <- 0

  for (i in 0:(reps[1] - 1)) {
    for (j in 0:(reps[2] - 1)) {
      for (k in 0:(reps[3] - 1)) {
        newpp3[(ind + 1):(ind + n), 1] <- ogCoords[, 1] +
          i * (pp3.domain$xrange[2] - pp3.domain$xrange[1])
        newpp3[(ind + 1):(ind + n), 2] <- ogCoords[, 2] +
          j * (pp3.domain$yrange[2] - pp3.domain$yrange[1])
        newpp3[(ind + 1):(ind + n), 3] <- ogCoords[, 3] +
          k * (pp3.domain$zrange[2] - pp3.domain$zrange[1])
        ind <- ind + n
      }
    }
  }

  newpp3 <- as.data.frame(newpp3)
  colnames(newpp3) <- c("x", "y", "z")

  toReturn <- createSpat(newpp3, win = pp3.box)

  toReturn2 <- subSquare(toReturn, win = boxSize)

  return(toReturn2)
}
