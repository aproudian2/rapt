#
# This file contains functions relating to the manipulation of pp3 objects.
#

#### subSquare ####
#' Select a subsection from the center of a \code{\link[spatstat]{pp3}} data
#' file.
#'
#' Given an original \code{\link[spatstat]{pp3}} object, \code{subSquare} will
#' select a rectangular prism centered at the center of the original point
#' pattern, and return a \code{\link[spatstat]{pp3}} object of the subsection.
#'
#' @param orig The original \code{\link[spatstat]{pp3}} object
#' @param win Numerical vector containing the dimensions for the box that you
#'   would like to select: c(xdim, ydim, zdim) (e.g. c(10,10,10)).
#' @return Returns a \code{\link[spatstat]{pp3}} object of the selected box,
#'   shifted so that the origin is still at (0,0,0).
#' @export

subSquare <- function(orig, win) {

  orig.domain <- domain(orig)
  orig.center <- c(mean(orig.domain$xrange),
                   mean(orig.domain$yrange),
                   mean(orig.domain$zrange))
  xs <-orig.center[1] - (win[1]/2)
  ys <-orig.center[2] - (win[2]/2)
  zs <-orig.center[3] - (win[3]/2)
  xb <-orig.center[1] + (win[1]/2)
  yb <-orig.center[2] + (win[2]/2)
  zb <-orig.center[3] + (win[3]/2)

  xr <- c(xs, xb)
  yr <- c(ys, yb)
  zr <- c(zs, zb)
  sub.box <- box3(xrange = xr, yrange = yr, zrange = zr)

  tflist <- inside.boxx(orig, w = sub.box)

  sub <- orig[tflist]

  xrn <- c(0, xb-xs)
  yrn <- c(0, yb-ys)
  zrn <- c(0, zb-zs)
  sub.box.new <- box3(xrange = xrn, yrange = yrn, zrange = zrn)

  coo <- coords(sub)

  coo$x <- coo$x - xs
  coo$y <- coo$y - ys
  coo$z <- coo$z - zs

  sub.new <- createSpat(coo,win=sub.box.new)

  return(sub.new)
}

#### percentSelect ####
#' Randomly select a percent of the points in a \code{\link[spatstat]{pp3}}
#' object.
#'
#' Function randomly selects a certain percent of points within an original
#' \code{\link[spatstat]{pp3}} object. This function was created to be used in
#' random relabeling of point patterns.
#'
#' @param perc The fraction of points from the original pattern that are to be
#'   selected. A value between 0 and 1.
#' @param pattern The original \code{\link[spatstat]{pp3}} object to be selected
#'   from.
#' @param s Seed for the random selection
#' @return A \code{\link[spatstat]{pp3}} object containing only the selected
#'   points.
#' @export

percentSelect <- function(perc, pattern) {

  reLabel <- rlabel(pattern,
                    labels = c(rep("A", round(npoints(pattern) * perc)),
                               rep("B", npoints(pattern) - round(npoints(pattern) * perc))))
  inds <- which(marks(reLabel) == "A")

  newPattern <- reLabel[inds]
  return(newPattern)
}
