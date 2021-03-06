% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rapt_cluster-sim.R
\name{rgblur}
\alias{rgblur}
\title{Helper for \code{\link{makecluster}} to 3D gaussian blur cluster positions.}
\usage{
rgblur(n = 1, mean = 0, sd = 1, coords = "rec", method = "r")
}
\arguments{
\item{n}{The number of points you want returned}

\item{mean}{Mean of the regular normal distribution for radius.}

\item{sd}{Standard deviation for the normal distribution for radius.}

\item{coords}{Return command. Either "rec" or "sph". See below for more.}

\item{method}{One of "r" and "xyz". "r" uniformly distributes direction and
normally distributes r, while "xyz" normally distributes x, y and z.}
}
\value{
If \code{coords = "rec"}, returns a vector of cartesian coordinates.
If \code{coords = "sph"}, returns a vector of spherical coordinates.
}
\description{
Function which returns a list of coordinates of radial gaussian blurred
points with uniform probability to be sent in any direction. This function
was written to apply a blur in positions of RCP spaced clusters. The radial
distribution is actually "folded normal" due to the fact that r > 0.
}
\references{
\url{https://en.wikipedia.org/wiki/Folded_normal_distribution}
}
\seealso{
\code{\link{makecluster}}
}
