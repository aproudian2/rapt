% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rapt_extend.R
\name{quantess.pp3}
\alias{quantess.pp3}
\title{Quantile Tessellation}
\usage{
\method{quantess}{pp3}(M, Z, n, ..., type = 2, origin = c(0, 0), eps = NULL)
}
\arguments{
\item{M}{A pp3}

\item{Z}{A spatial covariate (a pixel image or a function(x,y,z)) or one of
the strings \code{"x"}, \code{"y"}, or \code{"z"} indicating the Cartesian coordinates
\emph{x}, \emph{y}, or \emph{z} or one of the strings \code{"rad"} or \code{"ang"} indicating polar
coordinates. The range of values of \code{Z} will be broken into n bands
containing equal amounts of stuff.}

\item{n}{Number of bands. A positive integer.}
}
\description{
Divide space into tiles which contain equal amounts of stuff.
}
