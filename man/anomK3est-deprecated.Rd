% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rapt_envelopes.R
\name{anomK3est-deprecated}
\alias{anomK3est-deprecated}
\title{Perfrom anomaly K3est on a \code{\link[spatstat]{pp3}} object.}
\arguments{
\item{pattern}{The \code{\link[spatstat]{pp3}} object to analyze.}

\item{toSub}{Returned from \code{\link{pK3est}} with \code{anom = TRUE}.
Second item in the returned list. The data to subtract from the results.}

\item{rmax}{Max r value. See \code{\link[spatstat]{K3est}}. Should be the
same as the envelopes that you are comparing to.}

\item{nrval}{Number of r values. See \code{\link[spatstat]{K3est}}. Should be
the same as the envelopes that you are comparing to.}

\item{correction}{See \code{\link{pK3est}}.}
}
\value{
Returns data fram containing r values and associated anomaly K3est
values.
}
\description{
See \code{\link[spatstat]{K3est}}. Performs the anomaly K3est on a set of
point cloud data. This means taking the square root, and subtracting the 50th
percentile from the results. This centers the curve around zero, and
standardizeds the variance at different radii. Used for comparing data to
envelopes from \code{\link{pK3est}} where \code{anom = TRUE}. Will subtract
the same values used in the pK3est test that is being compared to.
}
\seealso{
\code{\link{rapt-deprecated}}
}
\keyword{internal}
