% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rapt_spec.R
\name{rangePOS}
\alias{rangePOS}
\title{Extract Hits in Range}
\usage{
rangePOS(pos, start, end)
}
\arguments{
\item{pos}{The \code{POS} or \code{ATO} to be ranged. A \code{data.frame} created by
\code{\link{readPOS}} or \code{\link{readATO}}}

\item{start}{Numeric. The start of the mass range}

\item{end}{Numeric. The end of the mass range}
}
\value{
A \code{data.frame} of the same structure as \code{pos} containing only hits
in the provided range
}
\description{
\code{rangePOS} extracts the rows of a \code{POS} or \code{ATO} object whose
mass is within the provided range. To count the number of hits within a mass
range, use \code{\link{rangeCount}}.
}
\seealso{
Other ranging functions: 
\code{\link{rangeCount}()},
\code{\link{rngCount}()},
\code{\link{rngPOS}()}
}
\concept{ranging functions}
