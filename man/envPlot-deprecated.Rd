% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rapt_envelopes.R
\name{envPlot-deprecated}
\alias{envPlot-deprecated}
\title{Plot envelopes of K3est test}
\arguments{
\item{tests}{The return file from \code{p(K/G/F)3est} or the first, [\link{1}],
entry in the list returned by \code{p(K/G/F)3est} with \code{anom = TRUE}.}

\item{percentiles}{Numerical vector of percentiles that you want to see the
envelopes for. Each between 0 and 1.}

\item{ylim}{Numerical vector containing the min and max values for the y axis
on the plot.}

\item{xlim}{Numerical vector containing the min and max values for the x axis
on the plot.}

\item{ylab}{Y axis label.}

\item{xlab}{X axis label.}

\item{leg}{True or false whether to show the automatically generated legend.}

\item{colors}{List of color names to make the envelopes.}

\item{...}{Arguments to be passed into \code{plot()}.}
}
\value{
Nothing, just produces a plot.
}
\description{
Plot the results of envelope calculations from the \code{\link{pK3est}},
\code{\link{pG3est}}, or \code{\link{pF3est}} functions, with the ability to
choose the percentiles for plotting.
}
\seealso{
\code{\link{rapt-deprecated}}
}
\keyword{internal}
