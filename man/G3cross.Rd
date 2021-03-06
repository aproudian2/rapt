% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rapt_extend.R
\name{G3cross}
\alias{G3cross}
\title{Multitype Nearest Neighbour Distance Function (i-to-j)}
\usage{
G3cross(X, i, j, rmax = NULL, nrval = 128, correction = c("rs", "km", "han"))
}
\arguments{
\item{X}{The observed point pattern, from which an estimate of the cross type
distance distribution function \eqn{G[3ij](r)} will be computed. It must
be a multitype point pattern (a marked point pattern whose marks are a
factor). See Details.}

\item{i}{The type (mark value) of the points in \code{X} from which distances
are measured. A character string (or something that will be converted to a
character string). Defaults to the first level of \code{marks(X)}.}

\item{j}{The type (mark value) of the points in \code{X} to which distances
are measured. A character string (or something that will be converted to a
character string). Defaults to the second level of \code{marks(X)}.}

\item{rmax}{Optional. Maximum value of argument \emph{r} for which
\eqn{G[3ij](r)} will be estimated.}

\item{nrval}{Optional. Number of values of \emph{r} for which
\eqn{G[3ij](r)} will be estimated. A large value of \code{nrval} is
required to avoid discretisation effects.}

\item{correction}{Optional. Character string specifying the edge
correction(s) to be used. Options are \code{"none"}, \code{"rs"},
\code{"km"}, \code{"hanisch"}, and \code{"best"}. Alternatively
\code{correction="all"} selects all options.}
}
\value{
An object of class "fv" (see \code{\link[spastat]{fv.object}}).
}
\description{
For a multitype point pattern, estimate the distribution of the distance from
a point of type i to the nearest point of type j.
}
\details{
The function \code{G3cross} and its companions \code{G3dot} (unimplemented)
and \code{\link{G3multi}} are generalisations of the function
\code{\link[spatstat]{G3est}} to multitype point patterns.

A multitype point pattern is a spatial pattern of points classified into a
finite number of possible "colors" or "types." In the \strong{spatstat}
package, a multitype pattern is represented as a single point pattern object
in which the points carry marks, and the mark value attached to each point
determines the type of that point.

The argument \code{X} must be a point pattern (object of class "pp3"). It
must be a marked point pattern, and the mark vector \code{X$marks} must be a
factor. The arguments \code{i} and \code{j} will be interpreted as levels of
the factor \code{X$marks}. (\strong{Warning:} this means that an integer value
\code{i=3} will be interpreted as the number 3, \emph{not} the 3rd smallest
level).

The "cross-type" (type \emph{i} to type \emph{j}) nearest neighbour
distance distribution function of a multitype point process is the cumulative
distribution function \eqn{G[3ij](r)} of the distance from a typical random
point of the process with type \emph{i} the nearest point of type
\emph{j}.

An estimate of \eqn{G[3ij](r)} is a useful summary statistic in exploratory
data analysis of a multitype point pattern. If the process of type \emph{i}
points were independent of the process of type \emph{j} points, then
\eqn{G[3ij](r)} would equal \eqn{F[3j](r)}, the empty space function of the
type \emph{j} points. For a multitype Poisson point process where the type
\emph{i} points have intensity \eqn{\lambda[i]}, we have

\deqn{G[3ij](r) = 1 - exp( - \lambda[j] * (4/3) * pi * r^3)}

Deviations between the empirical and theoretical \eqn{G[3ij](r)} curves may
suggest dependence between the points of types \emph{i} and \emph{j}.

This algorithm estimates the distribution function \eqn{G[3ij](r)} from the
point pattern \code{X}. It assumes that \code{X} can be treated as a
realisation of a stationary (spatially homogeneous) random spatial point
process in the plane, observed through a bounded window. The window (which is
specified in \code{X} as \code{Domain(X)}) may have arbitrary shape. Biases
due to edge effects are treated in the same manner as in
\code{\link[spatstat]{G3est}}.

The argument \code{rmax} is the maximum value of the distance \emph{r} at
which \eqn{G[3ij](r)} should be evaluated. It is also used to determine (in
combination with \code{nrval}) the breakpoints (in the sense of
\code{\link[graphics]{hist}}) for the computation of histograms of distances.
The reduced-sample and Kaplan-Meier estimators are computed from histogram
counts. In the case of the Kaplan-Meier estimator this introduces a
discretisation error which is controlled by the fineness of the breakpoints.

The algorithm also returns an estimate of the hazard rate function,
\eqn{lambda(r)}, of \eqn{G[3ij](r)}. This estimate should be used with
caution as \eqn{G[3ij](r)} is not necessarily differentiable.

The naive empirical distribution of distances from each point of the pattern
\code{X} to the nearest other point of the pattern, is a biased estimate of
\eqn{G[3ij](r)}. However this is also returned by the algorithm, as it is
sometimes useful in other contexts. Care should be taken not to use the
uncorrected empirical \eqn{G[3ij](r)} as if it were an unbiased estimator of
\eqn{G[3ij](r)}.
}
\seealso{
\code{\link{G3multi}}, \code{\link[spatstat]{G3est}},
\code{\link[spatstat]{marks}}

Other spatstat extensions: 
\code{\link{G3multi}()},
\code{\link{Tstat.pp3}()},
\code{\link{bdist.points}()},
\code{\link{marktable.pp3}()},
\code{\link{marktable}()},
\code{\link{quadratcount.pp3}()},
\code{\link{quadrats.pp3}()},
\code{\link{rPoissonCluster3}()},
\code{\link{rjitter.pp3}()},
\code{\link{rjitter.ppp}()},
\code{\link{rjitter}()},
\code{\link{rpoint3}()},
\code{\link{sample.pp3}()},
\code{\link{sample.ppp}()},
\code{\link{shift.pp3}()},
\code{\link{studpermu.test.pp3}()},
\code{\link{studpermu.test}()},
\code{\link{superimpose.pp3}()}
}
\concept{spatstat extensions}
