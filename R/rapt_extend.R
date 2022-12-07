#
# This file contains extensions to the spatstat package.
#

#### marktable ####
#' Tabulate Marks in Neighbourhood of Every Point in a Point Pattern
#'
#' @description This is an S3 generic that extends the use of
#'   \code{\link[spatstat.core]{marktable}} beyond "ppp" objects.
#'
#' @param X A marked point pattern. An object of class "ppp".
#' @param R Neighbourhood radius. Incompatible with `N`.
#' @param N Number of neighbours of each point. Incompatible with `R`.
#' @param exclude Logical. If `exclude=TRUE`, the neighbours of a point do
#'   not include the point itself. If `exclude=FALSE`, a point belongs to
#'   its own neighbourhood.
#' @param collapse Logical. If `collapse=FALSE` (the default) the results
#'   for each point are returned as separate rows of a table. If
#'   `collapse=TRUE`, the results are aggregated according to the type of
#'   point.
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat.core]{marktable.ppp}}, \code{\link{marktable.pp3}}
#' @export
marktable <- function(X, ...) UseMethod("marktable")

### marktable.ppp ###
#' Tabulate Marks in Neighbourhood of Every Point in a Point Pattern
#'
#' @seealso \code{\link[spatstat.core]{marktable}}
#' @export
marktable.ppp <- spatstat.core::marktable

### marktable.pp3 ###
#' Tabulate Marks in Neighbourhood of Every Point in a Point Pattern
#'
#' Visit each point in a point pattern, find the neighbouring points, and
#' compile a frequency table of the marks of these neighbour points.
#'
#' @param X A marked point pattern. An object of class "ppp".
#' @param R Neighbourhood radius. Incompatible with `N`.
#' @param N Number of neighbours of each point. Incompatible with `R`.
#' @param exclude Logical. If `exclude=TRUE`, the neighbours of a point do
#'   not include the point itself. If `exclude=FALSE`, a point belongs to
#'   its own neighbourhood.
#' @param collapse Logical. If `collapse=FALSE` (the default) the results
#'   for each point are returned as separate rows of a table. If
#'   `collapse=TRUE`, the results are aggregated according to the type of
#'   point.
#'
#' @family spatstat extensions
#'
#' @export
marktable.pp3 <- function(X, R, N, exclude = TRUE, collapse = FALSE) {
  spatstat.geom::verifyclass(X, "pp3")
  if (!spatstat.geom::is.marked(X, dfok = FALSE)) {
    stop("point pattern has no marks")
  }
  gotR <- !missing(R) && !is.null(R)
  gotN <- !missing(N) && !is.null(N)
  if (gotN == gotR) {
    stop("Exactly one of the arguments N and R should be given")
  }
  stopifnot(is.logical(exclude) && length(exclude) == 1)
  m <- spatstat.geom::marks(X)
  if (!is.factor(m)) {
    stop("marks must be a factor")
  }
  if (gotR) {
    stopifnot(is.numeric(R) && length(R) == 1 && R > 0)
    p <- spatstat.geom::closepairs(X, R, what = "indices")
    pi <- p$i
    pj <- p$j
    if (!exclude) {
      n <- spatstat.geom::npoints(X)
      pi <- c(pi, 1:n)
      pj <- c(pj, 1:n)
    }
  } else {
    stopifnot(is.numeric(N) && length(N) == 1)
    ii <- seq_len(spatstat.geom::npoints(X))
    nn <- spatstat.geom::nnwhich(X, k = 1:N)
    if (N == 1) {
      nn <- matrix(nn, ncol = 1)
    }
    if (!exclude) {
      nn <- cbind(ii, nn)
    }
    pi <- as.vector(row(nn))
    pj <- as.vector(nn)
  }
  if (!collapse) {
    i <- factor(pi, levels = seq_len(spatstat.geom::npoints(X)))
    mj <- m[pj]
    mat <- table(point = i, mark = mj)
  } else {
    mi <- m[pi]
    mj <- m[pj]
    mat <- table(point = mi, neighbour = mj)
  }
  return(mat)
}

#### rjitter ####
#' Random Perturbation of a Point Pattern
#'
#' @description This is an S3 generic that extends the use of
#'   \code{\link[spatstat.geom]{rjitter}} beyond "ppp" objects.
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat.geom]{rjitter}}, \code{\link{rjitter.ppp}},
#'   \code{\link{rjitter.pp3}}
#'
#' @export
rjitter <- function(X, ...) UseMethod("rjitter")

### rjitter.ppp ###
#' Random Perturbation of a Point Pattern
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat.geom]{rjitter}}, \code{\link{rjitter.pp3}}
#' @export
rjitter.ppp <- spatstat.geom::rjitter

### rjitter.pp3 ###
#' Random Perturbation of a Point Pattern
#'
#' Applies independent random displacements to each point in a point pattern.
#' Extends \code{\link[spatstat.geom]{rjitter}} to \code{\link[spatstat.geom]{pp3}}.
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat.geom]{rjitter}}, \code{\link{rjitter.ppp}}
#'
#' @export
rjitter.pp3 <- function(X, radius, retry = TRUE, giveup = 10000, ...,
                        nsim = 1, drop = TRUE) {
  verifyclass(X, "pp3")
  if (!missing(nsim)) {
    check.1.integer(nsim)
    stopifnot(nsim >= 0)
  }
  nX <- npoints(X)
  W <- domain(X)
  if (nX == 0) {
    result <- rep(list(X), nsim)
    result <- simulationresult(result, nsim, drop)
    return(result)
  }
  if (missing(radius) || is.null(radius)) {
    ## Stoyan rule of thumb
    bws <- 0.15 / sqrt(5 * nX / volume(W))
    radius <- min(bws, shortside(W)) # had to get rid of "Frame" function call
    sameradius <- TRUE
  } else {
    ## either one radius, or a vector of radii
    check.nvector(radius, nX, oneok = TRUE, vname = "radius")
    sameradius <- (length(radius) == 1)
  }
  #'
  result <- vector(mode = "list", length = nsim)
  for (isim in seq_len(nsim)) {
    if (!retry) {
      ## points outside window are lost
      rD <- radius * sqrt(runif(nX))
      thetaD <- runif(nX, max = pi)
      phiD <- runif(nX, max = 2 * pi)
      xnew <- X$data$x + rD * sin(thetaD) * cos(phiD)
      ynew <- X$data$y + rD * sin(thetaD) * sin(phiD)
      znew <- X$data$z + rD * cos(thetaD)
      ok <- inside.boxx(xnew, ynew, znew, w = W)
      result[[isim]] <- pp3(xnew[ok], ynew[ok], znew[ok], W)
    } else {
      ## retry = TRUE: condition on points being inside window
      undone <- rep.int(TRUE, nX)
      triesleft <- giveup
      Xshift <- X
      while (any(undone)) {
        triesleft <- triesleft - 1
        if (triesleft <= 0) {
          break
        }
        Y <- Xshift[undone]
        nY <- npoints(Y)
        RY <- if (sameradius) radius else radius[undone]
        rD <- RY * sqrt(runif(nY))
        thetaD <- runif(nY, max = pi)
        phiD <- runif(nY, max = 2 * pi)

        ## if there is only one point, then data structure is screwed up
        if (nY == 1) {
          xnew <- Y$data$x[[1]] + rD * sin(thetaD) * cos(phiD)
          ynew <- Y$data$y[[1]] + rD * sin(thetaD) * sin(phiD)
          znew <- Y$data$z[[1]] + rD * cos(thetaD)
          ok <- inside.boxx(xnew, ynew, znew, w = W)
        } else {
          xnew <- Y$data$x + rD * sin(thetaD) * cos(phiD)
          ynew <- Y$data$y + rD * sin(thetaD) * sin(phiD)
          znew <- Y$data$z + rD * cos(thetaD)
          ok <- inside.boxx(xnew, ynew, znew, w = W)
        }

        # result[[isim]] <- pp3(xnew[ok], ynew[ok], znew[ok], W)



        if (any(ok)) {
          changed <- which(undone)[ok]
          Xshift$data$x[changed] <- xnew[ok]
          Xshift$data$y[changed] <- ynew[ok]
          Xshift$data$z[changed] <- znew[ok]
          undone[changed] <- FALSE
        }
      }
      result[[isim]] <- Xshift
    }
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}

#### superimpose.pp3 ####
#' Superimpose Several Geometric Patterns
#'
#' Superimposes any number of 3D point patterns
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat.geom]{superimpose}}
#'
#' @export
superimpose.pp3 <- function(..., W = NULL, check = FALSE) {
  input.list <- list(...)
  df.list <- lapply(input.list, as.data.frame)
  df.comb <- do.call(rbind, df.list)
  out.pp3 <- createSpat(df.comb, win = W)
  if (!is.null(df.comb$marks)) {
    spatstat.geom::marks(out.pp3) <- df.comb$marks
  }
  return(out.pp3)
}

#### shift.pp3 ####
#' Apply Vector Translation
#'
#' Applies a vector shift to a 3D point pattern
#'
#' @param X Point pattern(object of class"pp3")
#' @param vec Numeric. A vector of length 3 specifying the translation
#' @param origin Location that will be shifted to the origin. Either a vector of
#'   length 3 specifying the new origin or one of the character strings
#'   "midpoint" or "bottomleft" (will be extended)
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat.geom]{shift}}
#'
#' @export
shift.pp3 <- function(X, vec = c(0, 0, 0), ..., origin = NULL) {
  spatstat.geom::verifyclass(X, "pp3")
  if (!is.null(origin)) {
    if (!missing(vec)) {
      warning("argument vec ignored; overruled by argument origin")
    }
    if (is.numeric(origin) & length(origin) == 3) {
      locn <- origin
    } else if (is.character(origin)) {
      origin <- spatstat.geom::pickoption(
        "origin", origin,
        c(midpoint = "midpoint", bottomleft = "bottomleft")
      )
      W <- X$domain
      locn <- switch(origin,
        midpoint = {
          c(mean(W$domain$xrange), mean(W$domain$yrange), mean(W$domain$zrange))
        },
        bottomleft = {
          c(W$domain$xrange[1], W$domain$yrange[1], W$domain$zrange[1])
        }
      )
    } else {
      stop("origin must be a character string or a numeric vector")
    }
    vec <- -locn
  }
  Y <- spatstat.geom::pp3(X$data$x + vec[1], X$data$y + vec[2], X$data$z + vec[3],
    xrange = X$domain$xrange + vec[1],
    yrange = X$domain$yrange + vec[2],
    zrange = X$domain$zrange + vec[3],
    marks = marks(X)
  )
  attr(Y, "lastshift") <- vec
  return(Y)
}

#### sample ####
### sample.ppp ###
#' Sample a Planar Point Pattern
#'
#' @param X An object of class "ppp". The point pattern from which to sample.
#' @param size A numeric. The number of points to sample.
#' @return A "ppp". The sampled point pattern.
#'
#' @family spatstat extensions
#' @seealso \code{\link{sample.pp3}}, \code{\link[base]{sample}}
#'
#' @export
sample.ppp <- function(X, size) {
  sam.n <- seq_len(spatstat.geom::npoints(X))
  sam.pts <- sample(sam.n, size, replace = FALSE)
  sam.dat <- X[sam.pts]
  return(sam.dat)
}

### sample.pp3 ###
#' Sample A 3D Point Pattern
#'
#' @param X An object of class "pp3". The point pattern from which to sample.
#' @param size A numeric. The number of points to sample.
#' @return A "pp3". The sampled point pattern.
#'
#' @family spatstat extensions
#' @seealso \code{\link{sample.ppp}}, \code{\link[base]{sample}}
#'
#' @export
sample.pp3 <- function(X, size) {
  sam.n <- seq_len(spatstat.geom::npoints(X))
  sam.pts <- sample(sam.n, size, replace = FALSE)
  sam.dat <- X[sam.pts]
  return(sam.dat)
}

#### findClusters.pp3 ####
#' Find Clusters by NN Adjacency Marks
findClusters.pp3 <- function(X, mark, k = 1) {
  if (!(mark %in% marks(X))) {
    stop("The specified mark does not exist in the pattern")
  } else if (length(mark) != 1) {
    stop("Only one mark may be tested")
  }
  spl <- split(X)[[mark]]
  # required because pp3 functions don't work on marked pp3 by default
  class(spl) <- c("pp3", "ppx")
  ind <- rownames.pp3(spl)
  fCl <- lapply(ind, function(x) {
    clus <- x
    rem <- rownames.pp3(X)
    nn <- 0
    while (length(nn) > 0) {
      rem <- setdiff(rem, nn)
      nn <- nncross(spl[clus], X[rem], k = 1:k, what = "which")
      nn <- unlist(nn)
      nn <- rownames.pp3(X[rem])[nn]
      nn <- nn[marks(X[nn]) == mark]
      clus <- append(clus, nn)
    }
    clus <- sort(unique(clus))
    return(clus)
  })
  return(fCl)
}

#### intensity.pp3 ####
#' Extends \code{\link[spatstat.geom]{intensity}} to \code{\link[spatstat.geom]{pp3}}.
#'
#' @name intensity.pp3-deprecated
#' @seealso \code{\link{rapt-deprecated}}
#' @keywords internal
NULL
#' @rdname rapt-deprecated
#' @section \code{intensity.pp3}:
#'   For \code{intensity.pp3}, use \code{\link[spatstat.geom]{intensity.ppx}}
#'
#' @export
intensity.pp3 <- function(X, weights = NULL) {
  n <- spatstat.geom::npoints(X)
  a <- spatstat.geom::volume(domain(X))
  if (is.null(weights)) {
    if (spatstat.geom::is.multitype(X)) {
      mks <- spatstat.geom::marks(X)
      answer <- as.vector(table(mks)) / a
      names(answer) <- levels(mks)
    } else {
      answer <- n / a
    }
    return(answer)
  } else {
    warning("Weights is not yet implemented for pp3")
  }
}

#### rownames.pp3 ####
#' Extends \code{\link[base:row+colnames]{rownames}} to
#' \code{\link[spatstat.geom]{pp3}}.
#'
#' @param pat An object of class "pp3". The point pattern from which to extract
#'   rownames.
#' @return Character vector. The rownames of the point pattern.
#' @seealso \code{\link[base:row+colnames]{rownames}}
rownames.pp3 <- function(pat) {
  dat <- rownames(as.data.frame(pat))
  return(dat)
}
#### rownames.pp3<- ####
#' Extends \code{\link[base:row+colnames]{rownames}} to
#' \code{\link[spatstat.geom]{pp3}}.
#'
#' @param pat An object of class "pp3". The point pattern to assign rownames.
#' @param lab character. The new rownames.
#' @seealso \code{\link[base:row+colnames]{rownames}}
"rownames.pp3<-" <- function(pat, lab) {
  rownames(as.data.frame(pat)) <- lab
  return(pat)
}

#### plot3d.pp3 ####
#' Plot a \code{\link[spatstat.geom]{pp3}} in a manipulatable 3D plot.
#'
#' (requires the rgl library)
#' @param X An object of class "pp3". The point pattern to visualize
#' @param ... Other arguments to pass to \code{\link[rgl]{plot3d}} from the
#'   `rgl` library.
#'
#' @family visualization functions
#' @seealso \code{\link[rgl]{plot3d}}
#'
#' @export
plot3d.pp3 <- function(X, ...) {
  rgl::plot3d(coords(X), ...)
}

#### quadrats.pp3 ####
#' Divide a 3D Point Pattern into Sub-Volumes
#'
#' Divides volume into quadrats and returns them.
#'
#' @param X The \code{\link[spatstat.geom]{pp3}} object to split up.
#' @param nx,ny,nz Number of rectangular quadrats in the x, y, and z directions,
#'   if you wish to split up your point patthern by number of boxes.
#' @param box.dims Vector containing the dimensions of the subsetted 3D boxxes,
#'   if you wish to define the individual boxx size. Use either `nx`, `ny`,
#'   `nz` or `box.dims`, but not both. See Details.
#'
#' @return A list containing the split up "pp3" objects.
#'
#' @details
#' If `box.dims` is not commensurate with the size of the object, the
#' quadrats are justified toward the smallest boundary in each dimension
#' (*i.e.* (`c(min(x), min(y), min(z))`)).
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat.geom]{quadrats}}
#'
#' @export
# Should return an object of class "tess" to provide consistent behaviour
quadrats.pp3 <- function(X, nx = 5, ny = nx, nz = ny, box.dims = NULL) {
  # Rewrite to use n = c(nx,ny,nz) with a default of c(1,1,1).
  spatstat.geom::verifyclass(X, c("pp3", "box3"))
  w <- as.box3(X)
  xlim <- w$xrange
  ylim <- w$yrange
  zlim <- w$zrange

  # create box3objects for each quadrat
  if (is.null(box.dims)) {
    xbreaks <- seq(xlim[1], xlim[2], length.out = (nx + 1))
    ybreaks <- seq(ylim[1], ylim[2], length.out = (ny + 1))
    zbreaks <- seq(zlim[1], zlim[2], length.out = (nz + 1))

    ntot <- nx * ny * nz
  } else {
    xbreaks <- seq(xlim[1], xlim[2], by = box.dims[1])
    ybreaks <- seq(ylim[1], ylim[2], by = box.dims[2])
    zbreaks <- seq(zlim[1], zlim[2], by = box.dims[3])
  }


  gridvals <- list()
  cnt <- 1

  for (i in 1:(length(xbreaks) - 1)) {
    for (j in 1:(length(ybreaks) - 1)) {
      for (k in 1:(length(zbreaks) - 1)) {
        gridvals[[cnt]] <- box3(
          xrange = xbreaks[i:(i + 1)],
          yrange = ybreaks[j:(j + 1)],
          zrange = zbreaks[k:(k + 1)]
        )
        cnt <- cnt + 1
      }
    }
  }

  # browser()
  coo <- coords(X)
  if (!is.null(marks(X))) {
    marks <- marks(X)
  }

  boxes <- lapply(gridvals, function(x) {
    res.coo <- coo[inside.boxx(X, w = x), ]
    if (!is.null(marks(X))) {
      res.marks <- marks[inside.boxx(X, w = x)]
      return(pp3(res.coo$x, res.coo$y, res.coo$z, x, marks = res.marks))
    } else {
      return(pp3(res.coo$x, res.coo$y, res.coo$z, x))
    }
  })

  return(boxes)
}

#### quadratcount.pp3 ####
#' Count Points in Sub-Volumes of a 3D Point Pattern
#'
#' Divides volume into quadrats and counts the number of points in each quadrat.
#'
#' @param X The \code{\link[spatstat.geom]{pp3}} object to split up.
#' @param nx,ny,nz Number of ractangular quadrats in the x, y, and z directions.
#'
#' @return A `data.frame` object containing the number of counts in each
#'   quadrat.
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat.geom]{quadratcount}}
#'
#' @export
quadratcount.pp3 <- function(X, nx = 5, ny = 5, nz = 5) {
  spatstat.geom::verifyclass(X, "pp3")
  w <- domain(X)

  # create box3objects for each quadrat
  xlim <- w$xrange
  ylim <- w$yrange
  zlim <- w$zrange

  xbreaks <- seq(xlim[1], xlim[2], length.out = (nx + 1))
  ybreaks <- seq(ylim[1], ylim[2], length.out = (ny + 1))
  zbreaks <- seq(zlim[1], zlim[2], length.out = (nz + 1))

  ntot <- nx * ny * nz
  gridvals <- list()
  cnt <- 1

  for (i in 1:nx) {
    for (j in 1:ny) {
      for (k in 1:nz) {
        gridvals[[cnt]] <- box3(
          xrange = xbreaks[i:(i + 1)],
          yrange = ybreaks[j:(j + 1)],
          zrange = zbreaks[k:(k + 1)]
        )
        cnt <- cnt + 1
      }
    }
  }

  inside.tf <- lapply(gridvals, function(x) {
    inside.boxx(X, w = x)
  })
  counts <- lapply(inside.tf, function(x) {
    sum(x)
  })
  counts <- unlist(counts)
  return(data.frame(quad.no = seq(1, ntot), count = counts))
}

#### quantess.pp3 ####
#' Quantile Tessellation
#'
#' Divide space into tiles which contain equal amounts of stuff.
#'
#' @param M A pp3
#' @param Z A spatial covariate (a pixel image or a function(x,y,z)) or one of
#'   the strings `"x"`, `"y"`, or `"z"` indicating the Cartesian coordinates
#'   *x*, *y*, or *z* or one of the strings `"rad"` or `"ang"` indicating polar
#'   coordinates. The range of values of `Z` will be broken into n bands
#'   containing equal amounts of stuff.
#' @param n Number of bands. A positive integer.
#'
quantess.pp3 <- function(M, Z, n, ..., type = 2, origin = c(0, 0), eps = NULL) {}

#### G3multi ####
#' Marked Nearest Neighbour Distance Function
#'
#' For a marked point pattern, estimate the distribution of the distance from a
#' typical point in subset `I` to the nearest point of subset `J`.
#'
#' @param X The observed point pattern, from which an estimate of the multitype
#'   distance distribution function \eqn{G[3IJ](r)} will be computed. It must be
#'   a marked point pattern. See Details.
#' @param I Subset of points of `X` from which distances are measured.
#' @param J Subset of points in `X` to which distances are measured.
#' @param rmax Optional. Maximum value of argument *r* for which \eqn{G[3IJ](r)}
#'   will be estimated.
#' @param nrval Optional. Number of values of *r* for which \eqn{G[3IJ](r)} will
#'   be estimated. A large value of `nrval` is required to avoid discretisation
#'   effects.
#' @param disjoint Optional flag indicating whether the subsets `I` and `J` are
#'   disjoint. If missing, this value will be computed by inspecting the vectors
#'   `I` and `J`.
#' @param correction Optional. Character string specifying the edge
#'   correction(s) to be used. Options are `"none"`, `"rs"`, `"km"`,
#'   `"hanisch"`, and `"best"`. Alternatively `correction="all"` selects all
#'   options.
#'
#' @details
#'
#' The function `G3multi` generalises \code{\link[spatstat.core]{G3est}} (for
#' unmarked point patterns) and \code{G3dot} (unimplmented) and
#' \code{\link{G3cross}} (for multitype point patterns) to arbitrary marked
#' point patterns.
#'
#' @family spatstat extensions
#'
#' @seealso \code{\link{G3cross}}, \code{\link[spatstat.core]{G3est}}
#'
#' @export
G3multi <- function(X, I, J, rmax = NULL, nrval = 128, disjoint = NULL,
                    correction = c("rs", "km", "han")) {
  spatstat.geom::verifyclass(X, "pp3")
  W <- X$domain
  npts <- spatstat.geom::npoints(X)
  volW <- spatstat.geom::volume(W)
  if (is.null(rmax)) {
    rmax <- spatstat.geom::diameter(W) / 2
  }
  I <- spatstat.geom::ppsubset(X, I)
  J <- spatstat.geom::ppsubset(X, J)
  if (is.null(I) || is.null(J)) {
    stop("I and J must be valid subset indices")
  }
  nI <- sum(I)
  nJ <- sum(J)
  if (nI == 0) {
    stop("No points satisfy condition I")
  }
  if (nJ == 0) {
    stop("No points satisfy condition J")
  }
  if (is.null(disjoint)) {
    disjoint <- !any(I & J)
  }
  if (is.null(correction)) {
    correction <- c("rs", "km", "han")
  }
  correction <- pickoption("correction", correction,
    c(
      none = "none",
      border = "rs", rs = "rs",
      KM = "km", km = "km", Kaplan = "km",
      han = "han", Hanisch = "han", hanisch = "han",
      best = "km"
    ),
    multi = TRUE
  )
  lamJ <- nJ / volW

  XI <- X[I]
  XJ <- X[J]
  if (disjoint) {
    nnd <- spatstat.geom::nncross(XI, XJ, what = "dist")
  } else {
    seqnp <- seq_len(npts)
    iX <- seqnp[I]
    iY <- seqnp[J]
    nnd <- spatstat.geom::nncross(XI, XJ, iX, iY, what = "dist")
  }
  bdry <- bdist.points(XI)
  d <- (nnd <= bdry)
  r <- seq(0, rmax, length.out = nrval)
  breaks <- c(r[1L] - r[2L], r)
  df <- data.frame(r = r, theo = 1 - exp(-lamJ * (4 / 3) * pi * r^3))
  fname <- c("G", "list(I,J)")
  Z <- fv(df, "r", quote(G[I, J](r)), "theo", . ~ r, c(0, rmax),
    c("r", spatstat.core::makefvlabel(NULL, NULL, fname, "pois")),
    c("distance argument r", "theoretical Poisson %s"),
    fname = fname, yexp = quote(G[list(I, J)](r))
  )
  if ("none" %in% correction) {
    if (npts == 0) {
      edf <- zeroes
    } else {
      hh <- hist(nnd[nnd <= rmax], breaks = breaks, plot = FALSE)$counts
      edf <- cumsum(hh) / length(nnd)
    }
    Z <- bind.fv(
      Z, data.frame(raw = edf),
      spatstat.core::makefvlabel(NULL, "hat", fname, "raw"),
      "uncorrected estimate of %s", "raw"
    )
  }
  if ("han" %in% correction) {
    if (npts == 0) {
      G <- zeroes
    } else {
      x <- nnd[d]
      a <- spatstat.geom::eroded.volumes(W, r)
      h <- hist(x[x <= rmax], breaks = breaks, plot = FALSE)$counts
      G <- cumsum(h / a)
      G <- G / max(G[is.finite(G)])
    }
    Z <- bind.fv(
      Z, data.frame(han = G),
      spatstat.core::makefvlabel(NULL, "hat", fname, "han"),
      "Hanisch estimate of %s", "han"
    )
    attr(Z, "alim") <- range(r[G <= 0.9])
  }
  if (any(correction %in% c("rs", "km"))) {
    if (npts == 0) {
      result <- data.frame(rs = zeroes, km = zeroes, hazard = zeroes)
    } else {
      o <- pmin.int(nnd, bdry)
      result <- spatstat.core::km.rs(o, bdry, d, breaks)
      result <- as.data.frame(result[c("rs", "km", "hazard")])
    }
    Z <- bind.fv(
      Z, result, c(
        spatstat.core::makefvlabel(NULL, "hat", fname, "bord"),
        spatstat.core::makefvlabel(NULL, "hat", fname, "km"), "hazard(r)"
      ),
      c(
        "border corrected estimate of %s", "Kaplan-Meier estimate of %s",
        "Kaplan-Meier estimate of hazard function lambda(r)"
      ),
      "km"
    )
    attr(Z, "alim") <- range(r[result$km <= 0.9])
  }
  nama <- names(Z)
  spatstat.core::fvnames(Z, ".") <- rev(nama[!(nama %in% c("r", "hazard"))])
  formula(Z) <- . ~ r
  spatstat.geom::unitname(Z) <- spatstat.geom::unitname(X)
  return(Z)
}

#### G3cross ####
#' Multitype Nearest Neighbour Distance Function (i-to-j)
#'
#' For a multitype point pattern, estimate the distribution of the distance from
#' a point of type i to the nearest point of type j.
#'
#' @param X The observed point pattern, from which an estimate of the cross type
#'   distance distribution function \eqn{G[3ij](r)} will be computed. It must
#'   be a multitype point pattern (a marked point pattern whose marks are a
#'   factor). See Details.
#' @param i The type (mark value) of the points in `X` from which distances
#'   are measured. A character string (or something that will be converted to a
#'   character string). Defaults to the first level of `marks(X)`.
#' @param j The type (mark value) of the points in `X` to which distances
#'   are measured. A character string (or something that will be converted to a
#'   character string). Defaults to the second level of `marks(X)`.
#' @param rmax Optional. Maximum value of argument *r* for which
#'   \eqn{G[3ij](r)} will be estimated.
#' @param nrval Optional. Number of values of *r* for which
#'   \eqn{G[3ij](r)} will be estimated. A large value of `nrval` is
#'   required to avoid discretisation effects.
#' @param correction Optional. Character string specifying the edge
#'   correction(s) to be used. Options are `"none"`, `"rs"`,
#'   `"km"`, `"hanisch"`, and `"best"`. Alternatively
#'   `correction="all"` selects all options.
#'
#' @return An object of class "fv" (see \code{\link[spastat]{fv.object}}).
#'
#' @family spatstat extensions
#'
#' @details
#' The function `G3cross` and its companions \code{G3dot} (unimplemented)
#' and \code{\link{G3multi}} are generalisations of the function
#' \code{\link[spatstat.core]{G3est}} to multitype point patterns.
#'
#' A multitype point pattern is a spatial pattern of points classified into a
#' finite number of possible "colors" or "types." In the **spatstat**
#' package, a multitype pattern is represented as a single point pattern object
#' in which the points carry marks, and the mark value attached to each point
#' determines the type of that point.
#'
#' The argument `X` must be a point pattern (object of class "pp3"). It
#' must be a marked point pattern, and the mark vector `X$marks` must be a
#' factor. The arguments `i` and `j` will be interpreted as levels of
#' the factor `X$marks`. (**Warning:** this means that an integer value
#' `i=3` will be interpreted as the number 3, *not* the 3rd smallest
#' level).
#'
#' The "cross-type" (type *i* to type *j*) nearest neighbour
#' distance distribution function of a multitype point process is the cumulative
#' distribution function \eqn{G[3ij](r)} of the distance from a typical random
#' point of the process with type *i* the nearest point of type
#' *j*.
#'
#' An estimate of \eqn{G[3ij](r)} is a useful summary statistic in exploratory
#' data analysis of a multitype point pattern. If the process of type *i*
#' points were independent of the process of type *j* points, then
#' \eqn{G[3ij](r)} would equal \eqn{F[3j](r)}, the empty space function of the
#' type *j* points. For a multitype Poisson point process where the type
#' *i* points have intensity \eqn{\lambda[i]}, we have
#'
#' \deqn{G[3ij](r) = 1 - exp( - \lambda[j] * (4/3) * pi * r^3)}
#'
#' Deviations between the empirical and theoretical \eqn{G[3ij](r)} curves may
#' suggest dependence between the points of types *i* and *j*.
#'
#' This algorithm estimates the distribution function \eqn{G[3ij](r)} from the
#' point pattern `X`. It assumes that `X` can be treated as a
#' realisation of a stationary (spatially homogeneous) random spatial point
#' process in the plane, observed through a bounded window. The window (which is
#' specified in `X` as `Domain(X)`) may have arbitrary shape. Biases
#' due to edge effects are treated in the same manner as in
#' \code{\link[spatstat.core]{G3est}}.
#'
#' The argument `rmax` is the maximum value of the distance *r* at
#' which \eqn{G[3ij](r)} should be evaluated. It is also used to determine (in
#' combination with `nrval`) the breakpoints (in the sense of
#' \code{\link[graphics]{hist}}) for the computation of histograms of distances.
#' The reduced-sample and Kaplan-Meier estimators are computed from histogram
#' counts. In the case of the Kaplan-Meier estimator this introduces a
#' discretisation error which is controlled by the fineness of the breakpoints.
#'
#' The algorithm also returns an estimate of the hazard rate function,
#' \eqn{lambda(r)}, of \eqn{G[3ij](r)}. This estimate should be used with
#' caution as \eqn{G[3ij](r)} is not necessarily differentiable.
#'
#' The naive empirical distribution of distances from each point of the pattern
#' `X` to the nearest other point of the pattern, is a biased estimate of
#' \eqn{G[3ij](r)}. However this is also returned by the algorithm, as it is
#' sometimes useful in other contexts. Care should be taken not to use the
#' uncorrected empirical \eqn{G[3ij](r)} as if it were an unbiased estimator of
#' \eqn{G[3ij](r)}.
#'
#' @seealso \code{\link{G3multi}}, \code{\link[spatstat.core]{G3est}},
#'   \code{\link[spatstat.geom]{marks}}
#'
#' @export
G3cross <- function(X, i, j, rmax = NULL, nrval = 128,
                    correction = c("rs", "km", "han")) {
  spatstat.geom::verifyclass(X, "pp3")
  if (!spatstat.geom::is.marked(X, dfok = FALSE)) {
    stop(paste("point pattern has no", sQuote("marks")))
  }
  stopifnot(spatstat.geom::is.multitype(X))
  marx <- spatstat.geom::marks(X, dfok = FALSE)
  if (missing(i)) {
    i <- levels(marx)[1]
  }
  if (missing(j)) {
    j <- levels(marx)[2]
  }
  I <- (marx == i)
  if (sum(I) == 0) {
    stop("No points are of type i")
  }
  if (i == j) {
    result <- spatstat.core::G3est(X[I],
      rmax = rmax, nrval = nrval,
      correction = correction
    )
  } else {
    J <- (marx == j)
    if (sum(J) == 0) {
      stop("No points are of type j")
    }
    result <- G3multi(X, I, J,
      rmax = rmax, nrval = nrval, disjoint = FALSE,
      correction = correction
    )
  }
  iname <- spatstat.utils::make.parseable(paste(i))
  jname <- spatstat.utils::make.parseable(paste(j))
  result <- spatstat.core::rebadge.fv(result,
    substitute(G[i, j](r), list(i = iname, j = jname)),
    c("G", paste0("list(", iname, ",", jname, ")")),
    new.yexp = substitute(
      G[list(i, j)](r),
      list(i = iname, j = jname)
    )
  )
  return(result)
}

#### K3scaled ####
#' Locally Scaled K-function in Three Dimensions
#'
#' Estimates the locally-rescaled K-function of a 3D point pattern
#'
#' @param X A pp3
#' @param lambda An estimate of the intensity function f(x,y,z). Defaults to a
#'   uniform intensity, estimated with \code{\link[spatstat.geom]{intensity}}; this
#'   results in the same behavior as \code{\link[spatstat.core]{K3est}}.
#' @param correction The correction to be used. Currently only the uncorrected
#'   estimate (`correction = "none"`) is implemented.
#'
#' In 3D the cube root is used instead of the square root to estimate the
#' correction of the scaled intensity. Thas is, the eucledian distance between
#' two points \eqn{u[i]} and \eqn{u[j]}, \eqn{d[ij]}, is
#' multiplied by the mean of the cube roots of the local intensities
#' \eqn{lambda[c](u[i])} and \eqn{lambda[c](u[j])}.
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat.core]{Kscaled}}
#'
K3scaled <- function(X, lambda = NULL, ...,
                     rmax = NULL, nrval = 128,
                     correction = "none", # c("isotropic", "translate"),
                     renormalise = FALSE,
                     normpower = 1) # , sigma = NULL, varcov = NULL)
{
  spatstat.geom::verifyclass(X, "pp3")
  B <- X$domain
  npts <- spatstat.geom::npoints(X)
  volW <- spatstat.geom::volume(B)
  halfdiameter <- spatstat.geom::diameter(B) / 2
  correction.given <- !missing(correction) && !is.null(correction)
  correction <- spatstat.geom::pickoption("correction", correction,
    c(
      none = "none",
      border = "border",
      isotropic = "isotropic",
      Ripley = "isotropic",
      trans = "translate",
      translate = "translate",
      translation = "translate",
      best = "best"
    ),
    multi = TRUE
  )
  if (missing(lambda)) {
    lambda <- rep(spatstat.geom::intensity.ppx(X), npts)
  } else {
    if (is.im(lambda)) {
      # lambda <- safelookup(lambda, X)
      stop("No im3 implemented for lambda; supply a function")
    } else if (is.function(lambda)) {
      lambda <- lambda(X$data$x, X$data$y, X$data$z)
    } else if (is.ppm(lambda)) {
      # lambda <- safelookup(predict(lambda, type = "trend"),
      #                      X)
      stop("No p3m implemented for lambda; supply a function")
    } else if (!is.numeric(lambda) || !is.null(dim(lambda))) {
      # stop(paste(sQuote("lambda"),
      #            "should be a vector, a pixel image, a function or a ppm"))
      stop(paste(
        sQuote("lambda"),
        "should be a function"
      ))
    }
    check.nvector(lambda, npts)
  }
  # need to confirm that this renormalization is correct in 3D...
  if (renormalise) {
    spatstat.utils::check.1.real(normpower)
    stopifnot(normpower %in% 1:2)
    renorm.factor <- (volW / sum(1 / lambda))^(normpower / 2)
    lambda <- lambda / renorm.factor
  }
  cra <- (range(lambda))^(1 / 3)
  minrescale <- cra[1]
  maxrescale <- cra[2]
  absrmaxdefault <- halfdiameter / maxrescale
  if (is.null(rmax)) {
    absrmax <- absrmaxdefault
  } else {
    absrmax <- rmax / maxrescale
  }
  absr <- seq(0, to = absrmax, length.out = nrval)
  r <- absr * maxrescale
  breaks <- breakpts.from.r(r)$val
  alim <- c(0, min(rmax, maxrescale * absrmaxdefault))
  rthresh <- minrescale * halfdiameter
  maxabsdist <- min(rmax / minrescale, halfdiameter)
  K <- data.frame(r = r, theo = (4 / 3) * pi * r^3)
  desc <- c("distance argument r", "theoretical Poisson %s")
  fname <- c("K", "list(3,scaled)")
  yexp <- quote(K[list(3, scaled)](r))
  K <- fv(K, "r", quote(K[3, scaled](r)), "theo",
    alim = alim,
    labl = c("r", "{%s[%s]^{pois}}(r)"),
    desc = desc, fname = fname, yexp = yexp
  )
  needXI <- any(correction %in% c("translate", "isotropic"))
  close <- closepairs(X, maxabsdist, what = if (needXI) {
    "all"
  } else {
    "ijd"
  })
  I <- close$i
  J <- close$j
  cbrtLambda <- (lambda)^(1 / 3)
  lamIJ <- (cbrtLambda[I] + cbrtLambda[J]) / 2
  absDIJ <- close$d
  DIJ <- absDIJ * lamIJ
  XI <- if (needXI) {
    pp3(close$xi, close$yi, close$zi, B)
  } else {
    NULL
  }
  if (any(correction == "none")) {
    wh <- whist(DIJ, breaks)
    Kun <- cumsum(wh) / npts
    K <- bind.fv(
      K, data.frame(un = Kun), "{hat(%s)[%s]^{un}}(r)",
      "uncorrected estimate of %s", "un"
    )
  }
  ### ***** corrections not implemented *****
  # if (any(correction == "border")) {
  #   b <- bdist.points(X) * cbrtLambda
  #   bI <- b[I]
  #   RS <- Kount(DIJ, bI, b, breaks)
  #   Kb <- RS$numerator/RS$denom.count
  #   Kb[r > rthresh] <- NA
  #   K <- bind.fv(K, data.frame(border = Kb), "{hat(%s)[%s]^{bord}}(r)",
  #                "border-corrected estimate of %s", "border")
  # }
  # if (any(correction == "translate")) {
  #   XJ <- pp3(close$xj, close$yj, close$zj, B)
  #   edgewt <- edge.Trans(XI, XJ, paired = TRUE)
  #   wh <- whist(DIJ, breaks, edgewt)
  #   Ktrans <- cumsum(wh)/npts
  #   Ktrans[r >= rthresh] <- NA
  #   K <- bind.fv(K, data.frame(trans = Ktrans), "{hat(%s)[%s]^{trans}}(r)",
  #                "translation-corrected estimate of %s", "trans")
  # }
  # if (any(correction == "isotropic")) {
  #   edgewt <- edge.Ripley(XI, matrix(absDIJ, ncol = 1))
  #   wh <- whist(DIJ, breaks$val, edgewt)
  #   Kiso <- cumsum(wh)/npts
  #   Kiso[r >= rthresh] <- NA
  #   K <- bind.fv(K, data.frame(iso = Kiso), "{hat(%s)[%s]^{iso}}(r)",
  #                "Ripley isotropic correction estimate of %s", "iso")
  # }
  formula(K) <- . ~ r
  nama <- rev(colnames(K))
  fvnames(K, ".") <- nama[!(nama %in% c("r", "bord", "iso", "trans"))]
  unitname(K) <- c("normalised unit", "normalised units")
  return(K)
}

#### K3multi ####
# barely works... Needs corrections and inferface streamlining
K3multi <- function(X, I, J, r, breaks,
                    correction = c("none", "isotropic", "translation"),
                    ..., ratio = FALSE) {
  spatstat.geom::verifyclass(X, "pp3")
  npts <- npoints(X)
  W <- X$domain
  volW <- volume(W)
  I <- ppsubset(X, I)
  J <- ppsubset(X, J)
  if (is.null(I) || is.null(J)) {
    stop("I and J must be valid subset indices")
  }
  if (!any(I)) {
    stop("no points belong to subset I")
  }
  if (!any(J)) {
    stop("no points belong to subset J")
  }
  nI <- sum(I)
  nJ <- sum(J)
  lambdaI <- nI / volW
  lambdaJ <- nJ / volW
  rmax <- max(r)
  alim <- c(0, rmax)
  K <- data.frame(r = r, theo = (4 / 3) * pi * r^3)
  desc <- c("distance argument r", "theoretical Poisson %s")
  K <- fv(K, "r", quote(K[IJ](r)), "theo", , alim, c("r", "{%s[%s]^{pois}}(r)"),
    desc,
    fname = c("K", "list(I,J)"), yexp = quote(K[list(I, J)](r))
  )
  if (ratio) {
    denom <- lambdaI * lambdaJ * volW
    numK <- eval.fv(denom * K)
    denK <- eval.fv(denom + K * 0)
    attributes(numK) <- attributes(denK) <- attributes(K)
    attr(numK, "desc")[2] <- "numerator for theoretical Poisson %s"
    attr(denK, "desc")[2] <- "denominator for theoretical Poisson %s"
  }
  XI <- X[I]
  XJ <- X[J]
  close <- crosspairs(XI, XJ, max(r))
  orig <- seq_len(npts)
  imap <- orig[I]
  jmap <- orig[J]
  iX <- imap[close$i]
  jX <- jmap[close$j]
  if (any(I & J)) {
    ok <- (iX != jX)
    if (!all(ok)) {
      close$i <- close$i[ok]
      close$j <- close$j[ok]
      close$d <- close$d[ok]
    }
  }
  dcloseIJ <- close$d
  icloseI <- close$i
  jcloseJ <- close$j
  if (any(correction == "none")) {
    wh <- whist(dcloseIJ, breaks)
    numKun <- cumsum(wh)
    denKun <- lambdaI * lambdaJ * volW
    Kun <- numKun / denKun
    K <- bind.fv(
      K, data.frame(un = Kun), "{hat(%s)[%s]^{un}}(r)",
      "uncorrected estimate of %s", "un"
    )
    if (ratio) {
      numK <- bind.fv(
        numK, data.frame(un = numKun), "{hat(%s)[%s]^{un}}(r)",
        "numerator of uncorrected estimate of %s", "un"
      )
      denK <- bind.fv(
        denK, data.frame(un = denKun), "{hat(%s)[%s]^{un}}(r)",
        "denominator of uncorrected estimate of %s",
        "un"
      )
    }
  }
  formula(K) <- . ~ r
  unitname(K) <- unitname(X)
  if (ratio) {
    formula(numK) <- formula(denK) <- . ~ r
    unitname(numK) <- unitname(denK) <- unitname(K)
    K <- rat(K, numK, denK, check = FALSE)
  }
  return(K)
}

#### studpermu.test ####
#' Studentised Permutation Test
#'
#' @description This is an S3 generic that extends the use of
#'   \code{\link[spatstat.core]{studpermu.test}} beyond "ppp" objects.
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat.core]{studpermu.test}},
#' \code{\link{studpermu.test.pp3}}
#'
#' @export
studpermu.test <- function(X, ...) UseMethod("studpermu.test")

### studpermu.test.list ###
#' Studentised Permutation Test
#'
#' @seealso \code{\link[spatstat.core]{studpermu.test}}
#' @export
studpermu.test.list <- spatstat.core::studpermu.test


### studpermu.test.hyperframe ###
#' Studentised Permutation Test
#'
#' @seealso \code{\link[spatstat.core]{studpermu.test}}
#' @export
studpermu.test.hyperframe <- function(X, ...) {
  h.class <- unclass(X)$vclass
  if (any(h.class == "ppp")) {
    studpermu.test.ppp(X, ... = ...)
  } else if (any(h.class == "pp3")) {
    studpermu.test.pp3(X, ... = ...)
  } else {
    stop("Unknown type for studpermu.test()")
  }
}

### studpermu.test.ppp ###
#' Studentised Permutation Test
#'
#' @seealso \code{\link[spatstat.core]{studpermu.test}}
#' @export
studpermu.test.ppp <- spatstat.core::studpermu.test

### studpermu.test.pp3 ###
#' Studentised Permutation Test
#'
#' Perform a studentised permutation test for a difference between groups of
#' point patterns
#'
#' @param X A hyperframe containing at least the point patterns and groups
#' @param formula Formula describing the grouping. The left side of the formula
#'   identifies which column of `X` contains the point patterns. The right
#'   side identifies the grouping factor
#' @param summaryfunction Summary function applicable to pp3. Defaults to
#'   \code{link[spatstat.core]{K3est}}
#' @param ... Additional arguments passed to `summaryfunction`
#' @param rinterval Numeric of length 2. Experimental
#' @param nperm Number of random permutations for the test; defaults to 999
#' @param use.Tbar Logical value indicating choice of test statistic. If TRUE,
#'   use the alternative test statistic, which is appropriate for summary
#'   functions with roughly constant variance, such as \eqn{K(r)/r} or
#'   \eqn{L(r)}. Defaults to FALSE
#' @param minpoints Minimum permissible number of points in a point pattern for
#'   inclusion in the test calculation
#'
#' @return An object of class "`studpermutest`".
#'
#' @family spatstat extensions
#'
#' @details
#'
#' This function performs the studentized permutation test of Hahn (2012) for a
#' difference between groups of point patterns.
#'
#' A group needs to contain at least two point patterns with at least
#' `minpoints` points in each pattern.
#'
#' The function returns an object of class "`htest`" and "`studpermutest`" that
#' can be printed and plotted. The printout shows the test result and *p*-value.
#' The plot shows the summary functions for the groups (and the group means if
#' requested).
#'
#' @references Hahn, U. (2012) A studentized permutation test for the comparison
#'   of spatial point patterns.
#'   *Journal of the American Statistical Association* **107** (498), 754-764.
#' @seealso \code{\link[spatstat.core]{studpermu.test}},
#'   \code{link[spatstat.core]{plot.studpermutest}}
#'
#' @export
# Add ability to supply a summary function directly...
studpermu.test.pp3 <- function(X, formula,
                               summaryfunction = K3est, ...,
                               nperm = 999, use.Tbar = FALSE, rinterval = NULL,
                               minpoints = 20, rmax = NULL, nrval = 128) {
  if (!is.hyperframe(X)) {
    stop(
      paste(
        "X needs to be a hyperframe",
        "if arguments for summary function are to be retrieved"
      ),
      call. = FALSE
    )
  }
  stopifnot(is.function(summaryfunction))
  if (is.hyperframe(X)) {
    if (dim(X)[2] < 2) {
      stop(paste(
        "Hyperframe X needs to contain at least 2 columns,",
        "one for patterns, one indicating groups"
      ), call. = FALSE)
    }
    data <- X
    Xclass <- unclass(X)$vclass
    factorcandidate <- Xclass %in% c(
      "integer", "numeric",
      "character", "factor"
    )
    ppcandidate <- Xclass == "pp3"
    Xnames <- names(X)
    names(factorcandidate) <- names(ppcandidate) <- names(Xclass) <- Xnames
    if (all(!factorcandidate) || all(!ppcandidate)) {
      stop(
        paste(
          "Hyperframe X needs to contain at least a column",
          "with point patterns, and one indicating groups"
        ),
        call. = FALSE
      )
    }
    if (!missing(formula)) {
      if (!inherits(formula, "formula")) {
        stop(paste("Argument", dQuote("formula"), "should be a formula"))
      }
      if (length(formula) < 3) {
        stop(paste("Argument", sQuote("formula"), "must have a left hand side"))
      }
      rhs <- spatstat.utils::rhs.of.formula(formula)
      ppname <- formula[[2]]
      if (!is.name(ppname)) {
        stop("Left hand side of formula should be a single name")
      }
      ppname <- paste(ppname)
      if (!ppcandidate[ppname]) {
        stop(
          paste(
            "Left hand side of formula",
            "should be the name of a column of point patterns"
          ),
          call. = FALSE
        )
      }
      groupvars <- all.vars(as.expression(rhs))
      if (!all(groupvars %in% Xnames) || any(!factorcandidate[groupvars])) {
        stop(paste(
          "Not all variables on right hand side of formula",
          "can be interpreted as factors"
        ), call. = FALSE)
      }
      group <- interaction(lapply(
        as.data.frame(data[, groupvars, drop = FALSE]), factor
      ))
      newnames <- Xnames
      newnames[Xnames == ppname] <- "pp"
      names(data) <- newnames
      data$group <- group
    } else {
      thepp <- which.max(ppcandidate)
      thegroup <- which.max(factorcandidate)
      formula <- as.formula(paste(Xnames[thepp], "~", Xnames[thegroup]))
      newnames <- Xnames
      newnames[thepp] <- "pp"
      newnames[thegroup] <- "group"
      names(data) <- newnames
      data$group <- as.factor(data$group)
    }
  } else {
    if (!is.list(X)) {
      stop("X should be a hyperframe or a list of lists of point patterns")
    }
    if (!is.list(X[[1]]) || !is.pp3(X[[1]][[1]])) {
      stop("X is a list, but not a list of lists of point patterns")
    }
    nams <- names(X)
    if (is.null(nams)) {
      nams <- paste("group", seq_along(X))
    }
    pp <- list()
    group <- NULL
    for (i in seq_along(X)) {
      pp <- c(pp, X[[i]])
      group <- c(group, rep(nams[i], length(X[[i]])))
    }
    group <- as.factor(group)
    data <- hyperframe(pp = pp, group = group)
    ppname <- "pp"
  }
  framename <- deparse(substitute(X))
  fooname <- deparse(substitute(summaryfunction))
  OK <- sapply(data$pp, npoints) >= minpoints
  if ((nbad <- sum(!OK)) > 0) {
    warning(paste(
      nbad, "patterns have been discarded",
      "because they contained fewer than",
      minpoints, "points"
    ), call. = FALSE)
  }
  data <- data[OK, , drop = FALSE]
  pp <- data$pp
  groupi <- as.integer(data$group)
  ngroups <- max(groupi)
  if (ngroups < 2) {
    stop(
      paste(
        "Sorry, after discarding patterns with fewer than",
        minpoints, "points,", if (ngroups < 1) {
          "nothing"
        } else {
          "only one group"
        }, "is left over.",
        "\n- nothing to compare, take a break!"
      ),
      call. = FALSE
    )
  }
  lev <- 1:ngroups
  m <- as.vector(table(groupi))
  if (any(m < 3)) {
    stop(
      paste(
        "Data groups need to contain at least two patterns;",
        "\nafter discarding those with fewer than", minpoints,
        "points, the remaining group sizes are",
        spatstat.utils::commasep(m)
      ),
      call. = FALSE
    )
  }
  npossible <- factorial(sum(m)) / prod(factorial(m)) / prod(factorial(table(m)))
  if (is.nan(npossible)) {
    npossible <- Inf
  }
  if (npossible < max(100, nperm)) {
    warning("Don't expect exact results - group sizes are too small")
  }
  if (is.null(rmax)) {
    rmax <- sapply(pp, function(P) {
      spatstat.geom::diameter(P$domain) / 2
    })
    rmax <- min(rmax)
  }
  if (is.null(rinterval)) {
    rinterval <- c(0, rmax)
  }
  ranger <- diff(range(rinterval))
  rr <- seq(0, rmax, length.out = nrval)
  taker <- rr >= rinterval[1] & rr <= rinterval[2]
  needcorx <- "correction" %in% names(formals(summaryfunction))
  gavecorx <- "correction" %in% names(list(...))
  corx <- if (needcorx && !gavecorx) {
    "translation"
  } else {
    NULL
  }
  fvlist <-
    if (is.null(corx)) {
      with(data, summaryfunction(pp, rmax = rmax, ...))
    } else {
      with(data, summaryfunction(pp, rmax = rmax, ..., correction = corx))
    }
  # can skip most of the previous if summary functions are supplied...
  # need to extract rmax, nrval from functions
  fvtemplate <- fvlist[[1]]
  valu <- attr(fvtemplate, "valu")
  argu <- attr(fvtemplate, "argu")
  foar <- sapply(lapply(fvlist, "[[", valu), "[", taker)
  combs <- combn(lev, 2)
  predigested <- list(
    lev = lev, foar = foar, m = m, combs = combs,
    rrr = rr[taker], ranger = ranger
  )
  if (use.Tbar) {
    Tobs <- T.barstat(groupi, predigested)
    Tsim <- replicate(nperm, T.barstat(sample(groupi), predigested))
  } else {
    Tobs <- T.stat(groupi, predigested)
    Tsim <- replicate(nperm, T.stat(sample(groupi), predigested))
  }
  names(Tobs) <- if (use.Tbar) {
    "Tbar"
  } else {
    "T"
  }
  pval <- (1 + sum(Tobs < Tsim)) / (1 + nperm)
  method <- c(
    "Studentized permutation test for grouped point patterns",
    ifelse(is.hyperframe(X),
      spatstat.utils::pasteFormula(formula), NULL
    ),
    spatstat.utils::choptext(
      ngroups, "groups:",
      paste(levels(data$group),
        collapse = ", "
      )
    ),
    spatstat.utils::choptext(
      "summary function:",
      paste0(fooname, ","), "evaluated on r in",
      spatstat.utils::prange(rinterval)
    ),
    spatstat.utils::choptext(
      "test statistic:",
      ifelse(use.Tbar, "Tbar,", "T,"),
      nperm, "random permutations"
    )
  )
  fooshort <- switch(fooname,
    pcf = "pair correlation ",
    Kinhom = "inhomogeneous K-",
    Linhom = "inhomogeneous L-",
    Kscaled = "locally scaled K-",
    Lscaled = "locally scaled L-",
    paste(substr(fooname, 1, 1), "-", sep = "")
  )
  alternative <- c(paste("not the same ", fooshort, "function",
    sep = ""
  ))
  testerg <- list(
    statistic = Tobs, p.value = pval, alternative = alternative,
    method = method, data.name = framename
  )
  class(testerg) <- c("studpermutest", "htest")
  fvs <- lapply(fvlist, "[.fv", j = c(argu, valu))
  fvs <- lapply(fvs, "attr<-", which = "alim", value = rinterval)
  testerg$curves <- spatstat.geom::hyperframe(fvs = fvs, groups = data$group)
  fvtheo <- fvlist[[1]]
  spatstat.core::fvnames(fvtheo, ".y") <- "theo"
  attr(fvtheo, "alim") <- rinterval
  testerg$curvtheo <- fvtheo[, c(argu, "theo")]
  grmn <- lapply(lev, splitmean, ind = groupi, f = foar)
  testerg$groupmeans <- lapply(grmn, makefv,
    xvals = rr[taker],
    template = fvtheo
  )
  return(testerg)
}

# These helpers are undocumented, but located at
# https://rdrr.io/github/spatstat/spatstat.core/src/R/studpermutest.R

splitmean <- function(l, ind, f) {
  apply(f[, ind == l], 1, mean)
}
splitvarn <- function(l, ind, f, m) {
  apply(f[, ind == l], 1, var) / m[l]
}
studentstat <- function(i, grmean, grvar) {
  (grmean[, i[1]] - grmean[, i[2]])^2 / (grvar[i[1], ] + grvar[i[2], ])
}

T.stat <- function(ind = groupi, predigested) {
  # predigested should be a list with entries lev, foar, m, combs, rrr
  with(predigested, {
    grmean <- sapply(lev, splitmean, ind = ind, f = foar)
    grvar <- t(sapply(lev, splitvarn, ind = ind, f = foar, m = m))
    y <- apply(combs, 2, studentstat, grmean = grmean, grvar = grvar)
    sum(apply(y, 2, trapint, x = rrr))
  })
}

intstudent <- function(i, rrr, grmean, meangrvar) {
  trapint(rrr, (grmean[, i[1]] - grmean[, i[2]])^2 /
    (meangrvar[i[1]] + meangrvar[i[2]]))
}

T.barstat <- function(ind = groupi, predigested) {
  # predigested should be a list
  # with entries lev, foar, m, combs, rrr, ranger
  with(predigested, {
    grmean <- sapply(lev, splitmean, ind = ind, f = foar)
    grvar <- t(sapply(lev, splitvarn, ind = ind, f = foar, m = m))
    meangrvar <- apply(grvar, 1, trapint, x = rrr) / ranger
    sum(apply(combs, 2, intstudent,
      rrr = rrr, grmean = grmean, meangrvar = meangrvar
    ))
    # trapint(rr[taker], grvar[i[1],] + grvar[i[2], ]))))
  })
}

makefv <- function(yvals, xvals, template) {
  fdf <- data.frame(r = xvals, y = yvals)
  argu <- fvnames(template, ".x")
  valu <- fvnames(template, ".y")
  names(fdf) <- c(argu, valu)
  fv(fdf,
    argu = argu, ylab = attr(template, "ylab"), valu = valu,
    fmla = attr(template, "fmla"), alim = attr(template, "alim")
  )
}

# Trapezoidal rule approximation to integral
trapint <- function(x, y) {
  nonan <- !is.na(y)
  nn <- sum(nonan)
  if (nn < 2L) {
    return(0)
  }
  Y <- y[nonan]
  X <- x[nonan]
  0.5 * sum((Y[-1] + Y[-nn]) * diff(X))
}

# call foo(x, further arguments) repeatedly
# further arguments are taken from hyperframe H and ...
multicall <- function(foo, x, H, ...) {
  stopifnot(is.hyperframe(H))
  if (is.hyperframe(x)) {
    x <- as.list(x)[[1]]
  } else if (!is.list(x)) {
    stop("in multicall: x should be a hyperframe or list", call. = FALSE)
  }

  # check if same length
  nrows <- dim(H)[1]
  if (length(x) != nrows) {
    stop(
      paste(
        "in multicall: x and H need to have",
        "the same number of rows or list elements"
      ),
      call. = FALSE
    )
  }
  dotargs <- list(...)
  hnames <- names(H)
  argnames <- names(formals(foo)) # always assume first argument is given

  ppname <- argnames[1]
  argnames <- argnames[-1]
  dotmatch <- pmatch(names(dotargs), argnames)
  dotmatched <- dotmatch[!is.na(dotmatch)]
  dotuseargs <- dotargs[!is.na(dotmatch)]
  restargs <- if (length(dotmatched) > 0) argnames[-dotmatched] else argnames
  hmatch <- pmatch(hnames, restargs)
  huse <- !is.na(hmatch)
  lapply(seq_len(nrows), function(i) {
    do.call(foo, c(
      list(x[[i]]),
      as.list(H[i, huse, drop = TRUE, strip = FALSE]),
      dotargs
    ))
  })
}

#### Tstat.pp3 ####
#' Extends Tstat to pp3
#'
#' Tstat.pp3 extends the third-order summary statistic
#' \code{\link[spatstat.core]{Tstat}} to pp3
#'
#' @param X The observed point pattern, from which an estimate of \eqn{T(r)}
#' will be computed. A \code{\link[spatstat.geom]{pp3}} object.
#' @param rmax Optional. Maximum value of argument *r* for which
#'   \eqn{T(r)} will be estimated.
#' @param nrval Optional. Number of values of *r* for which
#'   \eqn{T(r)} will be estimated. A large value of `nrval` is
#'   required to avoid discretisation effects.
#' @param correction One of `"none"` or `"isotropic"`. `"translation"`
#'   correction is planned but not yet implemented.
#' @param ratio Logical. If `TRUE`, the numerator and denominator of each
#'   edge-corrected estimate will also be saved, for use in analysing replicated
#'   point patterns.
#' @param verbose Logical. If `TRUE`, an estimate of the computation time is
#'   printed.
#'
#' @return An object of class "`fv`", see \code{\link[spatstat]{fv.object}},
#'   which can be plotted directly using \code{\link[spatstat]{plot.fv}}.
#'
#' @family spatstat extensions
#'
#' @details
#'
#' This command calculates the third-order summary statistic \eqn{T(r)} for a
#' spatial point patterns, defined by Schladitz and Baddeley (2000).
#'
#' The definition of \eqn{T(r)} is similar to the definition of Ripley's *K*
#' function \eqn{K(r)}, except that \eqn{K(r)} counts pairs of points while
#' \eqn{T(r)} counts triples of points. Essentially \eqn{T(r)} is a rescaled
#' cumulative distribution function of the diameters of triangles in the point
#' pattern. The diameter of a triangle is the length of its longest side.
#'
#' @seealso \code{\link[spatstat]{Tstat}}
#' @references Schladitz, K. & Baddeley, A.
#' "A third order point process characteristic",
#' *Scandinavian Journal of Statistics*, **27**, 657-671 (2000).
#' @export
Tstat.pp3 <- function(X, rmax = NULL, nrval = 128,
                      correction = "border",
                      ratio = FALSE, verbose = TRUE) {
  spatstat.geom::verifyclass(X, "pp3")
  npts <- spatstat.geom::npoints(X)
  W <- spatstat.geom::domain(X)
  areaW <- spatstat.geom::volume(W)
  lambda <- npts / areaW
  lambda2 <- (npts * (npts - 1)) / (areaW^2)
  lambda3 <- (npts * (npts - 1) * (npts - 2)) / (areaW^3)
  if (is.null(rmax)) {
    rmax <- spatstat.geom::diameter(W) / 2
  }
  r <- seq(0, rmax, length.out = nrval)
  breaks <- spatstat.geom::breakpts.from.r(r)
  correction.given <- !missing(correction) && !is.null(correction)
  if (!correction.given) {
    correction <- c("border", "bord.modif")
  } # , "translate") not implemented yet
  correction <- spatstat.geom::pickoption(
    "correction", correction,
    c(
      none = "none",
      border = "border", bord.modif = "bord.modif", trans = "translate",
      translate = "translate", translation = "translate", best = "best"
    ),
    multi = TRUE
  )
  alim <- c(0, rmax)
  TT <- data.frame(r = r, theo = 5 / 12 * pi^2 * r^6)
  desc <- c("distance argument r", "theoretical Poisson %s")
  TT <- spatstat.core::fv(TT, "r",
    quote(T(r)), "theo", , alim, c("r", "%s[pois](r)"),
    desc,
    fname = "T"
  ) # blank is in ppp version...
  if (ratio) {
    denom <- lambda2 * areaW
    numT <- spatstat.core::eval.fv(denom * TT)
    denT <- spatstat.core::eval.fv(denom + TT * 0)
    attributes(numT) <- attributes(denT) <- attributes(TT)
    attr(numT, "desc")[2] <- "numerator for theoretical Poisson %s"
    attr(denT, "desc")[2] <- "denominator for theoretical Poisson %s"
  }
  close <- spatstat.geom::closepairs(X, rmax,
    what = "ijd", twice = FALSE,
    neat = FALSE
  )
  I <- close$i
  J <- close$j
  DIJ <- close$d
  nI <- length(I)
  continue <- TRUE
  if (verbose) {
    nTmax <- nI * (nI - 1) / 2
    esttime <- exp(1.25 * log(nTmax) - 21.5)
    message(paste(
      "Searching", nTmax, "potential triangles;",
      "estimated time", spatstat.geom::codetime(esttime)
    ))
    if (esttime > 60) {
      continue <- askYesNo("Estimated time greater than one minute. Continue?")
    }
  }
  if (!continue) {
    return(NULL)
  }
  tri <- spatstat.geom::trianglediameters(I, J, DIJ, nvert = npts)
  stopifnot(identical(colnames(tri), c("i", "j", "k", "diam")))
  II <- with(tri, c(i, j, k))
  DD <- with(tri, rep.int(diam, 3))
  if (any(correction == "none")) {
    wh <- whist(DD, breaks$val)
    numTun <- cumsum(wh)
    denTun <- lambda3 * areaW
    Tun <- numTun / denTun
    TT <- spatstat.core::bind.fv(
      TT, data.frame(un = Tun), "hat(%s)[un](r)",
      "uncorrected estimate of %s", "un"
    )
    if (ratio) {
      numT <- spatstat.core::bind.fv(
        numT, data.frame(un = numTun), "hat(%s)[un](r)",
        "numerator of uncorrected estimate of %s", "un"
      )
      denT <- spatstat.core::bind.fv(
        denT, data.frame(un = denTun), "hat(%s)[un](r)",
        "denominator of uncorrected estimate of %s",
        "un"
      )
    }
  }
  if (any(correction == "border" | correction == "bord.modif")) {
    b <- bdist.points.pp3(X)
    bI <- b[II]
    RS <- spatstat.core::Kount(DD, bI, b, breaks)
    if (any(correction == "bord.modif")) {
      denom.area <- eroded.volumes(W, r)
      numTbm <- RS$numerator
      denTbm <- lambda3 * denom.area
      Tbm <- numTbm / denTbm
      TT <- spatstat.core::bind.fv(
        TT, data.frame(bord.modif = Tbm),
        "hat(%s)[bordm](r)",
        "modified border-corrected estimate of %s",
        "bord.modif"
      )
      if (ratio) {
        numT <- spatstat.core::bind.fv(
          numT, data.frame(bord.modif = numTbm),
          "hat(%s)[bordm](r)",
          "numerator of modified border-corrected estimate of %s",
          "bord.modif"
        )
        denT <- spatstat.core::bind.fv(
          denT, data.frame(bord.modif = denTbm),
          "hat(%s)[bordm](r)",
          "denominator of modified border-corrected estimate of %s",
          "bord.modif"
        )
      }
    }
    if (any(correction == "border")) {
      numTb <- RS$numerator
      denTb <- lambda2 * RS$denom.count
      Tb <- numTb / denTb
      TT <- spatstat.core::bind.fv(
        TT, data.frame(border = Tb), "hat(%s)[bord](r)",
        "border-corrected estimate of %s", "border"
      )
      if (ratio) {
        numT <- spatstat.core::bind.fv(
          numT, data.frame(border = numTb),
          "hat(%s)[bord](r)",
          "numerator of border-corrected estimate of %s",
          "border"
        )
        denT <- spatstat.core::bind.fv(
          denT, data.frame(border = denTb),
          "hat(%s)[bord](r)",
          "denominator of border-corrected estimate of %s",
          "border"
        )
      }
    }
  }
  if (any(correction == "translate")) {
    stop("Translation correction not yet implemented.")
    edgewt <- edgetri.Trans(X, tri[, 1:3])
    wh <- whist(tri$diam, breaks$val, edgewt)
    numTtrans <- 3 * cumsum(wh)
    denTtrans <- lambda3 * areaW
    Ttrans <- numTtrans / denTtrans
    h <- spatstat.geom::diameter(W) / 2
    Ttrans[r >= h] <- NA
    TT <- spatstat.core::bind.fv(
      TT, data.frame(trans = Ttrans), "hat(%s)[trans](r)",
      "translation-corrected estimate of %s", "trans"
    )
    if (ratio) {
      numT <- spatstat.core::bind.fv(
        numT, data.frame(trans = numTtrans),
        "hat(%s)[trans](r)",
        "numerator of translation-corrected estimate of %s",
        "trans"
      )
      denT <- spatstat.core::bind.fv(
        denT, data.frame(trans = denTtrans),
        "hat(%s)[trans](r)",
        "denominator of translation-corrected estimate of %s",
        "trans"
      )
    }
  }
  formula(TT) <- . ~ r
  spatstat.geom::unitname(TT) <- spatstat.geom::unitname(X)
  if (ratio) {
    formula(numT) <- formula(denT) <- . ~ r
    spatstat.geom::unitname(denT) <- spatstat.geom::unitname(TT)
    spatstat.geom::unitname(numT) <- spatstat.geom::unitname(TT)
    TT <- spatstat.core::rat(TT, numT, denT, check = FALSE)
  }
  return(TT)
}

#### bdist.points ####
#' Distance to Boundary of Domain
#'
#' @description This is an S3 generic that extends the use of
#'   \code{\link[spatstat.geom]{bdist.points}} beyond "ppp" objects.
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat.geom]{bdist.points}}, \code{\link{bdist.points.pp3}}
#'
#' @export
bdist.points <- function(X, ...) UseMethod("bdist.points")

### bdist.points.ppp ###
#' Distance to Boundary of Window
#'
#' @seealso \code{\link[spatstat.geom]{bdist.points}}
#' @export
bdist.points.ppp <- spatstat.geom::bdist.points

### bdist.points.pp3 ###
#' Distance to Boundary of Domain
#'
#' Finds the smallest distance to a boundary for each point in a point pattern.
#'
#' @param X The point pattern for analysis. A \code{\link[spatstat.geom]{pp3}}
#'   object.
#'
#' @return An object containing the shortest distance to the boundary for each
#'   point in the pattern X.
#' @seealso \code{\link[spatstat]{bdist.points}}
#'
#' @export
bdist.points.pp3 <- function(X) {
  spatstat.geom::verifyclass(X, "pp3")

  x <- X$data$x
  y <- X$data$y
  z <- X$data$z
  d <- X$domain

  xmin <- min(d$xrange)
  xmax <- max(d$xrange)
  ymin <- min(d$yrange)
  ymax <- max(d$yrange)
  zmin <- min(d$zrange)
  zmax <- max(d$zrange)
  result <- pmin.int(
    x - xmin, xmax - x,
    y - ymin, ymax - y,
    z - zmin, zmax - z
  )

  return(result)
}

#### bdist.points3.multi ####
#' Cardinal Direction Distances to Boundary of Domain
#'
#' Returns the shortest distances to boundaries in the x, y, and z directions
#' separately.
#'
#' @param X The point pattern for analysis. A \code{\link[spatstat.geom]{pp3}}
#'   object.
#' @return A `data.frame` containing the shortest distance to the closest three
#'   boundaries for each point in the pattern `X`.
#' @seealso \code{\link{bdist.points.pp3}}
#' @export
bdist.points3.multi <- function(X) {
  spatstat.geom::verifyclass(X, "pp3")

  x <- X$data$x
  y <- X$data$y
  z <- X$data$z
  d <- X$domain

  xmin <- min(d$xrange)
  xmax <- max(d$xrange)
  ymin <- min(d$yrange)
  ymax <- max(d$yrange)
  zmin <- min(d$zrange)
  zmax <- max(d$zrange)

  result <- data.frame(
    x = pmin.int(x - xmin, xmax - x),
    y = pmin.int(y - ymin, ymax - y),
    z = pmin.int(z - zmin, zmax - z)
  )

  return(result)
}
