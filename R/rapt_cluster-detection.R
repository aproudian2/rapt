#### msa ####
#' Identify Clusters in a Marked Point Pattern Using MSA
#'
#' `msa` segments a marked \code{\link[spatstat.geom]{pp3}} into clusters and
#' background matrix using the maximum separation algorithm (MSA). The marks can
#' have more than two types, but MSA requires that each type is categorized as
#' either a cluster or non-cluster species.
#'
#' @param X A marked \code{\link[spatstat.geom]{pp3}} object on which MSA will be
#'   performed.
#' @param dmax The maximum distance two points can be separated by and still be
#'   considered part of the same cluster.
#' @param Nmin The minimum number of points needed to classify a grouping of
#'   points as a cluster.
#' @param denv Any points within this distance of a point residing in a cluster
#'   will be included in the cluster (this controls the addition of background
#'   points to a cluster).
#' @param der Any points within this distance of a background matrix point
#'   (after enveloping) will be removed from the cluster.
#' @param clust.mark Vector containing the names of the marks in `X` that
#'   should be considered as cluster type points. All points with marks not
#'   included in this vector will be considered background points.
#'
#' @return A list of:
#' * `radius` - A vector containing estimated radius of each cluster found
#' * `den` - A vector containing estimated intra-cluster concentration of
#'   cluster type points in each cluster found
#' * `bgnd.den` - The estimated background concentration of cluster type points
#'   within the background (*i.e.* the entire domain minus the clusters)
#' * `cluster` - A list of indices from the original pattern of cluster type
#'   points that reside in clusters
#' * `bgnd` - A list of indices from the original pattern of background type
#'   points that reside in clusters.
#'
#' @family cluster identification functions
#'
#' @references
#' Marquis, E.A. & Hyde, J.M.,
#' "Applications of atom-probe tomography to the characterisation of solute
#' behaviours,"
#' *Materials Science and Engineering: R: Reports*, **69** (4-5), 37-62 (2010):
#' <https://doi.org/10.1016/j.mser.2010.05.001>
#'
#' @export
# Add a new marked pp3 with marks corresponding to ID'ed background and
# individual cluster points
msa <- function(X, dmax, Nmin, denv, der, clust.mark = c("A")) {
  verifyclass(X, "pp3")
  if (!is.marked(X)) {
    gripe <- paste(dQuote(X), "must be a marked point pattern")
    stop(gripe)
  }
  X.A <- X[marks(X) %in% clust.mark]
  X.B <- X[!(marks(X) %in% clust.mark)]

  marks(X.A) <- which(marks(X) %in% clust.mark)
  marks(X.B) <- which(!(marks(X) %in% clust.mark))

  # find the nns within dmax of all type A points:
  cp <- closepairs(X.A, rmax = dmax, twice = TRUE, what = "indices")
  cp <- data.frame("i" = cp$i[order(cp$i)], "j" = cp$j[order(cp$i)])

  # change results to list (nns): for each i, create an entry in a list
  # that contains a vector of its nns
  nns <- list()
  nns.gb <- dplyr::group_by(cp, i)
  nns.labs <- attr(nns.gb, "groups")$i
  nns.inds <- attr(nns.gb, "groups")$.rows
  nns.inds.no <- (1:npoints(X.A))[-unlist(nns.labs)]

  for (k in 1:length(nns.labs)) {
    nns[[nns.labs[k]]] <- cp$j[nns.inds[[k]]]
  }
  nns[nns.inds.no] <- c(0)

  # Find individual clusers with more points than Nmax
  diveDeep <- function(is, inds) {
    inds.all <- c()
    for (j in is) {
      if (nns[[j]][1] == 0) {
        next
      }
      inds.all <- append(inds.all, nns[[j]])
    }
    inds.all <- unique(inds.all)
    inds.new <- inds.all[!inds.all %in% inds]

    # base case
    if (length(inds.new) == 0) {
      return(inds)
    }

    # recursive case
    inds <- append(inds, inds.new)
    return(diveDeep(inds.new, inds))
  }

  clusters <- list()
  # ind.list <- 1:npoints(X.A)
  ind.list <- nns.labs
  cnt <- 1

  while (length(ind.list) > 0) {
    to.do <- ind.list[1]
    clusters[[cnt]] <- diveDeep(to.do, to.do)

    ind.list <- ind.list[!ind.list %in% clusters[[cnt]]]

    cnt <- cnt + 1
  }

  cluster.sizes <- sapply(clusters, length)
  clusters.pass <- clusters[cluster.sizes >= Nmin]
  if (length(clusters.pass) == 0) {
    print("No clusters found")
    return(NA)
  }

  # Get background points in clusters
  X.clusters.A <- X.A[unlist(clusters.pass)]
  mks <- rep(1:length(clusters.pass), sapply(clusters.pass, length))
  marks(X.clusters.A) <- mks

  B.in.clusters <- list()

  cp.AB <- crosspairs(X.clusters.A, X.B, rmax = denv, what = "indices")
  cp.AB <- data.frame(
    "i" = cp.AB$i[order(cp.AB$i)],
    "j" = cp.AB$j[order(cp.AB$i)]
  )

  AB.gb <- dplyr::group_by(cp.AB, i)
  AB.labs <- attr(AB.gb, "groups")$i
  AB.inds <- attr(AB.gb, "groups")$.rows

  for (k in 1:length(clusters.pass)) {
    clust.inds <- which(marks(X.clusters.A) == k)
    clust.inds <- clust.inds[clust.inds %in% AB.labs]
    data.inds <- which((AB.labs %in% clust.inds) == TRUE)
    to.pull <- unlist(lapply(data.inds, function(l) {
      AB.inds[[l]]
    }))
    B.in.clusters[[k]] <- unique(cp.AB$j[to.pull])
  }

  # crosspairs to erode
  matrix.bgnd <- X.B[-unique(unlist(B.in.clusters))]
  matrix.A <- X.A[-unique(unlist(clusters.pass))]
  # combine into full matrix
  coo.mat <- rbind(coords(matrix.bgnd), coords(matrix.A))
  matrix.all <- pp3(coo.mat$x, coo.mat$y, coo.mat$z, domain(X))

  # make pp3s of the cluster stuff so far:
  marks(X.clusters.A) <- unlist(clusters.pass)

  X.clusters.B <- X.B[unlist(B.in.clusters)]
  marks(X.clusters.B) <- unlist(B.in.clusters)

  # Check these
  cp.AM <- crosspairs(X.clusters.A, matrix.all, rmax = der, what = "indices")
  cp.BM <- crosspairs(X.clusters.B, matrix.all, rmax = der, what = "indices")

  A.remove <- marks(X.clusters.A[unique(cp.AM$i)])
  B.remove <- marks(X.clusters.B[unique(cp.BM$i)])


  A.clusters.eroded <- lapply(clusters.pass, function(x) {
    return(x[!x %in% A.remove])
  })

  B.clusters.eroded <- lapply(B.in.clusters, function(x) {
    return(x[!x %in% B.remove])
  })

  nclusters <- length(A.clusters.eroded)

  # Intra-cluster density
  cluster.den <- sapply(1:nclusters, function(x) {
    length(A.clusters.eroded[[x]]) /
      (length(B.clusters.eroded[[x]]) + length(A.clusters.eroded[[x]]))
  })

  # background density
  bgnd.total <- npoints(X) - length(unlist(A.clusters.eroded)) -
    length(unlist(B.clusters.eroded))
  bgnd.A <- npoints(X.A) - length(unlist(A.clusters.eroded))
  bgnd.den <- bgnd.A / bgnd.total

  # Guinier Radius
  cluster.Rg <- sapply(A.clusters.eroded, function(x) {
    coo <- coords(X.A[x])
    com <- apply(coo, 2, mean)
    rs2 <- apply(t(t(coo) - com), 1, function(y) {
      sum(y^2)
    })
    Rg <- sqrt(sum(rs2) / nrow(coo))
    Dg <- 2 * sqrt(5 / 3) * Rg
    return(Dg / 2)
  })

  A.cluster.inds.orig <- lapply(A.clusters.eroded, function(x) {
    marks(X.A)[x]
  })
  B.cluster.inds.orig <- lapply(B.clusters.eroded, function(x) {
    marks(X.B)[x]
  })

  return(list(
    "radius" = cluster.Rg, "den" = cluster.den, "bgnd.den" = bgnd.den,
    "cluster" = A.cluster.inds.orig, "bgnd" = B.cluster.inds.orig
  ))
}

#### gema ####
#' Identify Clusters in a Point Pattern Using GEMA
#' @family cluster identification functions
gema <- function(X, ...) UseMethod("gema")

### gema.pp3 ###
#' Identify Clusters in a Point Pattern Using GEMA
#'
#' @param X The point pattern (object of class \code{\link[spatstat.geom]{ppp}} or
#'   \code{\link[spatstat.geom]{pp3}}) in which to identify clusters.
#' @param cluster The marks of `X` that are cluster-type points. If
#'   `cluster = NULL` (the default) or `X` is unmarked, all points are assumed
#'   to be cluster-type points.
#' @param kde.n Number of sampling grid lines for the 3D kernel density estimate
#' @param max.clusters The maximum number of clusters to define the cluster
#'   search initialization. See Details.
#' @param threshold The probability threshold at which to assign a point to a
#'   cluster or to the background. Defaults to 0.5, which is usually Bayes
#'   optimal.
#' @return A marked point pattern (of the same class as `X`) with marks
#'   corresponding to the identified background and each identified cluster
#'   based on `threshold`.
#'
#' @details
#' A `data.frame` is returned as an attribute (`attr("prob")`) that contains the
#' probabilities of assignment to the background and each cluster.
#'
#' @family cluster identification functions
#'
#' @references Zelenty, J. *et al.*,
#' "Detecting Clusters in Atom Probe Data with Gaussian Mixture Models",
#' *Microscopy and Microanalysis*, **23** (2), 269-278 (2017):
#' <https://doi.org/10.1017/S1431927617000320>
gema.pp3 <- function(X, cluster = NULL,
                     kde.n = 20, max.clusters = 30,
                     threshold = 0.5) {
  stopifnot(inherits(X, c("ppp", "pp3")))
  stopifnot(threshold >= 0 & threshold <= 1)
  if (!is.null(cluster) & is.marked(X)) {
    Y <- unmark(X[marks(X) %in% cluster])
  } else {
    Y <- unmark(X)
  }

  kd0 <- misc3d::kde3d(Y$data$x, Y$data$y, Y$data$z, n = kde.n)
  sn <- runifpoint3(npoints(Y), domain = domain(Y))
  ks <- misc3d::kde3d(sn$data$x, sn$data$y, sn$data$z, n = kde.n)
  ks$d <- kd0$d - ks$d

  pts <- clustInit3(ks)
  pt3 <- pp3(pts$x, pts$y, pts$z, domain(Y))
  sinit <- crosspairs(pt3, Y, rmax = 4, what = "indices")
  sinit <- unique(unlist(sinit))
  ninit <- setdiff(1:npoints(Y), sinit)
  init <- list(subset = sinit, noise = ninit)
  gmlist <- vector("list", max.clusters)

  # cat('1, ')
  gmlist[[1]] <- Mclust(as.data.frame(Y$data),
    G = 1, modelNames = "VII", initialization = init,
    verbose = FALSE
  )
  for (j in 2:max.clusters) {
    # cat(paste0(j, ', '))
    # Assumes continuous position and uniform background
    sn <- runifpoint3(npoints(Y) * tail(gmlist$parameters$pro, n = 1),
      domain = domain(Y)
    )
    sg <- lapply(seq_len(j - 1), function(k) {
      MASS::mvrnorm(
        ceiling(npoints(Y) * gmlist[[j - 1]]$parameters$pro[k]),
        gmlist[[j - 1]]$parameters$mean[, k],
        gmlist[[j - 1]]$parameters$variance$sigma[, , k]
      )
    })
    sg <- do.call(rbind, sg)
    sg <- pp3(sg[, 1], sg[, 2], sg[, 3], domain(Y))
    sg <- sg[inside.boxx(sg, w = domain(sg))]
    scl <- superimpose(sn, sg)
    ks <- misc3d::kde3d(scl$data$x, scl$data$y, scl$data$z, n = kde.n)
    ks$d <- kd0$d - ks$d
    ks$d[ks$d < 0] <- 0
    pts <- clustInit3(ks, pts)
    pt3 <- pp3(pts$x, pts$y, pts$z, domain(Y))
    sinit <- crosspairs(pt3, Y, rmax = 4, what = "indices")
    sinit <- unique(unlist(sinit))
    ninit <- setdiff(1:npoints(Y), sinit)
    init <- list(subset = sinit, noise = ninit)
    gmlist[[j]] <- Mclust(as.data.frame(Y$data),
      G = j, modelNames = "VII", initialization = init,
      verbose = FALSE
    )
  }
  cat("Done.", fill = TRUE)

  attr(gmlist, "seed.pts") <- pts
  return(gmlist)
}

### gema.ppp ###
#' @inheritParams gema.pp3
#' @family cluster identification functions
gema.ppp <- gema.pp3

### gema.gema ###
#' Reprocess Clusters Identified Using GEMA
#'
#' @param X An object of class "gema" (as produced by \code{\link{gema.pp3}})
#' @param threshold The probability threshold at which to assign a point to a
#'   cluster or to the background.
#' @return An
#'
#' @family cluster identification functions
gema.gema <- function(X, threshold = 0.5) {}

clustInit3 <- function(kd, pts = NULL) {
  ind <- which(kd$d == max(kd$d), arr.ind = TRUE)
  p <- data.frame(x = kd$x[ind[1]], y = kd$y[ind[2]], z = kd$z[ind[2]])
  if (is.null(pts)) {
    dat <- p
  } else {
    dat <- rbind(pts, p)
  }
  return(dat)
}
clustInit2 <- function(kd, pts = NULL) {
  ind <- which(kd$z == max(kd$z), arr.ind = TRUE)
  p <- data.frame(x = kd$x[ind[1]], y = kd$y[ind[2]])
  if (is.null(pts)) {
    dat <- p
  } else {
    dat <- rbind(pts, p)
  }
  return(dat)
}
