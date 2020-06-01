#
# This file contains methods for simulating APT data.
#

#### rpoint3 ####
#' Generate N Random Points
#'
#' \code{rpoint3} extends \code{\link[spatstat]{rpoint}} to
#' \code{\link[spatstat]{pp3}}.
#' @seealso \code{\link[spatstat]{rpoint}}
#' @export
rpoint3 <- function (n, f, fmax = 1,  win = box3(), ...,
          giveup = 1000, verbose = FALSE, nsim = 1, drop = TRUE)
{
  if (missing(f) || (is.numeric(f) && length(f) == 1))
    return(runifpoint3(n, domain = win, nsim = nsim, drop = drop))
  if (!is.function(f) && !is.im(f))
    stop(paste(sQuote("f"), "must be a function"))
  verifyclass(win, "box3")
  if (n == 0) {
    emp <- pp3(numeric(0), numeric(0), numeric(0), window = win)
    if (nsim == 1 && drop)
      return(emp)
    result <- rep(list(emp), nsim)
    names(result) <- paste("Simulation", 1:nsim)
    return(as.anylist(result))
  }
  result <- vector(mode = "list", length = nsim)
  for (isim in 1:nsim) {
    X <- pp3(numeric(0), numeric(0), numeric(0), window = win)
    pbar <- 1
    nremaining <- n
    totngen <- 0
    ntries <- 0
    repeat {
      ntries <- ntries + 1
      ngen <- nremaining/pbar + 10
      totngen <- totngen + ngen
      prop <- runifpoint3(ngen, domain = win)
      if (npoints(prop) > 0) {
        fvalues <- f(prop$data$x, prop$data$y, prop$data$z, ...)
        paccept <- fvalues/fmax
        u <- runif(npoints(prop))
        Y <- prop[u < paccept]
        if (npoints(Y) > 0) {
          X <- superimpose(X, Y, W = win)
          nX <- npoints(X)
          pbar <- nX/totngen
          nremaining <- n - nX
          if (nremaining <= 0) {
            if (verbose)
              splat("acceptance rate = ", round(100 * pbar, 2), "%")
            result[[isim]] <- if (nX == n)
              X
            else X[1:n]
            break
          }
        }
      }
      if (ntries > giveup)
        stop(paste("Gave up after", giveup * n, "trials with",
                   npoints(X), "points accepted"))
    }
  }
  if (nsim == 1 && drop)
    return(result[[1]])
  names(result) <- paste("Simulation", 1:nsim)
  return(as.anylist(result))
}
#### rPoissonCluster3 ####
#' Simulate 3D Poisson Cluster Process
#'
#' \code{rPoissonCluster3} extends \code{\link[spatstat]{rPoissonCluster}} to
#' \code{\link[spatstat]{pp3}}.
#' @seealso \code{\link[spatstat]{rpoint}}
#' @export
rPoissonCluster3 <- function(kappa, expand, rcluster, win = box3(), ...,
                             nsim = 1, drop = T)
{
  if (missing(expand) && !is.null(rmax <- list(...)$rmax)) {
    expand <- rmax
    f <- rcluster
    rcluster <- function(..., rmax) f(...)
  }
  verifyclass(win, "box3")
  dilated <- box3(win$xrange + c(-expand, expand),
                  win$yrange + c(-expand, expand),
                  win$yrange + c(-expand, expand))
  parentlist <- rpoispp3(kappa, domain = dilated,
                        nsim = nsim)
  if (nsim == 1)
    parentlist <- list(parentlist)
  resultlist <- vector(mode = "list", length = nsim)
  for (isim in 1:nsim) {
    parents <- parentlist[[isim]]
    result <- NULL
    res.full <- NULL
    np <- npoints(parents)
    if (np > 0) {
      xparent <- parents$data$x
      yparent <- parents$data$y
      zparent <- parents$data$z
      for (i in seq_len(np)) {
        cluster <- rcluster(...)
        if (npoints(cluster) > 0) {
          cluster <- shift(cluster, vec = c(xparent[i], yparent[i], zparent[i]))
          cluster <- pp3(cluster$data$x, cluster$data$y, cluster$data$z,
                         xrange = win$xrange,
                         yrange = win$yrange,
                         zrange = win$zrange)
          clus.trunc <- subset(cluster, subset =
                                 (x > win$xrange[1] & x < win$xrange[2]) &
                                 (y > win$yrange[1] & y < win$yrange[2]) &
                                 (z > win$zrange[1] & z < win$zrange[2])
                               )
          if (is.null(result)) {
            result <- clus.trunc
            res.full <- cluster
            parentid <- rep.int(1, npoints(clus.trunc))
          }
          else {
            result <- superimpose(result, clus.trunc, W = win)
            res.full <- superimpose(res.full, cluster, W = win)
            parentid <- c(parentid, rep.int(i, npoints(clus.trunc)))
          }
        }
      }
    }
    else {
      result <- pp3(numeric(0), numeric(0), numeric(0), window = win)
      parentid <- integer(0)
    }
    attr(result, "parents") <- parents
    attr(result, "parentid") <- parentid
    attr(result, "expand") <- expand
    attr(result, "full") <- res.full
    resultlist[[isim]] <- result
  }
  if (nsim == 1 && drop)
    return(resultlist[[1]])
  names(resultlist) <- paste("Simulation", 1:nsim)
  return(as.anylist(resultlist))
}
#### latticeVectors ####
#' Transform a matrix of integer indices to spatial coordinates for a specified
#' lattice type
#' @export
latticeVectors <- function(indices, a = 1, lattice = "sc") {
  indices <- as.matrix(indices)
  if (lattice == "sc") {
    lat.dat <- apply(indices, 1, function(vec) {
      x <- a * vec[1]
      y <- a * vec[2]
      z <- a * vec[3]
      return(c(x, y, z))
    })
  } else if (lattice == "bcc") {
    lat.dat <- apply(indices, 1, function(vec) {
      x <- a/2 * (-vec[1] + vec[2] + vec[3])
      y <- a/2 * (vec[1] - vec[2] + vec[3])
      z <- a/2 * (vec[1] + vec[2] - vec[3])
      return(c(x, y, z))
    })
  } else if (lattice == "fcc") {
    lat.dat <- apply(indices, 1, function(vec) {
      x <- a/2 * (vec[2] + vec[3])
      y <- a/2 * (vec[1] + vec[3])
      z <- a/2 * (vec[1] + vec[2])
      return(c(x, y, z))
    })
  } else {
    warning("Unrecognized lattice")
  }
  lat.dat <- t(lat.dat)
  return(lat.dat)
}
#### lattice ####
#' Generate Spatial Lattice
#'
#' \code{lattice} creates a spatial region filled with lattice points of the specified type
#' @export
lattice <- function(domain = box3(), a = 1, lattice = "sc") {
  lat.xrange <- round(domain$xrange / a)
  lat.yrange <- round(domain$yrange / a)
  lat.zrange <- round(domain$zrange / a)
  lat.expand <- max(sapply(list(lat.xrange,lat.yrange,lat.zrange),
                           diff, simplify = T))
  lat.x <- seq(lat.xrange[1] - lat.expand, lat.xrange[2] + lat.expand)
  lat.y <- seq(lat.yrange[1] - lat.expand, lat.yrange[2] + lat.expand)
  lat.z <- seq(lat.zrange[1] - lat.expand, lat.zrange[2] + lat.expand)
  lat.index <- expand.grid(lat.x, lat.y, lat.z)
  names(lat.index) <- NULL
  lat.dat <- latticeVectors(lat.index, a = a, lattice = lattice)
  lat.pp3 <- pp3(x = lat.dat[,1], y = lat.dat[,2], z = lat.dat[,3], domain)
  ok <- inside.pp3(lat.pp3)
  lat.pp3 <- lat.pp3[ok]
  return(lat.pp3)
}

#### nmers ####
#' Creates an n-mer point pattern by selecting the n+1 nearest neighbors of a
#' sampled point pattern from its parent distribution.
nmers <- function(sam, parent, n = 2) {
  nmer.ind <- nncross(sam, parent, what = "which", k = 1:n)
  nmer.dat <- parent[unlist(nmer.ind)]
  return(nmer.dat)
}
