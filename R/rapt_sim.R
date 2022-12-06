#
# This file contains methods for simulating APT data.
#

#### rpoint3 ####
#' Generate N Random Points
#'
#' \code{rpoint3} extends \code{\link[spatstat.geom]{rpoint}} to
#' \code{\link[spatstat.geom]{pp3}}.
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat.geom]{rpoint}}
#'
#' @export
rpoint3 <- function(n, f, fmax = 1, win = box3(), ...,
                    giveup = 1000, verbose = FALSE, nsim = 1, drop = TRUE) {
  if (missing(f) || (is.numeric(f) && length(f) == 1)) {
    return(spatstat.random::runifpoint3(n, domain = win, nsim = nsim, drop = drop))
  }
  if (!is.function(f) && !is.im(f)) {
    stop(paste(sQuote("f"), "must be a function"))
  }
  spatstat.geom::verifyclass(win, "box3")
  if (n == 0) {
    emp <- pp3(numeric(0), numeric(0), numeric(0), window = win)
    if (nsim == 1 && drop) {
      return(emp)
    }
    result <- rep(list(emp), nsim)
    names(result) <- paste("Simulation", 1:nsim)
    return(spatstat.geom::as.anylist(result))
  }
  result <- vector(mode = "list", length = nsim)
  for (isim in 1:nsim) {
    X <- spatstat.geom::pp3(numeric(0), numeric(0), numeric(0), window = win)
    pbar <- 1
    nremaining <- n
    totngen <- 0
    ntries <- 0
    repeat {
      ntries <- ntries + 1
      ngen <- nremaining / pbar + 10
      totngen <- totngen + ngen
      prop <- spatstat.random::runifpoint3(ngen, domain = win)
      if (npoints(prop) > 0) {
        fvalues <- f(prop$data$x, prop$data$y, prop$data$z, ...)
        paccept <- fvalues / fmax
        u <- runif(spatstat.geom::npoints(prop))
        Y <- prop[u < paccept]
        if (spatstat.geom::npoints(Y) > 0) {
          X <- spatstat.geom::superimpose(X, Y, W = win)
          nX <- spatstat.geom::npoints(X)
          pbar <- nX / totngen
          nremaining <- n - nX
          if (nremaining <= 0) {
            if (verbose) {
              spatstat.utils::splat(
                "acceptance rate = ",
                round(100 * pbar, 2), "%"
              )
            }
            result[[isim]] <- if (nX == n) {
              X
            } else {
              X[1:n]
            }
            break
          }
        }
      }
      if (ntries > giveup) {
        stop(paste(
          "Gave up after", giveup * n, "trials with",
          spatstat.geom::npoints(X), "points accepted"
        ))
      }
    }
  }
  if (nsim == 1 && drop) {
    return(result[[1]])
  }
  names(result) <- paste("Simulation", 1:nsim)
  return(spatstat.geom::as.anylist(result))
}

#### rPoissonCluster3 ####
#' Simulate 3D Poisson Cluster Process
#'
#' \code{rPoissonCluster3} extends \code{\link[spatstat.random]{rPoissonCluster}} to
#' \code{\link[spatstat.geom]{pp3}}.
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat.random]{rpoint}}
#'
#'
#' @export
rPoissonCluster3 <- function(kappa, expand, rcluster, win = box3(), ...,
                             nsim = 1, drop = T) {
  if (missing(expand) && !is.null(rmax <- list(...)$rmax)) {
    expand <- rmax
    f <- rcluster
    rcluster <- function(..., rmax) f(...)
  }
  spatstat.geom::verifyclass(win, "box3")
  dilated <- spatstat.geom::box3(
    win$xrange + c(-expand, expand),
    win$yrange + c(-expand, expand),
    win$yrange + c(-expand, expand)
  )
  parentlist <- spatstat.random::rpoispp3(kappa,
    domain = dilated,
    nsim = nsim
  )
  if (nsim == 1) {
    parentlist <- list(parentlist)
  }
  resultlist <- vector(mode = "list", length = nsim)
  for (isim in 1:nsim) {
    parents <- parentlist[[isim]]
    result <- NULL
    res.full <- NULL
    np <- spatstat.geom::npoints(parents)
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
            zrange = win$zrange
          )
          clus.trunc <- subset(cluster,
            subset =
              (x > win$xrange[1] & x < win$xrange[2]) &
                (y > win$yrange[1] & y < win$yrange[2]) &
                (z > win$zrange[1] & z < win$zrange[2])
          )
          if (is.null(result)) {
            result <- clus.trunc
            res.full <- cluster
            parentid <- rep.int(1, npoints(clus.trunc))
          } else {
            result <- superimpose(result, clus.trunc, W = win)
            res.full <- superimpose(res.full, cluster, W = win)
            parentid <- c(parentid, rep.int(i, npoints(clus.trunc)))
          }
        }
      }
    } else {
      result <- pp3(numeric(0), numeric(0), numeric(0), window = win)
      parentid <- integer(0)
    }
    attr(result, "parents") <- parents
    attr(result, "parentid") <- parentid
    attr(result, "expand") <- expand
    attr(result, "full") <- res.full
    resultlist[[isim]] <- result
  }
  if (nsim == 1 && drop) {
    return(resultlist[[1]])
  }
  names(resultlist) <- paste("Simulation", 1:nsim)
  return(as.anylist(resultlist))
}

#### lattice ####
#' Generate a spatial lattice
#'
#' \code{lattice} creates a spatial region filled with lattice points of the
#' specified type. Currently, only simple cubic (`"sc"`), body-centered cubic
#' (`"bcc"`), and face-centered cubic (`"fcc"`) are implemented.
#'
#' @param domain A \code{\link[spatstat.geom]{box3}}. The domain in which to
#' generate the lattice
#' @param a Numeric. The lattice parameter of the lattice. For cubic lattices
#'   (*i.e.* `"sc"`, `"bcc"`, `"fcc"`), this is the length of the side of the
#'   cube; for `"hcp"`, this is the length of the side of the hexagon. See
#'   Details.
#' @param lattice character. The lattice to generate (one of `"sc"`, `"bcc"`,
#'   `"fcc"`, or `"hcp"`).
#' @return A \code{\link[spatstat.geom]{pp3}} with points at the lattice positions
#'
#' @details For a hexagonal close packed (hcp) lattice, there are two dimensions
#'   that are defined: the side of the hexagon, *a*, and the height of the
#'   hexagonal prism, *c*. For a perfect hcp crystal, the ratio of these two
#'   dimensions is exactly \eqn{c/a = \sqrt{8/3}}.
#'
#' @family simulation functions
#' @seealso \code{\link{hcp.gen}}, \code{\link[spatstat.geom]{pp3}}
#'
#' @export
lattice <- function(domain = box3(), a = 1, lattice = "sc") {
  if (lattice == "sc") {
    x <- seq(domain$xrange[1], domain$xrange[2], a)
    y <- seq(domain$yrange[1], domain$yrange[2], a)
    z <- seq(domain$zrange[1], domain$zrange[2], a)
    p <- expand.grid(x = x, y = y, z = x)
    pp3(p$x, p$y, p$z, domain)
  } else if (lattice == "bcc") {
    x <- seq(domain$xrange[1], domain$xrange[2], a)
    y <- seq(domain$yrange[1], domain$yrange[2], a)
    z <- seq(domain$zrange[1], domain$zrange[2], a)
    g1 <- expand.grid(x = x, y = y, z = x)
    x2 <- x + a / 2
    x2 <- x2[x2 <= domain$xrange[2]]
    y2 <- y + a / 2
    y2 <- y2[y2 <= domain$yrange[2]]
    z2 <- z + a / 2
    z2 <- z2[z2 <= domain$zrange[2]]
    g2 <- expand.grid(x = x2, y = y2, z = z2)
    p <- rbind(g1, g2)
    pp3(p$x, p$y, p$z, domain)
  } else if (lattice == "fcc") {
    x <- seq(domain$xrange[1], domain$xrange[2], a / 2)
    y <- seq(domain$yrange[1], domain$yrange[2], a / 2)
    z <- seq(domain$zrange[1], domain$zrange[2], a / 2)
    g <- expand.grid(x = x, y = y, z = x)
    ind <- expand.grid(i = seq_along(x), j = seq_along(y), k = seq_along(z)) - 1
    p <- g[(ind$i + ind$j + ind$k) %% 2 == 0, ]
    pp3(p$x, p$y, p$z, domain)
  } else if (lattice == "hcp") {
    df <- hcp.gen(a / 2, win = domain) # using this function for now; will revise
    pp3(df$x, df$y, df$z, domain)
    # x1 <- seq(domain$xrange[1], domain$xrange[2], a)
  } else {
    warning("Unrecognized lattice")
  }
}

#### nmers ####
#' Select n-mer Clusters
#'
#' Creates an n-mer point pattern by selecting the n+1 nearest neighbors of a
#' sampled point pattern from its parent distribution.
#'
#' @param sam pp3. The seed points for the n-mers
#' @param parent pp3. The point pattern from which to draw the n-mers
#' @param n numeric. The number of points in the n-mer. Default is 2
#' @param full logical. Return a marked pattern with both n-mers and
#'   "background" points? Default is `FALSE`
#'
#' @return A pp3 with the n-mers. If `full = FALSE` (the default) this is
#'   an unmarked pattern with only the positions of the n-mers. If
#'   `full = TRUE`, this is a marked pattern with marks "nmers" and "bkgd".
#'
#' @family simulation functions
#' @seealso \code{\link{clustersim}}
nmers <- function(sam, parent, n = 2, full = FALSE) {
  nmer.ind <- nncross(sam, parent, what = "which", k = 1:n)
  nmer.ind <- unlist(nmer.ind)
  if (full == FALSE) {
    nmer.dat <- parent[nmer.ind]
  } else {
    nmer.dat <- parent
    marks(nmer.dat) <- "bkgd"
    marks(nmer.dat)[nmer.ind] <- "nmer"
  }
  return(nmer.dat)
}
