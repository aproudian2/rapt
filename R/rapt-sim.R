#
# This file contains methods for simulating APT data
#

# Extends the superimpose function from "SpatStat" to handle pp3
superimpose.pp3 <- function(..., W = NULL, check = F) {
  input.list <- list(...)
  df.list <- lapply(input.list, as.data.frame)
  df.comb <- Reduce(rbind, df.list)
  out.pp3 <-  createSpat(df.comb, win = W)
  return(out.pp3)
}
# Extends the shift function from "SpatStat" to handle pp3
shift.pp3 <- function (X, vec = c(0, 0, 0), ..., origin = NULL)
{
  verifyclass(X, "pp3")
  if (!is.null(origin)) {
    if (!missing(vec))
      warning("argument vec ignored; overruled by argument origin")
    if (is.numeric(origin)) {
      locn <- origin
    }
    else if (is.character(origin)) {
      origin <- pickoption("origin", origin,
                           c(midpoint = "midpoint", bottomleft = "bottomleft"))
      W <- X$domain
      locn <- switch(origin, midpoint = {
        c(mean(W$domain$xrange), mean(W$domain$yrange), mean(W$domain$zrange))
      }, bottomleft = {
        c(W$domain$xrange[1], W$domain$yrange[1], W$domain$zrange[1])
      })
    }
    else stop("origin must be a character string or a numeric vector")
    vec <- -locn
  }
  Y <- pp3(X$data$x + vec[1], X$data$y + vec[2], X$data$z + vec[3],
           xrange = X$domain$xrange + vec[1],
           yrange = X$domain$yrange + vec[2],
           zrange = X$domain$zrange + vec[3])
  attr(Y, "lastshift") <- vec
  return(Y)
}
# Extends the rpoint function from "SpatStat" to handle pp3
rpoint3 <- function (n, f, fmax = 1,  win = box3(), ...,
          giveup = 1000, verbose = FALSE, nsim = 1, drop = TRUE)
{
  if (missing(f) || (is.numeric(f) && length(f) == 1))
    return(runifpoint(n, win, giveup, nsim = nsim, drop = drop))
  if (!is.function(f) && !is.im(f))
    stop(paste(sQuote("f"), "must be  a function"))
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
      if (dim(prop$data)[1] > 0) {
        fvalues <- f(prop$data$x, prop$data$y, prop$data$z, ...)
        paccept <- fvalues/fmax
        u <- runif(dim(prop$data)[1])
        Y <- prop[u < paccept]
        if (dim(Y$data)[1] > 0) {
          X <- superimpose(X, Y, W = win)
          nX <- dim(X$data)[1]
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
                   dim(X$data)[1], "points accepted"))
    }
  }
  if (nsim == 1 && drop)
    return(result[[1]])
  names(result) <- paste("Simulation", 1:nsim)
  return(as.anylist(result))
}
# Extends the rPoissonCluster function from "SpatStat" to pp3
rPoissonCluster3 <- function(kappa, expand, rcluster, win = box3(), ...,
                             nsim = 1, drop = T)
{
  if (missing(expand) && !is.null(rmax <- list(...)$rmax)) {
    expand <- rmax
    f <- rcluster
    rcluster <- function(..., rmax) f(...)
  }
  verifyclass(win, "box3")
  frame <- boundingbox(win)
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
    np <- dim(parents$data)[1]
    if (np > 0) {
      xparent <- parents$data$x
      yparent <- parents$data$y
      zparent <- parents$data$z
      for (i in seq_len(np)) {
        cluster <- rcluster(...)
        if (dim(cluster$data)[1] > 0) {
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
            parentid <- rep.int(1, dim(clus.trunc$data)[1])
          }
          else {
            result <- superimpose(result, clus.trunc, W = win)
            res.full <- superimpose(res.full, cluster, W = win)
            parentid <- c(parentid, rep.int(i, dim(clus.trunc$data)[1]))
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
