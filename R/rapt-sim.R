#
# This file contains methods for simulating APT data
#

# Extends the superimpose function from "SpatStat" to handle pp3
superimpose.pp3 <- function(..., W = NULL) {
  input.list <- list(...)
  df.list <- lapply(input.list, as.data.frame)
  df.comb <- Reduce(rbind, df.list)
  out.pp3 <-  createSpat(df.comb, win = W)
  return(out.pp3)
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
# Generalization of the "SpatStat" cluster generators to pp3
# Modeled off of rPoissonCluster
rCluster <- function(kappa, expand, rcluster, win = box3(), ...,
                     lmax = NULL, nsim = 1, drop = T)
{

}
