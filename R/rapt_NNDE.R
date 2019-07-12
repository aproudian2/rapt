#
# This file contains functions relating to the NNDE density estimation method.
#

#### nndensity.pp3 ####
#' Extension of \code{\link[spatstat]{nndensity}} to handle pp3 objects.
#'
#' Calculates the 3D nearest-neighbor intensity estimate of a point process at
#' either a grid of points or at the point locations in the data set. Utilizes
#' the volume weighted edge correction. See Statistics for Spatial Data by
#' Cressie pg. 654 for more info.
#'
#' @param X The point pattern to estimate the intensity of.
#' @param k Vector containing the nearest-neighbor #s that the estimate should
#'   be calculated for.
#' @param nx,ny,nz If estimating on a grid, the number of grid points in x, y,
#'   and z.
#' @param dz The z spacing for numeric integration for the edge correction
#'   (suggested ~ 0.1-0.5)
#' @param at.points \code{TRUE} or \code{FALSE}. Whether or not to estimate
#'   intensity at points in pattern. If \code{TRUE}, nx, ny, and nz are not
#'   used.
#' @param par \code{TRUE} or \code{FALSE}: whether or not to calculate in
#'   parallel
#' @param cores If \code{par = TRUE}, this is the number of cores to use for the
#'   parallel calculation.
#'
#' @return List containing: [[1]] A data frame of the intensity estimates for
#'   each nearest neighbor value. [[2]] The coordinates of the estimates. [[3]]
#'   The coordinates of the original points from the data set.
#' @export
nndensity.pp3 <- function(X, k, nx, ny, nz, dz,
                          at.points = FALSE, par = TRUE, cores = 7){

  if(at.points == FALSE){
    # set up grid of points and find nearest neighbors from grid to data set
    d <- domain(X)

    xsep <- (d$xrange[2]-d$xrange[1])/nx
    xstart <- d$xrange[1] + xsep/2
    xend <- d$xrange[2] - xsep/2

    ysep <- (d$yrange[2]-d$yrange[1])/ny
    ystart <- d$yrange[1] + ysep/2
    yend <- d$yrange[2] - ysep/2

    zsep <- (d$zrange[2]-d$zrange[1])/nz
    zstart <- d$zrange[1] + zsep/2
    zend <- d$zrange[2] - zsep/2

    x <- seq(xstart,xend, xsep)
    y <- seq(ystart,yend, ysep)
    z <- seq(zstart,zend, zsep)

    coo <- expand.grid(x,y,z)
    names(coo) <- c('x','y','z')
    grid <- pp3(coo$x, coo$y, coo$z, domain(X))

    nnk <- nncross(grid, X, what = "dist", k = k)
    bdist <- bdist.points3.multi(grid) # distances to nearest boundaries
    est.points <- coo

  }else{
    #find nearest neighbors from data set to itself
    nnk <- nndist(X, k = k)
    bdist <- bdist.points3.multi(X) # distances to nearest boundaries
    est.points <- coords(X)
  }

  lambda.est <- local.den.engine(bdist, nnk, k, dz, par, cores)

  res <- list(lambda.est = lambda.est,
              estimate.coords = est.points,
              x = coords(X))

  return(res)
}

#### nncrossden.pp3 ####
#' Similar to \code{\link{nndensity.pp3}}, but specifically for marked patterns
#' with a global intensity inhomogeneity.
#'
#' Calculates the 3D nearest-neighbor intensity estimate of a point process at
#' either a grid of points or at the point locations in the data set. Requires
#' two point patterns: One (\code{X}) that contains just points of the type you
#' want to estimate the intesnity of. A second (\code{Y}) that contains all
#' points from the sample (note that this includes the points from pattern
#' \code{X}) for removing global intensity. Note that this function runs faster
#' on linux and mac machines than windows, but will work on both.
#'
#' @param X The point pattern to estimate the intensity of.
#' @param Y The full point pattern to estimate global intensity.
#' @param k Vector containing the nearest-neighbor #s that the estimate should
#'   be calculated for.
#' @param nx,ny,nz If estimating on a grid, the number of grid points in x, y,
#'   and z.
#' @param at.points \code{TRUE} or \code{FALSE}. Whether or not to estimate
#'   intensity at points in pattern. If \code{TRUE}, nx, ny, and nz are not
#'   used.
#' @param nsplit In this function, large \code{pp3} objects are split into
#'   multiple smaller data sets so that the memory is not overloaded while doing
#'   computations. This parameter is the number of points per split set.
#' @param cores Number of cores to use for parallelization. Set to one for
#'   serial calculation.
#' @param os Either 'windows', 'mac', or 'linux'. Changes the parallelization
#'   method used.
#' @return List containing: [[1]] A data frame of the intensity estimates for
#'   each nearest neighbor value. [[2]] The coordinates of the estimates. [[3]]
#'   The coordinates of the original points from the data set. [[4]] A vector
#'   containing all k values tested.
#' @export
nncrossden.pp3 <- function(X, Y, k, nx, ny, nz,
                           at.points = FALSE, nsplit = 1000,
                           cores = 8, os = "linux") {
  # yes this is complicated as hell... but it works I promise -GV

  t1 <- Sys.time()
  if(at.points == FALSE){
    # set up grid of points
    d <- domain(X)

    xsep <- (d$xrange[2]-d$xrange[1])/nx
    xstart <- d$xrange[1] + xsep/2
    xend <- d$xrange[2] - xsep/2

    ysep <- (d$yrange[2]-d$yrange[1])/ny
    ystart <- d$yrange[1] + ysep/2
    yend <- d$yrange[2] - ysep/2

    zsep <- (d$zrange[2]-d$zrange[1])/nz
    zstart <- d$zrange[1] + zsep/2
    zend <- d$zrange[2] - zsep/2

    x <- seq(xstart,xend, xsep)
    y <- seq(ystart,yend, ysep)
    z <- seq(zstart,zend, zsep)

    coo <- expand.grid(x,y,z)
    names(coo) <- c('x','y','z')

    grid.n <- nrow(coo)
    if(grid.n > nsplit) {
      coo.split.ind <- split(1:nrow(coo), ceiling(1:nrow(coo)/nsplit))
      coo.split <- lapply(coo.split.ind, function(x){coo[x,]})
      grid.split <- lapply(coo.split, function(x){
        pp3(x$x, x$y, x$z, domain(X))
      })
    }else {
      grid.split <- list(pp3(coo$x, coo$y, coo$z, domain(X)))
      coo.split.ind <- list(1:grid.n)
    }

    est.points <- coo

  }else {
    grid.n <- npoints(X)
    if(grid.n > nsplit){
      coo <- coords(X)
      coo.split.ind <- split(1:nrow(coo), ceiling(1:nrow(coo)/nsplit))
      coo.split <- lapply(coo.split.ind, function(x){coo[x,]})
      grid.split <- lapply(coo.split, function(x) {
        pp3(x$x, x$y, x$z, domain(X))
      })
    }else {
      grid.split <- list(X)
      coo.split.ind <- list(1:grid.n)
    }

    est.points <- coords(X)
  }

  lambda.global.Y <- npoints(Y)/volume(domain(Y))

  if(length(k) > 1) {
    lambda.est <- matrix(NA, nrow = grid.n, ncol = length(k))

    if(os == "windows"){
      cl <- makePSOCKcluster(cores)
      clusterExport(cl, "nncross")
      clusterExport(cl, c("X", "k"), envir = environment())
    }else {
      cl <- makeForkCluster(cores)
    }

    nnk.X.split <- parallel::parLapply(cl, grid.split, function(x){
      nncross(x, X, what = "dist", k = k)
    })

    if(os == "windows") {
      clusterExport(cl, c("lambda.global.Y","nnk.X.split"),
                    envir = environment())
    }

    tot <- length(k) * length(grid.split)
    cnt <- 0

    for(i in 1:length(k)){
      for(j in 1:length(grid.split)){
        cp.Y <- crosspairs(grid.split[[j]], Y,
                           rmax = max(nnk.X.split[[j]]), what = "ijd")
        cp.Y <- data.frame(i = cp.Y[[1]], j = cp.Y[[2]], dist = cp.Y[[3]])
        cp.Y.list <- split(cp.Y, factor(cp.Y$i))
        xi <- 1:nrow(nnk.X.split[[j]])
        lambda.est[coo.split.ind[[j]],i] <- parSapply(
          cl, xi, function(xq,i,cp.Y.list,j) {
            (k[i] / sum(cp.Y.list[[xq]]$dist <
                          nnk.X.split[[j]][[xq,i]])) * lambda.global.Y
          },
          i,cp.Y.list,j
        )
        cnt <- cnt + 1
        print(paste(toString(round(100*cnt/tot, 1)), "%", sep = ""))
      }
    }
    stopCluster(cl)

  }else {
    lambda.est <- matrix(NA, nrow = grid.n, ncol = 1)

    if(os == "windows"){
      cl <- makePSOCKcluster(cores)
      clusterExport(cl, "nncross")
      clusterExport(cl, c("X", "k"), envir = environment())
    }else{
      cl <- makeForkCluster(cores)
    }

    nnk.X.split <- parallel::parLapply(cl, grid.split, function(x) {
      nncross(x, X, what = "dist", k = k)
    })

    if(os == "windows"){
      clusterExport(cl,c("lambda.global.Y","nnk.X.split"),
                    envir = environment())
    }

    for(i in 1:length(grid.split)){
      cp.Y <- crosspairs(grid.split[[i]], Y,
                         rmax = max(nnk.X.split[[i]]), what = "ijd")
      cp.Y <- data.frame(i = cp.Y[[1]], j = cp.Y[[2]], dist = cp.Y[[3]])
      cp.Y.list <- split(cp.Y, factor(cp.Y$i))
      xi <- 1:length(nnk.X.split[[i]])

      lambda.est[coo.split.ind[[i]],1] <- parSapply(
        cl, xi, function(xq, i, cp.Y.list) {
          (k/sum(cp.Y.list[[xq]]$dist <
                   nnk.X.split[[i]][[xq]])) * lambda.global.Y},
        i, cp.Y.list
      )

      print(paste(toString(round(100*i/length(grid.split), 1)), "%", sep = ""))
    }
    stopCluster(cl)
  }

  res <- list(lambda.est = lambda.est,
              estimate.coords = est.points,
              x = coords(X))

  res <- as.data.frame(lambda.est)
  names <- sapply(k,function(x){return(paste("nn",toString(x),sep = ""))})
  colnames(res) <- names

  t2 <- Sys.time()
  print("Time to complete:")
  print(t2-t1)
  return(list(lambda.est = res,
              estimate.coords = est.points,
              x = coords(X),
              k = k))
}

#### local.den.onevol ####
#' Helper for \code{\link{local.den.engine}}.
#'
#' Calculates the volume of a sphere that lies inside the domain given its
#' distance to the domain boundary in x, y, z, and its radius.
#'
#' @param x,y,z Shortest distance to boundary in the x, y, and z directions,
#'   respectively.
#' @param r Distance to nearest neighbor of interest (radius of sphere)
#' @param dz z spacing for numeric integration
local.den.onevol <- function(x, y, z, r, dz){

  if(x > r & y > r & z > r){ #If sphere lies fully inside, just return the full sphere volume
    return((4/3)*pi*r^3)
  }

  # If not, calculate the volume inside (SEE MATHEMATICA NOTEBOOK FOR DERIVATIONS)
  Rcalc <- function(z,r){return(sqrt(r^2 - abs(z)^2))} # Plane intersection circle radius at z value from center of sphere of radius r
  Cseg <- function(x,R){return(R^2 * acos(x/R) - x * sqrt(R^2 - x^2))} # Area of circular segment a distance x from center of circle with radius R
  Iseg <- function(x,y,R){return(0.5 * (sqrt((R - x)*(R + x)) - y) * (sqrt((R - y)*(R + y)) - x) -
                                   0.5 * sqrt(R^2 - sqrt((R - x) * (R + x)) * y - sqrt((R - y) * (R + y)) * x) * sqrt(R^2 + sqrt((R - x) * (R + x)) * y + sqrt((R - y) * (R + y)) * x) +
                                   R^2 * acos(sqrt(R^2 + sqrt(R^2 - x^2) * y + sqrt(R^2 - y^2) * x)/(R*sqrt(2))))} # Area of intersection of two circular segments
  Carea <- function(x,y,z,zmax,R){ # The outside area of a circle of radius R at some x, y in the xy plane a distance z from the center of the sphere
    if(z >= zmax){
      return(pi * R^2)
    }else if (sqrt(x^2 + y^2) < R){
      return(Cseg(x,R) + Cseg(y,R) - Iseg(x, y, R))
    }else if(x < R & y < R){
      return(Cseg(x,R) + Cseg(y,R))
    }else if(x < R & y >= R){
      return(Cseg(x,R))
    }else if(x >= R & y < R){
      return(Cseg(y,R))
    }else{
      return(0)
    }}

  # numeric integration
  zseq <- seq(-r,r,by = dz)

  return((4/3)*pi*r^3 - sum(sapply(zseq, function(q){Carea(x,y,q,z,Rcalc(q,r))*dz})))
}

#### local.den.engine ####
#' Find the local intensity density estimate based on nearest neighbors.
#'
#' Helper for \code{\link{nndensity.pp3}}.
#'
#' @param bdist Result from \code{\link{bdist.points3.multi}} giving shortest
#'   distance to boundary in the x, y, and z directions.
#' @param nnk Result from \code{\link[spatstat]{nndist}} or
#'   \code{\link[spatstat]{nncross}} containing nearest neighbor distances.
#' @param k Vector of nn#s to calculate estimate for.
#' @param dz The spacing for numeric integration over the z direction to
#'   determine edge corrections.
#' @param par \code{TRUE} or \code{FALSE}: whether or not to calculate in
#'   parallel
#' @param cores If \code{par = TRUE}, this is the number of cores to use for the
#'   parallel calculation.
#' @return A data.frame with intensity estimates for each value of k (# of
#'   nearest neighbors).
#' @seealso \code{\link{nndensity.pp3}}

local.den.engine <- function(bdist, nnk, k, dz, par = TRUE, cores = 7){
  x <- bdist$x
  y <- bdist$y
  z <- bdist$z

  r <- as.list(nnk)
  ind.x <- 1:length(x)
  ind.k <- 1:length(k)

  if(par == TRUE){
    cl <- makePSOCKcluster(cores)
    clusterExport(cl,"local.den.onevol")
    clusterExport(cl,c("x","y","z","dz","k","r","ind.x","ind.k"), envir = environment())

    lambda.est <- parallel::parLapply(cl, ind.k, function(i.k){
      vols <- sapply(ind.x, function(i.x,rn){local.den.onevol(x[i.x],y[i.x],z[i.x],rn[i.x],dz)}, r[[i.k]])
      l.est <- k[i.k]/vols
      return(l.est)
    })

    stopCluster(cl)

  }else{
    lambda.est <- lapply(ind.k, function(i.k){
      vols <- sapply(ind.x, function(i.x,rn){local.den.onevol(x[i.x],y[i.x],z[i.x],rn[i.x],dz)}, r[[i.k]])
      l.est <- k[i.k]/vols
      return(l.est)
    })
  }

  res <- as.data.frame(lambda.est)
  names <- sapply(k,function(x){return(paste("nn",toString(x),sep = ""))})
  colnames(res) <- names

  return(res)
}
