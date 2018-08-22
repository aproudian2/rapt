# Stitching and re-scaling Functions for Periodic Boundary Conditioned RCP files

#### scale ####
#' Scale RCP data to a desired radius
#'
#' Given an original raw RCP output file, scale the 3D point cloud as if the
#' spheres have a specific radius.
#'
#' Input the original radius, taken from the "system" file output from the
#' RCP generation, and the new desired radius, \code{scale} will scale the
#' RCP point cloud as if the original generation hapened with spheres of
#' radius that you desire.
#'
#' @param pp3file A \code{\link[spatstat]{pp3}} object containing the RCP
#' generated 3D point positions.
#' @param newRadius The radius that you would like to scale to. Default is
#' 0.5 units
#' @param oldRadius The radius of the original RCP generation, taken from
#' the "system" file output from the generation. Required input.
#' @param win A \code{\link[spatstat]{box3}} object indicating the size of
#' the original rcp generation. Can leave blank if these are integer value.
#' @return Will return a \code{\link[spatstat]{pp3}} object with the scaled
#' point positions.
scale <- function(pp3file, newRadius = .5, oldRadius = NULL, win = NULL){
  # pp3file is a pp3 point pattern that you would like to scale
  # newRadius is the radius that you would like your particles to have after scaling
  # oldRadius is the radius that your particles have in the original pp3
  # win in a box3 variable indicating the size of the rcp generation (leave null if integer values, the program will figure it out)

  if(is.null(oldRadius)){
    print("Old radius is required.")
    return(NULL)
  }

  s <- (newRadius/oldRadius)

  if(is.null(win)){
    pp3.domain <- domain(pp3file)
    pp3.box <- box3(xrange=s*round(pp3.domain$xrange),yrange=s*round(pp3.domain$yrange),zrange=s*round(pp3.domain$zrange))
  } else {
    pp3.box <- box3(xrange=s*win$xrange,yrange=s*win$yrange,zrange=s*win$zrange)
  }
  rcp_xyz <- coords(pp3file)
  toReturn <- createSpat(s*rcp_xyz,win=pp3.box)

  return(toReturn)
}

#### stitch ####
#' Stitches together multiple RCP point patterns
#'
#' Takes RCP point patterns with periodic boundary conditions and creates
#' a single, larger point cloud by stitching these patterns together.
#'
#' @param pp3file A \code{\link[spatstat]{pp3}} object containing the RCP
#' generated 3D point positions.
#' @param reps A numerical vector containing the number of repetitions you
#' would like in the x, y and z directions, respectively. c(x,y,z).
#' @param win A \code{\link[spatstat]{box3}} object indicating the size of
#' the original rcp generation. Can leave blank if these are integer value.
#' @return Will return a \code{\link[spatstat]{pp3}} object with the total
#' stiched point pattern.
stitch <- function(pp3file, reps = c(2,2,2), win=NULL){
  # pp3file is a pp3 point pattern with periodic boundary conditions which you want to stitch together
  # bounds is a vector containing the bounds of the box that the point pattern exists in: [xmin, xmax, ymin, ymax, zmin, zmax]
  # reps are the number of iterations you want in each direction c(xreps, yreps, zreps)
  # win in a box3 variable indicating the size of the rcp generation (leave null if integer values, the program will figure it out)

  if(is.null(win)){
    pp3.domain <- domain(pp3file)
    pp3.box <- box3(xrange=reps[1]*pp3.domain$xrange,yrange=reps[2]*pp3.domain$yrange,zrange=reps[3]*pp3.domain$zrange)
  } else {
    pp3.domain <- win
    pp3.box <- box3(xrange=reps[1]*win$xrange,yrange=reps[2]*win$yrange,zrange=reps[3]*win$zrange)
  }

  ogCoords <- coords(pp3file)
  newpp3 <- NULL

  for(i in 0:(reps[1]-1)){
    for(j in 0:(reps[2]-1)){
      for(k in 0:(reps[3]-1)){
        newCoords <- ogCoords
        newCoords$x <- newCoords$x + i*(pp3.domain$xrange[2]-pp3.domain$xrange[1])
        newCoords$y <- newCoords$y + j*(pp3.domain$yrange[2]-pp3.domain$yrange[1])
        newCoords$z <- newCoords$z + k*(pp3.domain$zrange[2]-pp3.domain$zrange[1])

        newpp3 <- rbind(newpp3,newCoords)
      }
    }
  }

  toReturn <- createSpat(newpp3, win=pp3.box)

  return(toReturn)
}

#### nncrossR ####
#' Find Nearest Neighbors within a given radius
#'
#' Finds all nearest neighbors of one pattern within a given input radius
#' of each point in aother pattern. This function is not needed and slow.
#' \code{\link[spatstat]{closepairs.pp3}} will do the same thing but more
#' efficiently.
#'
#' @param X \code{\link[spatstat]{pp3}} object 1. Base points.
#' @param Y \code{\link[spatstat]{pp3}} object 2. Find all points in Y with
#' in a given radius of points in X.
#' @param r Radius to find nearest neighbors within.
#'
#' @return Returns a list of:
#' [[1]] The coordinates of the nearest neighbor points
#' [[2]] The index positions of the nearest neighbors
nncrossR <- function(X, Y, r){

  Y.tf <- TRUE
  X.nn <- vector("numeric")
  i = 1

  while(any(Y.tf)){
    X.nntemp <- vector("numeric")
    a <- nncross(X,Y,k=i)
    Y.tf <- a$dist<r
    if(i == 1){
      Y.coords <- a$which[Y.tf]
      X.nn[Y.tf] <- a$which[Y.tf]
    }else {
      Y.coords <- c(Y.coords,a$which[Y.tf])
      X.nntemp[Y.tf] <- a$which[Y.tf]
      X.nn <- rbind(X.nn,X.nntemp)
    }
    i = i+1
  }
  rownames(X.nn)<-1:length(X.nn[,1])

  return(list(Y.coords,X.nn))
}

#### nnR  ####
#' Find Nearest Neighbors within a given radius
#'
#' Finds all nearest neighbors within a given input radius of each point
#' in the pattern. This function is not needed and slow.
#' \code{\link[spatstat]{closepairs.pp3}} will do the same thing but more
#' efficiently.
#'
#' @param X \code{\link[spatstat]{pp3}} object to find nearest points in.
#' @param r Radius to find nearest neighbors within.
#'
#' @return Returns a list of:
#' [[1]] The nearest neighbor distances.
#' [[2]] The nearest neighbor indixes.
nnR <- function(X,r){
  X.tf <- TRUE
  X.nnw <- vector("numeric")
  X.nnd <- vector("numeric")
  i = 1

  while(any(X.tf)){
    X.nnwtemp <- vector("numeric")
    X.nndtemp <- vector("numeric")
    nnd <- nndist(X,k=i)
    nnw <- nnwhich(X,k=i)
    X.tf <- nnd<r
    if(i==1){
      X.nnw[X.tf] <- nnw[X.tf]
      X.nnd[X.tf] <- nnd[X.tf]
    }else{
      X.nnwtemp[X.tf] <- nnw[X.tf]
      X.nndtemp[X.tf] <- nnd[X.tf]
      X.nnw <- rbind(X.nnw,X.nnwtemp)
      X.nnd <- rbind(X.nnd,X.nndtemp)
    }
    i = i+1
  }
  rownames(X.nnw) <- 1:length(X.nnw[,1])
  rownames(X.nnd) <- 1:length(X.nnd[,1])

  return(list(X.nnd,X.nnw))
}
