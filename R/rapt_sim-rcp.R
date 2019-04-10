# Stitching and re-scaling Functions for Periodic Boundary Conditioned RCP files

#' Scale an RCP pattern
#'
#' \code{scale} scales an RCP file to a desired radius, given the original
#' radius.
#' @param pp3file A \code{pp3}. A point pattern to scale.
#' @param newRadius A numeric. The radius that you would like your particles
#' to have after scaling.
#' @param oldRadius A numeric. The radius that your particles have in the
#' original \code{pp3}.
#' @param win A \code{box3}. A domain indicating the size of the rcp
#' generation (leave null if integer values, the program will figure it out).
scale <- function(pp3file, newRadius = 0.5, oldRadius = NULL, win = NULL) {

  if(is.null(oldRadius)) {
    print("Old radius is required.")
    return(NULL)
  }

  s <- (newRadius / oldRadius)

  if(is.null(win)) {
    pp3.domain <- domain(pp3file)
    pp3.box <- box3(xrange = s * round(pp3.domain$xrange),
                    yrange = s * round(pp3.domain$yrange),
                    zrange = s * round(pp3.domain$zrange))
  } else {
    pp3.box <- box3(xrange = s * win$xrange,
                    yrange = s * win$yrange,
                    zrange = s * win$zrange)
  }
  rcp_xyz <- coords(pp3file)
  toReturn <- createSpat(s * rcp_xyz, win = pp3.box)
  return(toReturn)
}

#' Stitch a \code{pp3} with periodic boundary conditions.
#'
#' \code{stitch} takes a \code{pp3} generated with periodic boundary conditions
#' (usually read in from an external RCP generator) and tiles it together.
#'
#' @param pp3file A \code{pp3}. A point pattern with periodic boundary
#' conditions which you want to stitch together.
#' @param reps A numeric vector. The number of iterations you want in each
#' direction: c(xreps, yreps, zreps).
#' @param win A \code{box3}. The size of the rcp generation (leave null if
#' integer values, the program will figure it out).
stitch <- function(pp3file, reps = c(2,2,2), win = NULL) {
  if(is.null(win)) {
    pp3.domain <- domain(pp3file)
    pp3.box <- box3(xrange = reps[1] * pp3.domain$xrange,
                    yrange = reps[2] * pp3.domain$yrange,
                    zrange = reps[3] * pp3.domain$zrange)
  } else {
    pp3.domain <- win
    pp3.box <- box3(xrange = reps[1] * win$xrange,
                    yrange = reps[2] * win$yrange,
                    zrange = reps[3] * win$zrange)
  }

  ogCoords <- coords(pp3file)
  newpp3 <- NULL

  for(i in 0:(reps[1]-1)) {
    for(j in 0:(reps[2]-1)) {
      for(k in 0:(reps[3]-1)) {
        newCoords <- ogCoords
        newCoords$x <- newCoords$x +
          i * (pp3.domain$xrange[2] - pp3.domain$xrange[1])
        newCoords$y <- newCoords$y +
          j * (pp3.domain$yrange[2] - pp3.domain$yrange[1])
        newCoords$z <- newCoords$z +
          k * (pp3.domain$zrange[2] - pp3.domain$zrange[1])

        newpp3 <- rbind(newpp3, newCoords)
      }
    }
  }
  toReturn <- createSpat(newpp3, win = pp3.box)
  return(toReturn)
}

# Duplicate. Remove.
#Function to find the indices of the nearest neighbors within a certain radius
#Finds the nearest neighbors in Y to each point in X within a radius of r (points can be shared)
# Now realize that these functions below are the same thing as closepairs.pp3... oh well
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

# Duplicate. Remove.
# Finds the nearest neighbors within a point pattern within a given radius of each point
nnR <- function(X,r) {
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
