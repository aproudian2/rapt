# Stitching and re-scaling Functions for Periodic Boundary Conditioned RCP files

#### scaleRCP ####
#' Scale RCP data to a desired radius
#'
#' Given an original raw RCP output file, scale the 3D point cloud as if the
#' spheres have a specific radius.
#'
#' Input the original radius, taken from the "system" file output from the RCP
#' generation, and the new desired radius, \code{scaleRCP} will scale the RCP point
#' cloud as if the original generation hapened with spheres of radius that you
#' desire.
#'
#' @param pp3file A \code{\link[spatstat]{pp3}} object containing the RCP
#'   generated 3D point positions.
#' @param newRadius The radius that you would like to scale to. Default is 0.5
#'   units
#' @param oldRadius The radius of the original RCP generation, taken from the
#'   "system" file output from the generation. Required input.
#' @param win A \code{\link[spatstat]{box3}} object indicating the size of the
#'   original rcp generation. Can leave blank if these are integer value.
#' @return Will return a \code{\link[spatstat]{pp3}} object with the scaled
#'   point positions.
scaleRCP <- function(pp3file, newRadius = .5, oldRadius = NULL, win = NULL){

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
  rcp_xyz <- coords(pp3file)*s
  toReturn <- createSpat(rcp_xyz,win=pp3.box)

  return(toReturn)
}

#### stitch ####
#' Stitches together a RCP point pattern with periodic boundary conditions
#'
#' Takes RCP point patterns with periodic boundary conditions and creates a
#' single, larger point cloud by stitching these patterns together.
#'
#' @param pp3file A \code{\link[spatstat]{pp3}} object containing the RCP
#' generated 3D point positions.
#' @param reps A numerical vector containing the number of repetitions you would
#'   like in the x, y and z directions, respectively. c(x,y,z).
#' @param win A \code{\link[spatstat]{box3}} object indicating the size of the
#'   original rcp generation. Can leave blank if these are integer value.
#'
#' @return Will return a \code{\link[spatstat]{pp3}} object with the total
#'   stiched point pattern.
stitch <- function(pp3file, reps = c(2,2,2), win=NULL){

  if(is.null(win)){
    pp3.domain <- domain(pp3file)
    pp3.box <- box3(xrange=reps[1]*pp3.domain$xrange,yrange=reps[2]*pp3.domain$yrange,zrange=reps[3]*pp3.domain$zrange)
  } else {
    pp3.domain <- win
    pp3.box <- box3(xrange=reps[1]*win$xrange,yrange=reps[2]*win$yrange,zrange=reps[3]*win$zrange)
  }

  ogCoords <- as.matrix(coords(pp3file))
  n <- nrow(ogCoords)
  newpp3 <- matrix(NaN,n*reps[1]*reps[2]*reps[3],3)
  ind <- 0

  for(i in 0:(reps[1]-1)){
    for(j in 0:(reps[2]-1)){
      for(k in 0:(reps[3]-1)){
        newpp3[(ind + 1):(ind + n), 1] <- ogCoords[,1] + i*(pp3.domain$xrange[2]-pp3.domain$xrange[1])
        newpp3[(ind + 1):(ind + n), 2] <- ogCoords[,2] + j*(pp3.domain$yrange[2]-pp3.domain$yrange[1])
        newpp3[(ind + 1):(ind + n), 3] <- ogCoords[,3] + k*(pp3.domain$zrange[2]-pp3.domain$zrange[1])
        ind <- ind + n
      }
    }
  }

  newpp3 <- as.data.frame(newpp3)
  colnames(newpp3) <- c('x','y','z')

  toReturn <- createSpat(newpp3, win=pp3.box)

  return(toReturn)
}

#### stitch.size ###
#' Stitches together a RCP point pattern with periodic boundary conditions
#'
#' Similar to \code{\link{stitch}}. Instead of only inputting the number of
#' repetitions in each dimension, \code{stitch.size} allows you to specify the
#' domain size that you want to return, even if it is not an integer multiple of
#' the original dimensions.
#'
#' @param pp3file A \code{\link[spatstat]{pp3}} object containing the RCP
#'   generated 3D point positions.
#' @param win A \code{\link[spatstat]{box3}} object indicating the size of the
#'   original rcp generation. Can leave blank if these are integer value.
#' @param boxSize A numeric vector of the dimensions that you want in the final
#'   \code{\link[spatstat]{pp3}} object: c(xmax,ymax,zmax). Assumes that
#'   (xmin,ymin,zmin) = (0,0,0).

stitch.size <- function(pp3file, win=NULL, boxSize){

  if(is.null(win)){
    pp3.domain <- domain(pp3file)
    pp3.domain <- box3(xrange=pp3.domain$xrange,yrange=pp3.domain$yrange,zrange=pp3.domain$zrange)
  } else {
    pp3.domain <- win
  }

  reps <- ceiling(c(boxSize[1]/pp3.domain$xrange[2],boxSize[2]/pp3.domain$yrange[2],boxSize[3]/pp3.domain$zrange[2]))
  pp3.box <- box3(xrange=reps[1]*pp3.domain$xrange,yrange=reps[2]*pp3.domain$yrange,zrange=reps[3]*pp3.domain$zrange)

  ogCoords <- as.matrix(coords(pp3file))
  n <- nrow(ogCoords)
  newpp3 <- matrix(NaN,n*reps[1]*reps[2]*reps[3],3)
  ind <- 0

  for(i in 0:(reps[1]-1)){
    for(j in 0:(reps[2]-1)){
      for(k in 0:(reps[3]-1)){
        newpp3[(ind + 1):(ind + n), 1] <- ogCoords[,1] + i*(pp3.domain$xrange[2]-pp3.domain$xrange[1])
        newpp3[(ind + 1):(ind + n), 2] <- ogCoords[,2] + j*(pp3.domain$yrange[2]-pp3.domain$yrange[1])
        newpp3[(ind + 1):(ind + n), 3] <- ogCoords[,3] + k*(pp3.domain$zrange[2]-pp3.domain$zrange[1])
        ind <- ind + n
      }
    }
  }

  newpp3 <- as.data.frame(newpp3)
  colnames(newpp3) <- c('x','y','z')

  toReturn <- createSpat(newpp3, win=pp3.box)

  toReturn2 <- subSquare(toReturn,win=boxSize)

  return(toReturn2)
}

### read.rcp ####
#' Upload an RCP object
#'
#' Function to read in an RCP code output file into R as a pp3 object.
#'
#' @param fpath_config The file path to the RCP FinalConfig file.
#' @param fpath_sys The file path to the associated RCP system file.
#' @param scaleUP Boolean. If \code{TRUE}, scales RCP so particles have new
#'   radius. If \code{FALSE}, RCP stays as generated.
#' @param newRadius If \code{scaleUP = TRUE}, this is the new radius that the
#'   RCP particles will be scaled to. Default is 0.5.
#'
#' @return A \code{\link[spatstat]{pp3}} object of the RCP pattern.

read.rcp <- function(fpath_config, fpath_sys, scaleUp, newRadius=0.5){
  temp_upload <- read.table(fpath_config,sep=" ",col.names=c("x","y","z","type"))

  if(scaleUp == TRUE){
    a <- read.table(fpath_sys)
    r<-as.numeric(levels(a$V1)[2])
    temp <- scaleRCP(createSpat(temp_upload[,c("x","y","z")]),newRadius = newRadius, oldRadius = r)
    return(temp)
  }
  temp <- createSpat(temp_upload[,c("x","y","z")])
  return(temp)
}
