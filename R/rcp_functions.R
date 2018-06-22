# Stitching and re-scaling Functions for Periodic Boundary Conditioned RCP files


# Function to scale an RCP file to a desired radius, given the original radius
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

#####################################

# Function to stitch an RCP pattern with periodic boundary conditions together
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

#############################################################################

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

######################################################################################

# Finds the nearest neighbors within a point pattern within a given radius of each point
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
