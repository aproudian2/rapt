# Functions having to do with the calculation of and display of envelopes ran on random re-samplings of pp3 patterns

#### subSquare ####
#' Select a subsection from the center of a \code{\link[spatstat]{pp3}} data
#' file
#'
#' Given an original \code{\link[spatstat]{pp3}} object, \code{subSquare} will
#' select a rectangular prism centered at the center of the original point
#' pattern, and return a \code{\link[spatstat]{pp3}} object of the subsection.
#'
#' @param orig The original \code{\link[spatstat]{pp3}} object
#' @param win Numerical vector containing the dimensions for the box that you
#'   would like to select: c(xdim, ydim, zdim) (e.g. c(10,10,10)).
#' @return Returns a \code{\link[spatstat]{pp3}} object of the selected box,
#'   shifted so that the origin is still at (0,0,0).

subSquare <- function(orig,win){

  orig.domain <- domain(orig)
  orig.center <- c(mean(orig.domain$xrange),mean(orig.domain$yrange),mean(orig.domain$zrange))
  xs <-orig.center[1]-(win[1]/2)
  ys <-orig.center[2]-(win[2]/2)
  zs <-orig.center[2]-(win[3]/2)
  xb <-orig.center[1]+(win[1]/2)
  yb <-orig.center[2]+(win[2]/2)
  zb <-orig.center[3]+(win[3]/2)

  xr <- c(xs,xb)
  yr <- c(ys,yb)
  zr <- c(zs,zb)

  sub.box <- box3(xrange = xr, yrange = yr, zrange = zr)

  tflist <- inside.boxx(orig,w=sub.box)

  sub <- orig[tflist]

  xrn <- c(0,xb-xs)
  yrn <- c(0,yb-ys)
  zrn <- c(0,zb-zs)
  sub.box.new <- box3(xrange = xrn, yrange = yrn, zrange = zrn)

  coo <- coords(sub)

  coo$x <- coo$x - xs
  coo$y <- coo$y - ys
  coo$z <- coo$z - zs

  sub.new <- createSpat(coo,win=sub.box.new)

  return(sub.new)
}

#### percentSelect ####
#' Randomly select a percent of the points in a \code{\link[spatstat]{pp3}}
#' object.
#'
#' Function randomly selects a certain percent of points within an original
#' \code{\link[spatstat]{pp3}} object. This function was created to be used in
#' random relabeling of point patterns.
#'
#' @param perc The fraction of points from the original pattern that are to be
#'   selected. A value between 0 and 1.
#' @param pattern The original \code{\link[spatstat]{pp3}} object to be selected
#'   from.
#' @return A \code{\link[spatstat]{pp3}} object containing only the selected
#'   points.

percentSelect <- function(perc,pattern){

  reLabel <- rlabel(pattern,labels = c(rep("A",round(npoints(pattern)*perc)),rep("B",round(npoints(pattern)*(1-perc)))))
  inds <- which(marks(reLabel)=="A")
  newPattern <- reLabel[inds]
  return(newPattern)
}

#### envPlot ####
#' Plot envelopes of K3est test
#'
#' Plot the results of envelope calculations from the \code{\link{pK3est}} or
#' \code{\link{panomK3est}}, with the ability to choose the percentiles for
#' plotting.
#'
#' @param tests The return file from \code{\link{pK3est}} or the first, [[1]],
#'   entry in the list returned by \code{\link{panomK3est}}. Contains results
#'   from many 3D K tests.
#' @param percentiles Numerical vector of percentiles that you want to see the
#'   envelopes for. Each between 0 and 1.
#' @param ylim Numerical vector containing the min and max values for the y axis
#'   on the plot.
#' @param xlim Numerical vector containing the min and max values for the x axis
#'   on the plot.
#' @return Nothing.

envPlot <- function(tests,percentiles=c(.999,.99,.97),ylim=c(-3,3),xlim=c(0,ceiling(max(tests[,1])))){

  color <- c("lightskyblue","mediumpurple","lightpink")

  # break up data into r values and test results
  rvals <- tests[,1]
  tvals <- tests[,2:ncol(tests)]

  nTests <- ncol(tvals) # number of tests done
  prange <- percentiles*nTests # get the range of indeces for which each percentile spans

  sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
  percentileIndicesBig <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
  percentileIndicesSmall <- round(nTests/2)-floor(prange/2) # do the same for the low end

  # grab out the columns from the sorted test results that we will plot
  toPlotBigs <- matrix(0,nrow=nrow(tvals),ncol=length(percentiles))
  toPlotSmalls <- matrix(0,nrow=nrow(tvals),ncol=length(percentiles))
  for(i in 1:length(percentiles)){
    toPlotBigs[,i] <- sortedtVals[,percentileIndicesBig[i]]
    toPlotSmalls[,i] <- sortedtVals[,percentileIndicesSmall[i]]
  }

  # plot the envelopes from the percentile data
  par(oma = c(0, 2, 0, 0))
  plot(rvals,tvals[,1],type="n",main="Envelopes for K Function",xlab="r",ylab="",ylim=ylim,xlim=xlim)
  mtext(text=expression(sqrt('K'[3]*'(r)')*'  Anomaly'),side=2,line=0,outer=TRUE)
  axis(1,at=0:xlim[2],labels=FALSE)
  axis(1,at=seq(0,xlim[2],by=2))
  for(i in 1:length(percentiles)){
    polygon(c(rvals,rev(rvals)),c(toPlotBigs[,i],rev(toPlotSmalls[,i])),col=color[i])#,border=color[i],lwd=2)
  }
  abline(h=0,lty=2,lwd=1,col="black")
  legend(0, ylim[2], legend=c(paste(toString(percentiles[1]*100),"% AI"), paste(toString(percentiles[2]*100),"% AI"),paste(toString(percentiles[3]*100),"% AI")),col=c(color[1],color[2],color[3]), lty=c(1,1,1), lwd=c(10,10,10))
}

#### pK3est ####
#' Perform K3est on random relabelings in parallel
#'
#' \code{pK3est} first randomly relabels a specified percentage of points from
#' the original \code{\link[spatstat]{pp3}} object. It then performs a 3D K
#' function test (\code{\link[spatstat]{K3est}}) on these selected points. It
#' repeats this as many times as specified. These tests are run in parallel to
#' increase computation speed.
#'
#' @param perc The fraction of points to select randomly each time out of the
#'   original \code{\link[spatstat]{pp3}} object. Number between 0 and 1.
#' @param pattern The original \code{\link[spatstat]{pp3}} object.
#' @param nEvals The number of random relabelings and  that should be performed.
#' @param rmax See \code{\link[spatstat]{K3est}}. Maximum radius to be
#'   calculated for \code{\link[spatstat]{K3est}}.
#' @param nrval See \code{\link[spatstat]{K3est}}. Number of radii that
#'   \code{\link[spatstat]{K3est}} should be calculated at.
#' @param correction Either "iso", "trans", or "bord" edge correction.
#' @section Edge Corrections: See \code{\link[spatstat]{Kest}} or book availible
#'   at \url{http://spatstat.org/book.html} for more info on these edge
#'   corrections.
#'
#' \subsection{Isotropic - "iso"}{Isotropic edge correction. Assumes point
#' pattern is isotropic, or that it can rotate in space without changing
#' statistics.}
#' \subsection{Translation - "trans"}{Translation edge correction.
#' Assumes translation of point pattern does not change statistics.}
#' \subsection{Border - "bord"}{Border edge correction. Makes no assumptions
#' about data. Uses only data provided in the original point pattern. Only
#' evaluates \code{\link[spatstat]{K3est}} when the radius of the search stays
#' within the domain of the point pattern itself.}
#'
#' @return Returns a matrix containing the data from all of the
#'   \code{\link[spatstat]{K3est}} runs on different re-labelings. Can plot data
#'   using \code{\link{envPlot}}.

pK3est <- function(perc, pattern, nEvals,rmax=NULL,nrval=128,correction="iso"){

  #find cores and initialize the cluster
  cores2use <- detectCores()-1
  cl <- makePSOCKcluster(cores2use)
  clusterExport(cl,"percentSelect")
  clusterExport(cl,c("pattern","rmax","nrval","correction"),envir = environment())
  clusterEvalQ(cl,library(spatstat))

  percents <- as.list(rep(perc, nEvals))

  toTest <- parLapply(cl,percents,function(x){
    percentSelect(x,pattern)
  })

  # apply K3est function to each of the pp3 patterns in parallel
  if(correction=="iso"){
    result <- parLapply(cl,toTest,function(x){
      K3est(x,rmax=rmax,nrval=nrval,correction = "isotropic")
    })
  }else if(correction=="trans"){
    result <- parLapply(cl,toTest,function(x){
      K3est(x,rmax=rmax,nrval=nrval,correction = "translation")
    })
  }else if(correction=="bord"){
    clusterExport(cl,"bK3est")
    clusterExport(cl,"bdist.points3")
    result <- parLapply(cl,toTest,function(x){
      bK3est(x,rmax=rmax,nrval=nrval)
    })
  }else{
    print("Please input valid correction argument.")
    return()
  }

  # stop the cluster and revert computer to normal
  stopCluster(cl)

  #fill matrix with results
  tst.length <- length(result[[1]]$r)
  tests <- matrix(0,nrow=tst.length,ncol=(nEvals+1))
  tests[,1] <- result[[1]]$r

  # convert the results into the matrix tests
  for(i in 1:length(result)){
    if(correction=="iso"){
      tests[,(1+i)] <- result[[i]]$iso
    }else if(correction == "trans"){
      tests[,(1+i)] <- result[[i]]$trans
    }else if(correction == "bord"){
      tests[,(1+i)] <- result[[i]]$bord
    }
  }

  return(tests)
}

#### panomK3est ####
#' Perform anomaly K3est envelope calculations
#'
#' See \code{\link{pK3est}}. Performs exactly the same as this function, except
#' returns the "anomaly K3est" results. This means that it returns the square
#' room of the \code{\link[spatstat]{K3est}} results, with the 50th percentile
#' subtracted out. This centers envelopes around zero, and the square root
#' standardized variance across all r values. See book at
#' \url{http://spatstat.org/book.html} for a good statistical reference.
#'
#' @param perc The fraction of points to select randomly each time out of the
#'   original \code{\link[spatstat]{pp3}} object. Number between 0 and 1.
#' @param pattern The original \code{\link[spatstat]{pp3}} object.
#' @param nEvals The number of random relabelings and  that should be performed.
#' @param rmax See \code{\link[spatstat]{K3est}}. Maximum radius to be
#'   calculated for \code{\link[spatstat]{K3est}}.
#' @param nrval See \code{\link[spatstat]{K3est}}. Number of radii that
#'   \code{\link[spatstat]{K3est}} should be calculated at.
#' @param correction Either "iso", "trans", or "bord" edge correction.
#' @section Edge Corrections: See \code{\link[spatstat]{Kest}} or book availible
#'   at \url{http://spatstat.org/book.html} for more info on these edge
#'   corrections.
#'
#' \subsection{Isotropic - "iso"}{Isotropic edge correction. Assumes point
#' pattern is isotropic, or that it can rotate in space without changing
#' statistics.}
#' \subsection{Translation - "trans"}{Translation edge correction.
#' Assumes translation of point pattern does not change statistics.}
#' \subsection{Border - "bord"}{Border edge correction. Makes no assumptions
#' about data. Uses only data provided in the original point pattern. Only
#' evaluates \code{\link[spatstat]{K3est}} when the radius of the search stays
#' within the domain of the point pattern itself.}
#'
#' @param toSub If NULL, use the 50th percentile of the calculated set of
#'   \code{\link[spatstat]{K3est}} envelopes to subtract off. Otherwise, the
#'   second, [[2]], entry in the list returned from this same function. This is
#'   how to compare envelope calculations from different point patterns. You
#'   must subtract the same values from both data sets. toSub allows you to
#'   input the values that were subtracted from a previous set of envelopes, for
#'   comparison.
#'
#' @return A list of: [[1]] Matrix of data for all relabelings. Can be plotted
#'   using \code{\link{envPlot}}. [[2]] Vector containing the values that were
#'   subtracted from the results at each r value. Can be used to subtract from
#'   another set of envelopes for comparison. [[3]] rmax used in the
#'   calculation. [[4]] nrval used in the calculation.

panomK3est <- function(perc, pattern, nEvals,rmax=NULL,nrval=128,correction="iso",toSub=NULL){

  #find cores and initialize the cluster
  cores2use <- detectCores()-1
  cl <- makePSOCKcluster(cores2use)
  clusterExport(cl,"percentSelect")
  clusterExport(cl,c("pattern","rmax","nrval","correction"),envir = environment())
  clusterEvalQ(cl,library(spatstat))

  percents <- as.list(rep(perc, nEvals))

  toTest <- parLapply(cl,percents,function(x){
    percentSelect(x,pattern)
  })

  # apply K3est function to each of the pp3 patterns in parallel
  if(correction=="iso"){
    result <- parLapply(cl,toTest,function(x){
      K3est(x,rmax=rmax,nrval=nrval,correction = "isotropic")
    })
  }else if(correction=="trans"){
    result <- parLapply(cl,toTest,function(x){
      K3est(x,rmax=rmax,nrval=nrval,correction = "translation")
    })
  }else if(correction=="bord"){
    clusterExport(cl,"bK3est")
    clusterExport(cl,"bdist.points3")
    result <- parLapply(cl,toTest,function(x){
      bK3est(x,rmax=rmax,nrval=nrval)
    })
  }else{
    print("Please input valid correction argument.")
    return()
  }

  # stop the cluster and revert computer to normal
  stopCluster(cl)

  #fill matrix with results
  tst.length <- length(result[[1]]$r)
  tests <- matrix(0,nrow=tst.length,ncol=(nEvals+1))
  tests[,1] <- result[[1]]$r

  # convert the results into the matrix tests
  for(i in 1:length(result)){
    if(correction=="iso"){
      tests[,(1+i)] <- result[[i]]$iso
    }else if(correction == "trans"){
      tests[,(1+i)] <- result[[i]]$trans
    }else if(correction == "bord"){
      tests[,(1+i)] <- result[[i]]$bord
    }
  }

  tvals <- tests[,2:ncol(tests)]
  tvals <- sqrt(tvals)
  tvals <- t(apply(tvals,1,sort))

  if(is.null(toSub)){
    if(nEvals%%2==0){
      top <- nEvals/2
      bot <- top+1
      toSub <- (tvals[,top]+tvals[,bot])/2
    }else {
      toSub <- tvals[,(round(nEvals/2))]
    }

    tvals <- apply(tvals,2,function(x){
      x-toSub
    })
  }else{
    tvals <- apply(tvals,2,function(x){
      x-toSub
    })
  }

  tests <- cbind(tests[,1],tvals)

  return(list(tests,toSub,rmax,nrval))
}

#### anomK3est ####
#' Perfrom anomaly K3est on a \code{\link[spatstat]{pp3}} object.
#'
#' See \code{\link[spatstat]{K3est}}. Performs the anomaly K3est on a set of
#' point cloud data. This means taking the square root, and subtracting the 50th
#' percentile from the results. This centers the curve around zero, and
#' standardizeds the variance aat different radii. Used for comparing data to
#' envelopes from \code{\link{panomK3est}}. Will subtract the same values used
#' in the panomK3est test that is being compared to.
#'
#' @param pattern The \code{\link[spatstat]{pp3}} object to analyze.
#' @param result List returned from \code{\link{panomK3est}}.
#' @param correction See \code{\link{panomK3est}} or \code{\link{pK3est}}.
#'
#' @return Returns data fram containing r values and associated anomaly K3est
#'   values.

anomK3est <- function(pattern,result,correction = "iso"){

  if(correction == "iso"){
    a <- K3est(pattern,rmax=result[[3]],nrval=result[[4]],correction="isotropic")
    tvals <- sqrt(a$iso) - result[[2]]
    b <- as.data.frame(cbind(a$r,tvals))
    colnames(b)<-c("r","iso")
    return(b)

  }else if(correction == "trans"){
    a <- K3est(pattern,rmax=result[[3]],nrval=result[[4]],correction="translation")
    tvals <- sqrt(a$trans) - result[[2]]
    b <- as.data.frame(cbind(a$r,tvals))
    colnames(b)<-c("r","trans")
    return(b)

  }else if(correction == "bord"){
    a <- bK3est(pattern,rmax=result[[3]],nrval=result[[4]])
    tvals <- sqrt(a$bord) - result[[2]]
    b <- as.data.frame(cbind(a$r,tvals))
    colnames(b)<-c("r","bord")
    return(b)

  }else if(correction == "all"){
    b <- matrix(0,nrow=result[[4]],ncol=3)

    a <- K3est(pattern,rmax=result[[3]],nrval=result[[4]])
    b[,2] <- sqrt(a$iso) - result[[2]]
    b[,3] <-  sqrt(a$trans) - result[[2]]
    b[,1] <- a$r

    b <- as.data.frame(b)
    colnames(b)<-c("r","iso","trans")

    return(b)
  }else{
    print("Please input valid correction argument.")
    return()
  }
}

#### bK3est ####
#' 3D Border correction for K3est
#'
#' Helper function for \code{\link{panomK3est}}, \code{\link{pK3est}}, and
#' \code{\link{anomK3est}}. This function is a hand written extension of the
#' border correction for 3D point patterns.
#'
#' @param X The point pattern for analysis. \code{\link[spatstat]{pp3}} object.
#' @param rmax See \code{\link[spatstat]{K3est}}. Maximum radius to be
#'   calculated for \code{\link[spatstat]{K3est}}.
#' @param nrval See \code{\link[spatstat]{K3est}}. Number of radii that
#'   \code{\link[spatstat]{K3est}} should be calculated at.
#' @return Border corrected \code{\link[spatstat]{K3est}} data for object X.

bK3est <- function(X,rmax=NULL,nrval=128){

  verifyclass(X,"pp3")

  bi <- bdist.points3(X)
  n <- npoints(X)
  lambda <- n/volume(domain(X))

  if(is.null(rmax)){
    rmax <- max(bi)
  }else if(rmax > max(bi)){
    Print("rmax is too large for data set")
    return()
  }

  cp <- closepairs(X,rmax,twice=FALSE,what="indices")
  cpm <- cbind(cp[[1]],cp[[2]])
  cpm<-cpm[order(cpm[,1]),]
  distmat <- as.matrix(dist(coords(X)))
  cpmdist <- rep(0,nrow(cpm))
  for(i in 1:nrow(cpm)){
    temp <- sort(cpm[i,])
    cpmdist[i] <- distmat[temp[2],temp[1]]
  }

  rlist <- seq(from=0,to=rmax,length.out=nrval)
  Kb <- rep(0,nrval)

  np <- 0
  for(i in 1:n){
    if(bi[i] >= rmax){
      np <- np + 1
    }
  }

  for(j in 1:length(rlist)){
    t <- 0
    r <- rlist[j]
    for(i in 1:nrow(cpm)){
      if(cpmdist[i] <= r){
        if((bi[cpm[i,1]] >= rmax) & (bi[cpm[i,2]] >= rmax)){
          t <- t + 2
        }else if((bi[cpm[i,1]] < rmax) & (bi[cpm[i,2]] < rmax)){
        }else{
          t <- t + 1
        }
      }
    }
    Kb[j] <- t/(lambda*np)
  }

  K <- as.data.frame(cbind(rlist,Kb))
  colnames(K)<-c("r","bord")

  return(K)
}

#### bdist.points3 ####
#' Helper function for border correction \code{\link{bK3est}}.
#'
#' Finds the smallest distance to a boundary for each point in a point pattern.
#'
#' @param X The point pattern for analysis. A \code{\link[spatstat]{pp3}} object.
#' @return An object containing the shortest distance to the boundary for each
#'   point in the pattern X.

bdist.points3 <- function (X) {

  verifyclass(X, "pp3")

  x <- X$data$x
  y <- X$data$y
  z <- X$data$z
  d <- X$domain

  xmin <- min(d$xrange)
  xmax <- max(d$xrange)
  ymin <- min(d$yrange)
  ymax <- max(d$yrange)
  zmin <- min(d$zrange)
  zmax <- max(d$zrange)
  result <- pmin.int(x - xmin, xmax - x, y - ymin, ymax - y , z - zmin , zmax - z)

  return(result)
}


