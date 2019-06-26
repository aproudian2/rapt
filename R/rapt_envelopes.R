#
# This file contains information relating to the calculation and display of envelopes.
#

#### envPlot ####
#' Plot envelopes of K3est test
#'
#' Plot the results of envelope calculations from the \code{\link{pK3est}},
#' \code{\link{pG3est}}, or \code{\link{pF3est}} functions, with the ability to
#' choose the percentiles for plotting.
#'
#' @param tests The return file from \code{p(K/G/F)3est} or the first, [[1]],
#'   entry in the list returned by \code{p(K/G/F)3est} with \code{anom = TRUE}.
#' @param percentiles Numerical vector of percentiles that you want to see the
#'   envelopes for. Each between 0 and 1.
#' @param ylim Numerical vector containing the min and max values for the y axis
#'   on the plot.
#' @param xlim Numerical vector containing the min and max values for the x axis
#'   on the plot.
#' @param leg True or falsel whether to show the automatically generated legend.
#' @param colors List of color names to make the envelopes.
#' @param ... Arguments to be passed into \code{plot()}.
#' @return Nothing, just produces a plot.
#' @seealso \code{\link{pK3est}}, \code{\link{pG3est}}, \code{\link{pF3est}}
#' @export

envPlot <- function(tests, percentiles = c(0.999, 0.99, 0.97),
                    ylim = c(-3,3), xlim = c(0, ceiling(max(tests[,1]))),
                    leg = TRUE, ...) {
  color <- c("lightskyblue", "mediumpurple", "lightpink")

  color <- colors
  # break up data into r values and test results
  rvals <- tests[,1]
  tvals <- tests[,2:ncol(tests)]

  nTests <- ncol(tvals) # number of tests done
  prange <- percentiles * nTests # get the range of indeces for which each percentile spans

  sortedtVals <- t(apply(tvals, 1, sort)) # sort the results at each r value from lowest to highest
  percentileIndicesBig <- round(nTests/2) + floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
  percentileIndicesSmall <- round(nTests/2) - floor(prange/2) # do the same for the low end

  # grab out the columns from the sorted test results that we will plot
  toPlotBigs <- matrix(0, nrow = nrow(tvals), ncol = length(percentiles))
  toPlotSmalls <- matrix(0, nrow = nrow(tvals), ncol = length(percentiles))
  for(i in 1:length(percentiles)) {
    toPlotBigs[,i] <- sortedtVals[,percentileIndicesBig[i]]
    toPlotSmalls[,i] <- sortedtVals[,percentileIndicesSmall[i]]
  }

  # plot the envelopes from the percentile data
  #par(oma = c(0, 2, 0, 0))
  toplt <- data.frame(rvals,tvals[,1])
  plot(toplt, type = "n", xlab = "r",
       ylab = expression(sqrt('K'[3]*'(r)')*'  Anomaly'),
       ylim = ylim, xlim = xlim, xaxt = "n", ...)
  axis(1, at = 0:xlim[2], labels=FALSE)
  axis(1, at = seq(0,xlim[2],by=2))
  a <- c(rvals$V1, rev(rvals$V1))
  for(i in 1:length(percentiles)) {
    polygon(a, c(toPlotBigs[,i], rev(toPlotSmalls[,i])), col = color[i])#,border=color[i],lwd=2)
  }
  abline(h = 0, lty = 2, lwd = 1, col="black")
  if(leg == TRUE) {
    legend(0, ylim[2], legend = c(paste(toString(percentiles[1]*100), "% AI"),
                                  paste(toString(percentiles[2]*100), "% AI"),
                                  paste(toString(percentiles[3]*100), "% AI")),
           col=c(color[1], color[2], color[3]),
           lty = c(1,1,1), lwd = c(10,10,10))
  }
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
#' @param anom Whether or not to retun the anomaly results. \code{TRUE} or
#'   \code{FALSE}. See section below for more info.
#' @param toSub The numeric vector of data to subtract for the "anom" pK3est.
#'   Only used when \code{anom = TRUE}. See below for more info.
#' @param sorted Whether to return a sorted table of RRLs (TRUE) or an unsorted
#'   one, where the RRLs are in their original rows.
#' @param cores Number of cores to use for parallel processing. If no parallel
#'   desired, set to 1. If NULL, will automatically set to the number of cores
#'   available on the machine.
#' @section Edge Corrections: See \code{\link[spatstat]{Kest}} or book availible
#'   at \url{http://spatstat.org/book.html} for more info on these edge
#'   corrections.
#'
#'   \subsection{Isotropic - "iso"}{Isotropic edge correction. Assumes point
#'   pattern is isotropic, or that it can rotate in space without changing
#'   statistics.} \subsection{Translation - "trans"}{Translation edge
#'   correction. Assumes translation of point pattern does not change
#'   statistics.} \subsection{Border - "bord"}{Border edge correction. Makes no
#'   assumptions about data. Uses only data provided in the original point
#'   pattern. Only evaluates \code{\link[spatstat]{K3est}} when the radius of
#'   the search stays within the domain of the point pattern itself.}
#'
#' @section Anomaly K3est: When \code{anom = TRUE}, the function returns the
#'   anomaly K3est.This means that it returns the square root of the
#'   \code{\link[spatstat]{K3est}} results, with the 50th percentile subtracted
#'   out. This centers envelopes around zero, and the square root standardized
#'   variance across all r values. See book at
#'   \url{http://spatstat.org/book.html} for a good statistical reference.
#'
#'   \code{toSub} is an argument to be paired with \code{anom = TRUE}. If NULL,
#'   use the 50th percentile of the calculated set of
#'   \code{\link[spatstat]{K3est}} envelopes to subtract off. Otherwise, use the
#'   second, [[2]], entry in the list returned from this same function. This is
#'   how to compare envelope calculations from different point patterns. You
#'   must subtract the same values from both data sets. toSub allows you to
#'   input the values that were subtracted from a previous set of envelopes, for
#'   comparison.
#'
#' @return \subsection{For \code{anom = FALSE}}{Returns a matrix containing the
#' data from the all of the \code{\link[spatstat]{K3est}} runs on different
#' re-labelings. Can plot data using \code{\link{envPlot}}.} \subsection{For
#' \code{anom = TRUE}}{A list of: [[1]] Matrix of data for all relabelings. Can
#' be plotted using \code{\link{envPlot}}. [[2]] Vector containing the values
#' that were subtracted from the results at each r value. Can be used to
#' subtract from another set of envelopes for comparison. [[3]] rmax used in the
#' calculation. [[4]] nrval used in the calculation.}
#' @export

pK3est <- function(perc, pattern, nEvals,rmax=NULL,nrval=128,
                   correction="trans",anom=FALSE,toSub=NULL, sorted=TRUE,
                   cores = NULL){

  #find cores and initialize the cluster
  if(is.null(cores)){
    cores2use <- detectCores()
  }else{
    cores2use <- cores
  }
  cl <- makePSOCKcluster(cores2use)
  clusterExport(cl,c("pattern","rmax","nrval","correction"),envir = environment())
  clusterExport(cl,"percentSelect")
  clusterEvalQ(cl,library(spatstat))

  percents <- as.list(rep(perc, nEvals))

  ### old, heavy memory usage
  #toTest <- parLapply(cl,percents,function(x){
  #  percentSelect(x,pattern)
  #})

  # apply K3est function to each of the pp3 patterns in parallel
  if(correction=="iso"){
    result <- parLapply(cl,percents,function(x){
      K3est(percentSelect(x,pattern),rmax=rmax,nrval=nrval,correction = "isotropic")
    })
  }else if(correction=="trans"){
    result <- parLapply(cl,percents,function(x){
      K3est(percentSelect(x,pattern),rmax=rmax,nrval=nrval,correction = "translation")
    })
  }else if(correction=="bord"){
    clusterExport(cl,"bK3est")
    clusterExport(cl,"bdist.points3")
    result <- parLapply(cl,percents,function(x){
      bK3est(percentSelect(x,pattern),rmax=rmax,nrval=nrval)
    })
    if(is.null(result[[1]])){
      print("rmax is too large for border correction.")
      stopCluster(cl)
      return()
    }
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

  # Convert to anomaly deviation or not
  if (anom == FALSE){
    # If not, just return the regular tests matrix
    return(tests)

  }else if (anom == TRUE){
    # If yes, sort the values, take the sqare root (to keep variance constant over r)
    # then subtract toSub from another pattern or from the 50th perentile of this
    # pattern. Return the results.
    tvals <- tests[,2:ncol(tests)]
    tvals <- sqrt(tvals)

    tvals_sorted <- t(apply(tvals,1,sort))

    #browser()

    if(is.null(toSub)){
      if(nEvals%%2==0){
        top <- nEvals/2
        bot <- top+1
        toSub <- (tvals_sorted[,top]+tvals_sorted[,bot])/2
      }else {
        toSub <- tvals_sorted[,(round(nEvals/2))]
      }
    }

    if(sorted == TRUE){
      tvals <- apply(tvals_sorted,2,function(x){
        x-toSub})
    }else{
      tvals <- apply(tvals,2,function(x){
        x-toSub})
    }

    tests <- cbind(tests[,1],tvals)

    return(list(tests,toSub,rmax,nrval))
  }
}

#### anomK3est ####
#' Perfrom anomaly K3est on a \code{\link[spatstat]{pp3}} object.
#'
#' See \code{\link[spatstat]{K3est}}. Performs the anomaly K3est on a set of
#' point cloud data. This means taking the square root, and subtracting the 50th
#' percentile from the results. This centers the curve around zero, and
#' standardizeds the variance at different radii. Used for comparing data to
#' envelopes from \code{\link{pK3est}} where \code{anom = TRUE}. Will subtract
#' the same values used in the pK3est test that is being compared to.
#'
#' @param pattern The \code{\link[spatstat]{pp3}} object to analyze.
#' @param toSub Returned from \code{\link{pK3est}} with \code{anom = TRUE}.
#'   Second item in the returned list. The data to subtract from the results.
#' @param rmax Max r value. See \code{\link[spatstat]{K3est}}. Should be the
#'   same as the envelopes that you are comparing to.
#' @param nrval Number of r values. See \code{\link[spatstat]{K3est}}. Should be
#'   the same as the envelopes that you are comparing to.
#' @param correction See \code{\link{pK3est}}.
#'
#' @return Returns data fram containing r values and associated anomaly K3est
#'   values.
#'
#' @export

anomK3est <- function(pattern,toSub,rmax,nrval,correction = "trans"){

  if(correction == "iso"){
    a <- K3est(pattern,rmax=rmax,nrval=nrval,correction="isotropic")
    tvals <- sqrt(a$iso) - toSub
    b <- as.data.frame(cbind(a$r,tvals))
    colnames(b)<-c("r","iso")
    return(b)

  }else if(correction == "trans"){
    a <- K3est(pattern,rmax=rmax,nrval=nrval,correction="translation")
    tvals <- sqrt(a$trans) - toSub
    b <- as.data.frame(cbind(a$r,tvals))
    colnames(b)<-c("r","trans")
    return(b)

  }else if(correction == "bord"){
    a <- bK3est(pattern,rmax=rmax,nrval=nrval)
    tvals <- sqrt(a$bord) - toSub
    b <- as.data.frame(cbind(a$r,tvals))
    colnames(b)<-c("r","bord")
    return(b)

  }else if(correction == "all"){
    b <- matrix(0,nrow=nrval,ncol=3)

    a <- K3est(pattern,rmax=rmax,nrval=nrval)
    b[,2] <- sqrt(a$iso) - toSub
    b[,3] <-  sqrt(a$trans) - toSub
    b[,1] <- a$r

    b <- as.data.frame(b)
    colnames(b)<-c("r","iso","trans")

    return(b)
  }else{
    print("Please input valid correction argument.")
    return()
  }
}

#### pG3est ####
#' Perform G3est on random relabelings in parallel
#'
#' \code{pG3est} first randomly relabels a specified percentage of points from
#' the original \code{\link[spatstat]{pp3}} object. It then performs a 3D G
#' function test (\code{\link[spatstat]{G3est}}) on these selected points. It
#' repeats this as many times as specified. These tests are run in parallel to
#' increase computation speed.
#'
#' @param perc The fraction of points to select randomly each time out of the
#'   original \code{\link[spatstat]{pp3}} object. Number between 0 and 1.
#' @param pattern The original \code{\link[spatstat]{pp3}} object.
#' @param nEvals The number of random relabelings and  that should be performed.
#' @param rmax See \code{\link[spatstat]{G3est}}. Maximum radius to be
#'   calculated for \code{\link[spatstat]{G3est}}.
#' @param nrval See \code{\link[spatstat]{G3est}}. Number of radii that
#'   \code{\link[spatstat]{G3est}} should be calculated at.
#' @param correction Either "rs", "km", or "Hanisch" edge correction.
#' @param anom Whether or not to retun the anomaly results. \code{TRUE} or
#'   \code{FALSE}. See section below for more info.
#' @param toSub The numeric vector of data to subtract for the "anom" pG3est.
#'   Only used when \code{anom = TRUE}. See below for more info.
#' @param cores Number of cores to use for parallel processing. If no parallel
#'   desired, set to 1. If NULL, will automatically set to the number of cores
#'   available on the machine.
#' @section Edge Corrections: See  book availible at
#'   \url{http://spatstat.org/book.html} for more info on these edge
#'   corrections.
#'
#' \subsection{Reduced Sample - "rs"}{}
#' \subsection{Kaplan-Meier - "km"}{}
#' \subsection{Hanisch - "Hanisch"}{}
#'
#' @section Anomaly G3est:
#' When \code{anom = TRUE}, the function returns the anomaly
#' \code{\link[spatstat]{G3est}}.This means that it returns the
#' \code{\link[spatstat]{G3est}} results with the 50th percentile subtracted
#' out. This centers envelopes around zero.
#'
#' \code{toSub} is an argumet to be paired with \code{anom = TRUE}. If NULL, use
#' the 50th percentile of the calculated set of \code{\link[spatstat]{G3est}}
#' envelopes to subtract off. Otherwise, use the second, [[2]], entry in the
#' list returned from this same function. This is how to compare envelope
#' calculations from different point patterns. You must subtract the same values
#' from both data sets. toSub allows you to input the values that were
#' subtracted from a previous set of envelopes, for comparison.
#'
#' @return
#' \subsection{For \code{anom = FALSE}}{Returns a matrix containing the data
#' from the all of the \code{\link[spatstat]{G3est}} runs on different
#' re-labelings. Can plot data using \code{\link{envPlot}}.}
#' \subsection{For \code{anom = TRUE}}{A list of: [[1]] Matrix of data for all
#' relabelings. Can be plotted using \code{\link{envPlot}}. [[2]] Vector
#' containing the values that were subtracted from the results at each r value.
#' Can be used to subtract from another set of envelopes for comparison. [[3]]
#' rmax used in the calculation. [[4]] nrval used in the calculation.}
#'
#' @export

pG3est <- function(perc, pattern, nEvals, rmax=NULL, nrval=128,
                   correction="rs", anom=FALSE, toSub=NULL,
                   cores=NULL){

  #find cores and initialize the cluster
  if(is.null(cores)){
    cores2use <- detectCores()
  } else {
    cores2use <- cores
  }
  cl <- makePSOCKcluster(cores2use)
  clusterExport(cl,"percentSelect")
  clusterExport(cl,c("pattern","rmax","nrval","correction"), envir = environment())
  clusterEvalQ(cl,library(spatstat))

  percents <- as.list(rep(perc, nEvals))

  #toTest <- parLapply(cl,percents,function(x){
  #  percentSelect(x,pattern)
  #})

  # apply G3est function to each of the pp3 patterns in parallel
  if(correction=="rs"){
    result <- parLapply(cl,percents,function(x){
      G3est(percentSelect(x,pattern),rmax=rmax,nrval=nrval,correction = "rs")
    })
  }else if(correction=="km"){
    result <- parLapply(cl,percents,function(x){
      G3est(percentSelect(x,pattern),rmax=rmax,nrval=nrval,correction = "km")
    })
  }else if(correction=="Hanisch"){
    result <- parLapply(cl,percents,function(x){
      G3est(percentSelect(x,pattern),rmax=rmax,nrval=nrval,correction = "Hanisch")
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
    if(correction=="rs"){
      tests[,(1+i)] <- result[[i]]$rs
    }else if(correction == "km"){
      tests[,(1+i)] <- result[[i]]$km
    }else if(correction == "Hanisch"){
      tests[,(1+i)] <- result[[i]]$han
    }
  }

  # Convert to anomaly deviation or not
  if (anom == FALSE){
    # If not, just return the regular tests matrix
    return(tests)

  }else if (anom == TRUE){
    # If yes, sort the values, subtract toSub from another pattern or from the 50th perentile of this
    # pattern. Return the results.
    tvals <- tests[,2:ncol(tests)]
    #tvals <- sqrt(tvals)
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
}

#### anomG3est ####
#' Perfrom anomaly G3est on a \code{\link[spatstat]{pp3}} object.
#'
#' See \code{\link[spatstat]{G3est}}. Performs the anomaly G3est on a set of
#' point cloud data. This means subtracting the 50th percentile from the
#' results. This centers the curve around zero. Used for comparing data to
#' envelopes from \code{\link{pG3est}} where \code{anom = TRUE}. Will subtract
#' the same values used in the pG3est test that is being compared to.
#'
#' @param pattern The \code{\link[spatstat]{pp3}} object to analyze.
#' @param toSub Returned from \code{\link{pK3est}} with \code{anom = TRUE}.
#'   Second item in the returned list. The data to subtract from the results.
#' @param rmax Max r value. See \code{\link[spatstat]{G3est}}. Should be the
#'   same as the envelopes that you are comparing to.
#' @param nrval Number of r values. See \code{\link[spatstat]{G3est}}. Should be
#'   the same as the envelopes that you are comparing to.
#' @param correction See \code{\link{pG3est}}. "rs", "km", "Hanisch", or "all".
#'
#' @return Returns data fram containing r values and associated anomaly G3est
#'   values.
#'
#' @export

anomG3est <- function(pattern,toSub,rmax,nrval,correction = "rs"){

  if(correction == "rs"){
    a <- G3est(pattern,rmax=rmax,nrval=nrval,correction="rs")
    tvals <- a$rs - toSub
    b <- as.data.frame(cbind(a$r,tvals))
    colnames(b)<-c("r","rs")
    return(b)

  }else if(correction == "km"){
    a <- G3est(pattern,rmax=rmax,nrval=nrval,correction="km")
    tvals <- a$km - toSub
    b <- as.data.frame(cbind(a$r,tvals))
    colnames(b)<-c("r","km")
    return(b)

  }else if(correction == "Hanisch"){
    a <- G3est(pattern,rmax=rmax,nrval=nrval,correction="Hanisch")
    tvals <- a$han - toSub
    b <- as.data.frame(cbind(a$r,tvals))
    colnames(b)<-c("r","han")
    return(b)

  }else if(correction == "all"){
    b <- matrix(0,nrow=nrval,ncol=4)

    a <- G3est(pattern,rmax=rmax,nrval=nrval)
    b[,2] <- a$rs - toSub
    b[,3] <- a$km - toSub
    b[,4] <- a$han - toSub
    b[,1] <- a$r

    b <- as.data.frame(b)
    colnames(b)<-c("r","rs","km","han")

    return(b)
  }else{
    print("Please input valid correction argument.")
    return()
  }
}


#### pF3est ####
#' Perform F3est on random relabelings in parallel
#'
#' \code{pF3est} first randomly relabels a specified percentage of points from
#' the original \code{\link[spatstat]{pp3}} object. It then performs a 3D F
#' function test (\code{\link[spatstat]{F3est}}) on these selected points. It
#' repeats this as many times as specified. These tests are run in parallel to
#' increase computation speed.
#'
#' @param perc The fraction of points to select randomly each time out of the
#'   original \code{\link[spatstat]{pp3}} object. Number between 0 and 1.
#' @param pattern The original \code{\link[spatstat]{pp3}} object.
#' @param nEvals The number of random relabelings and  that should be performed.
#' @param rmax See \code{\link[spatstat]{F3est}}. Maximum radius to be
#'   calculated for \code{\link[spatstat]{F3est}}.
#' @param nrval See \code{\link[spatstat]{F3est}}. Number of radii that
#'   \code{\link[spatstat]{F3est}} should be calculated at.
#' @param correction Either "rs", "km", or "cs" edge correction.
#' @param anom Whether or not to retun the anomaly results. \code{TRUE} or
#'   \code{FALSE}. See section below for more info.
#' @param toSub The numeric vector of data to subtract for the "anom" pF3est.
#'   Only used when \code{anom = TRUE}. See below for more info.
#' @param cores Number of cores to use for parallel processing. If no parallel
#'   desired, set to 1. If NULL, will automatically set to the number of cores
#'   available on the machine.
#' @section Edge Corrections: See  book availible at
#'   \url{http://spatstat.org/book.html} for more info on these edge
#'   corrections.
#'
#'   \subsection{Reduced Sample - "rs"}{} \subsection{Kaplan-Meier - "km"}{}
#'   \subsection{Chiu-Stoyan (aka Hanisch) - "cs"}{}
#'
#' @section Anomaly F3est: When \code{anom = TRUE}, the function returns the
#'   anomaly \code{\link[spatstat]{F3est}}.This means that it returns the
#'   \code{\link[spatstat]{F3est}} results with the 50th percentile subtracted
#'   out. This centers envelopes around zero.
#'
#'   \code{toSub} is an argumet to be paired with \code{anom = TRUE}. If NULL,
#'   use the 50th percentile of the calculated set of
#'   \code{\link[spatstat]{F3est}} envelopes to subtract off. Otherwise, use the
#'   second, [[2]], entry in the list returned from this same function. This is
#'   how to compare envelope calculations from different point patterns. You
#'   must subtract the same values from both data sets. toSub allows you to
#'   input the values that were subtracted from a previous set of envelopes, for
#'   comparison.
#'
#' @return \subsection{For \code{anom = FALSE}}{Returns a matrix containing the
#' data from the all of the \code{\link[spatstat]{F3est}} runs on different
#' re-labelings. Can plot data using \code{\link{envPlot}}.} \subsection{For
#' \code{anom = TRUE}}{A list of: [[1]] Matrix of data for all relabelings. Can
#' be plotted using \code{\link{envPlot}}. [[2]] Vector containing the values
#' that were subtracted from the results at each r value. Can be used to
#' subtract from another set of envelopes for comparison. [[3]] rmax used in the
#' calculation. [[4]] nrval used in the calculation.}
#'
#' @export

pF3est <- function(perc, pattern, nEvals, rmax=NULL, nrval=128,
                   correction="rs", anom=FALSE, toSub=NULL,
                   cores=NULL){

  #find cores and initialize the cluster
  if(is.null(cores)){
    cores2use <- detectCores()
  } else {
    cores2use <- cores
  }
  cl <- makePSOCKcluster(cores2use)
  clusterExport(cl,"percentSelect")
  clusterExport(cl,c("pattern","rmax","nrval","correction"),envir = environment())
  clusterEvalQ(cl,library(spatstat))

  percents <- as.list(rep(perc, nEvals))

  #toTest <- parLapply(cl,percents,function(x){
  #  percentSelect(x,pattern)
  #})

  # apply F3est function to each of the pp3 patterns in parallel
  if(correction=="rs"){
    result <- parLapply(cl,percents,function(x){
      F3est(percentSelect(x,pattern),rmax=rmax,nrval=nrval,correction = "rs")
    })
  }else if(correction=="km"){
    result <- parLapply(cl,percents,function(x){
      F3est(percentSelect(x,pattern),rmax=rmax,nrval=nrval,correction = "km")
    })
  }else if(correction=="cs"){
    result <- parLapply(cl,percents,function(x){
      F3est(percentSelect(x,pattern),rmax=rmax,nrval=nrval,correction = "cs")
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
    if(correction=="rs"){
      tests[,(1+i)] <- result[[i]]$rs
    }else if(correction == "km"){
      tests[,(1+i)] <- result[[i]]$km
    }else if(correction == "cs"){
      tests[,(1+i)] <- result[[i]]$cs
    }
  }

  # Convert to anomaly deviation or not
  if (anom == FALSE){
    # If not, just return the regular tests matrix
    return(tests)

  }else if (anom == TRUE){
    # If yes, sort the values, subtract toSub from another pattern or from the 50th perentile of this
    # pattern. Return the results.
    tvals <- tests[,2:ncol(tests)]
    #tvals <- sqrt(tvals)
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
}

#### anomF3est ####
#' Perfrom anomaly F3est on a \code{\link[spatstat]{pp3}} object.
#'
#' See \code{\link[spatstat]{F3est}}. Performs the anomaly F3est on a set of
#' point cloud data. This means subtracting the 50th percentile from the
#' results. This centers the curve around zero. Used for comparing data to
#' envelopes from \code{\link{pF3est}} where \code{anom = TRUE}. Will subtract
#' the same values used in the pF3est test that is being compared to.
#'
#' @param pattern The \code{\link[spatstat]{pp3}} object to analyze.
#' @param toSub Returned from \code{\link{pK3est}} with \code{anom = TRUE}.
#'   Second item in the returned list. The data to subtract from the results.
#' @param rmax Max r value. See \code{\link[spatstat]{F3est}}. Should be the
#'   same as the envelopes that you are comparing to.
#' @param nrval Number of r values. See \code{\link[spatstat]{F3est}}. Should be
#'   the same as the envelopes that you are comparing to.
#' @param correction See \code{\link{pF3est}}. "rs", "km", "cs", or "all".
#'
#' @return Returns data fram containing r values and associated anomaly F3est
#'   values.
#'
#' @export

anomF3est <- function(pattern,toSub,rmax,nrval,correction = "rs"){

  if(correction == "rs"){
    a <- F3est(pattern,rmax=rmax,nrval=nrval,correction="rs")
    tvals <- a$rs - toSub
    b <- as.data.frame(cbind(a$r,tvals))
    colnames(b)<-c("r","rs")
    return(b)

  }else if(correction == "km"){
    a <- F3est(pattern,rmax=rmax,nrval=nrval,correction="km")
    tvals <- a$km - toSub
    b <- as.data.frame(cbind(a$r,tvals))
    colnames(b)<-c("r","km")
    return(b)

  }else if(correction == "cs"){
    a <- F3est(pattern,rmax=rmax,nrval=nrval,correction="cs")
    tvals <- a$cs - toSub
    b <- as.data.frame(cbind(a$r,tvals))
    colnames(b)<-c("r","cs")
    return(b)

  }else if(correction == "all"){
    b <- matrix(0,nrow=nrval,ncol=4)

    a <- F3est(pattern,rmax=rmax,nrval=nrval)
    b[,2] <- a$rs - toSub
    b[,3] <- a$km - toSub
    b[,4] <- a$cs - toSub
    b[,1] <- a$r

    b <- as.data.frame(b)
    colnames(b)<-c("r","rs","km","cs")

    return(b)
  }else{
    print("Please input valid correction argument.")
    return()
  }
}

#### bK3est ####
#' 3D Border correction for K3est
#'
#' Helper function for \code{\link{pK3est}}. This function is a hand written
#' extension of the border correction for 3D point patterns.
#'
#' @param X The point pattern for analysis. \code{\link[spatstat]{pp3}} object.
#' @param rmax See \code{\link[spatstat]{K3est}}. Maximum radius to be
#'   calculated for \code{\link[spatstat]{K3est}}.
#' @param nrval See \code{\link[spatstat]{K3est}}. Number of radii that
#'   \code{\link[spatstat]{K3est}} should be calculated at.
#' @return Border corrected \code{\link[spatstat]{K3est}} data for object X.

bK3est <- function(X, rmax=NULL, nrval=128){

  verifyclass(X,"pp3")

  bi <- bdist.points3(X)
  n <- npoints(X)
  lambda <- n/volume(domain(X))

  if(is.null(rmax)){
    rmax <- max(bi)
  }else if(rmax > max(bi)){
    print("rmax is too large for data set")
    return()
  }

  cp <- closepairs(X, rmax, twice=FALSE, what="indices")
  cpm <- cbind(cp[[1]], cp[[2]])
  cpm<-cpm[order(cpm[,1]),]
  distmat <- as.matrix(dist(coords(X)))
  cpmdist <- rep(0,nrow(cpm))
  for(i in 1:nrow(cpm)){
    temp <- sort(cpm[i,])
    cpmdist[i] <- distmat[temp[2], temp[1]]
  }

  rlist <- seq(from=0, to=rmax, length.out=nrval)
  Kb <- rep(0, nrval)

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
