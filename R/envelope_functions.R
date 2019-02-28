# Functions having to do with the calculation of and display of envelopes ran
# on random re-samplings of pp3 patterns

# Function to grab a cubic subsection from the center of the given data, given
# dimensions
subSquare <- function(orig, win) {
 # win = list of dimensions of the box you want to pull: c(xdim,ydim,zdim)
 # orig = the original point pattern (pp3) - needs to be scaled in the same way
 # that your window is scaled
  orig.domain <- domain(orig)
  orig.center <- c(mean(orig.domain$xrange),
                   mean(orig.domain$yrange),
                   mean(orig.domain$zrange))
  xs <-orig.center[1] - (win[1]/2)
  ys <-orig.center[2] - (win[2]/2)
  zs <-orig.center[2] - (win[3]/2)
  xb <-orig.center[1] + (win[1]/2)
  yb <-orig.center[2] + (win[2]/2)
  zb <-orig.center[3] + (win[3]/2)

  xr <- c(xs, xb)
  yr <- c(ys, yb)
  zr <- c(zs, zb)

  sub.box <- box3(xrange = xr, yrange = yr, zrange = zr)

  tflist <- inside.boxx(orig, w = sub.box)

  sub <- orig[tflist]

  xrn <- c(0, xb-xs)
  yrn <- c(0, yb-ys)
  zrn <- c(0, zb-zs)
  sub.box.new <- box3(xrange = xrn, yrange = yrn, zrange = zrn)

  coo <- coords(sub)

  coo$x <- coo$x - xs
  coo$y <- coo$y - ys
  coo$z <- coo$z - zs

  sub.new <- createSpat(coo,win=sub.box.new)

  return(sub.new)
}

######################################################################

# Function to randomly relabel and select perc percent of the given point
# pattern
percentSelect <- function(perc, pattern) {
  # perc = fraction of points you want to select (0-1)
  # pattern = the point pattern you want to draw from
  reLabel <- rlabel(pattern,
                    labels = c(rep("A", round(npoints(pattern) * perc)),
                               rep("B", round(npoints(pattern) * (1-perc))))
                    )
  inds <- which(marks(reLabel) == "A")
  newPattern <- reLabel[inds]
  return(newPattern)
}

######################################################################

# Function to plot envelope results from the pK3est or panomK3est functions
# below, based on percentiles
envPlot <- function(tests, percentiles = c(0.999, 0.99, 0.97),
                    ylim = c(-3, 3), xlim = c(0, ceiling(max(tests[,1])))
                    ) {
  # tests = array of values returned from the rrK3est function above,
  # percentiles = vector including the different percentiles you would like to
  #               see on the plot (0-1)
  # do these in descending order please

  color <- c("lightskyblue", "mediumpurple", "lightpink")

  # break up data into r values and test results
  rvals <- tests[,1]
  tvals <- tests[,2:ncol(tests)]

  nTests <- ncol(tvals) # number of tests done
  prange <- percentiles * nTests # get the range of indeces for which each percentile spans

  sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
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
  par(oma = c(0, 2, 0, 0))
  plot(rvals, tvals[,1],
       type = "n", main = "Envelopes for K Function", xlab = "r", ylab = "",
       ylim = ylim, xlim = xlim)
  mtext(text = expression(sqrt('K'[3]*'(r)')*'  Anomaly'),
        side = 2, line = 0, outer = TRUE)
  axis(1 ,at = 0:xlim[2], labels = FALSE)
  axis(1, at = seq(0, xlim[2], by = 2))
  for(i in 1:length(percentiles)) {
    polygon(c(rvals, rev(rvals)), c(toPlotBigs[,i], rev(toPlotSmalls[,i])),
            col = color[i])#,border=color[i],lwd=2)
  }
  abline(h = 0, lty = 2, lwd = 1, col = "black")
  legend(0, ylim[2], legend = c(paste(toString(percentiles[1]*100),"% AI"),
                                paste(toString(percentiles[2]*100),"% AI"),
                                paste(toString(percentiles[3]*100),"% AI")),
         col = c(color[1], color[2], color[3]),
         lty = c(1, 1, 1), lwd = c(10, 10, 10))
}

# Parellelized version of rrK3est, written above. Perfeorms the K3est function
# calculations in parallel
pK3est <- function(perc, pattern, nEvals,rmax=NULL,nrval=128,correction="iso"){
  # perc = percent of original pattern to sample
  # pattern = the original pattern
  # nEvals = number of times to run K3est

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

# Same as the parrrK3est function above, but returns envelopes centered around zero
panomK3est <- function(perc, pattern, nEvals,rmax=NULL,nrval=128,correction="iso",toSub=NULL){
  # perc = percent of original pattern to sample
  # pattern = the original pattern
  # nEvals = number of times to run K3est

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

# Normal K3 anomoly function (no envelopes)
# Need to pass in the result of an envelope generation, to know what to subtract from the result

anomK3est <- function(pattern,result,correction = "iso") {
  #pattern is the pp3 pattern you would like to run a test on
  #result is the list returned by "panomK3est", which holds the information needed to perform the same test
  #correction is the type of edge correction to use - "iso", "trans", "bord", or "all"

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

#######################################
# Border correction for K3est

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

#############################################
# Find distance to boundary for each point in pattern
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
