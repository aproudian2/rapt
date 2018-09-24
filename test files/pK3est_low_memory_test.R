under <- read.rcp('~/Research/point_patterns/Final/FinalConfig1','~/Research/point_patterns/Final/system1',scaleUp = TRUE,newRadius = 0.5)

t1 <- Sys.time()
a <- pK3est(0.06,under,2000,nrval=200,correction="trans",anom=TRUE)
t2 <- Sys.time()
print(t2-t1)

t1 <- Sys.time()
b <- pK3est.new(0.06,under,2000,nrval=200,correction="trans",anom=TRUE)
t2 <- Sys.time()
print(t2-t1)

envPlot(a[[1]])

pK3est.new <- function(perc, pattern, nEvals,rmax=NULL,nrval=128,correction="iso",anom=FALSE,toSub=NULL){

  #find cores and initialize the cluster
  cores2use <- detectCores()-1
  cl <- makePSOCKcluster(cores2use)
  clusterExport(cl,"percentSelect")
  clusterExport(cl,c("pattern","rmax","nrval","correction"),envir = environment())
  clusterEvalQ(cl,library(spatstat))

  percents <- as.list(rep(perc, nEvals))

  #toTest <- parLapply(cl,percents,function(x){
  #  percentSelect(x,pattern)
  #})

  # apply K3est function to each of the pp3 patterns in parallel
  if(correction=="iso"){
    result <- parLapply(cl,toTest,function(x){
      K3est(x,rmax=rmax,nrval=nrval,correction = "isotropic")
    })
  }else if(correction=="trans"){
    result <- parLapply(cl,percents,function(x){
      K3est(percentSelect(x,pattern),rmax=rmax,nrval=nrval,correction = "translation")
    })
  }else if(correction=="bord"){
    clusterExport(cl,"bK3est")
    clusterExport(cl,"bdist.points3")
    result <- parLapply(cl,toTest,function(x){
      bK3est(x,rmax=rmax,nrval=nrval)
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
