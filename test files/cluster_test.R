#source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/cluster_functions.R')
#source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/rcp_functions.R')
#source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/envelope_functions.R')
#source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/rapt-file.R')
#source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/rapt-extend.R')
#library(spatstat)
#library(rgl)

#r1 <- as.numeric(levels(read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/system_test_1')$V1)[2])
#r2 <- as.numeric(levels(read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/system_test_2')$V1)[2])
#rcp1 <- read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/rcp_test_1',sep=" ",col.names=c("x","y","z","type"))
#rcp1 <- scaleRCP(createSpat(rcp1[,c("x","y","z")]),newRadius = 0.5,oldRadius = r1)
#rcp2 <- read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/rcp_test_2',sep=" ",col.names=c("x","y","z","type"))
#rcp2 <- scaleRCP(createSpat(rcp2[,c("x","y","z")]),newRadius = 0.5,oldRadius = r2)

#pcp <- 0.06
#pic <- 1
#radius1 <- 0.5
#radius2 <- 0.5
#under <- rcp1
#over <- rcp2
#cr <- 2
#toPlot <- TRUE
#showOverPts <- TRUE

makecluster.fast <- function(under,over,radius1,radius2,type = "ppc",ppc=NULL,cr=NULL,d=NULL,pic = 1,pcp = 0.06,toPlot=FALSE,showOverPts=FALSE){

  #real cluster percent
  rcp <- pcp*pic

  under.r <- radius1
  over.r <- radius2
  under.vol <- volume(trueBox(under))

  over.rf <- under.r*cr*((4*pi*npoints(under))/(3*under.vol*rcp))^(1/3)

  over.scaled <- scaleRCP(over,newRadius = over.rf, oldRadius = over.r)
  over.scaledf <- subSample(under,over.scaled)

  cluster.nnR.new <- crosspairs.pp3(over.scaledf,under,cr,what="indices",twice=FALSE,distinct=TRUE,neat=TRUE)
  cluster.ind <- cluster.nnR.new[[2]]
  cluster.info <- factor(cluster.nnR.new[[1]])
  diff <- round(rcp*npoints(under)-length(cluster.ind))
  cluster.adj <- crAdjust.new(cluster.ind,cluster.info,diff,over.scaledf,under)
  cluster.ind <- cluster.adj[[1]]
  cluster.info <- cluster.adj[[2]]

  more <- npoints(under)*pcp-npoints(under)*rcp
  if(more==0){

  }else{
    cluster.ind <- randomInsert(cluster.ind,more,npoints(under))
  }

  cluster.xyz <- coords(under)[cluster.ind,]
  cluster.xyz <- na.omit(cluster.xyz)
  cluster <- createSpat(cluster.xyz)

  if(toPlot==TRUE){
    plot3d.pp3(cluster,col="red",size=5)
    plot3d.pp3(under,col="lightgray",add=TRUE)
    if(showOverPts==TRUE){
      plot3d.pp3(over.scaledf,size= 6,col="black",add=TRUE)
    }
  }

  return(list(cluster,over.scaledf,cluster.info))

}

##################
#X <- over.scaledf
#Y <- under

crAdjust.new <- function(cluster.ind,cluster.info, diff, X, Y){
  # mat is the matrix filled with points and associated cluster index values
  # diff is the difference between the number of values needed and the number in mat
  # X is the over point pattern
  # Y is the under point pattern

  if(diff > 0){
    over.n <- length(levels(cluster.info))
    a <- floor(diff/over.n)
    b <- (diff %% over.n)

    if(a>0){
      nPoints <- as.numeric(summary(cluster.info))
      maxn <- max(nPoints)
      minn <- min(nPoints)
      d <- nncross(X,Y,k = (minn + 1):(maxn + a),what = "which")
      cluster.info <- as.numeric(cluster.info)
      for(i in 1:over.n){
        n <- nPoints[i]
        cluster.ind <- c(cluster.ind,as.numeric(d[i,(n-minn+1):(n-minn+a)]))
        cluster.info <- c(cluster.info,rep(i,length(as.numeric(d[i,(n-minn+1):(n-minn+a)]))))
      }
      cluster.info <- factor(cluster.info)
    }
    if(b > 0){
      nPoints <- as.numeric(summary(cluster.info))[1:b]
      maxn <- max(nPoints)
      minn <- min(nPoints)
      d <- nncross(X,Y,k = (minn + 1):(maxn + 1), what = "which")
      cluster.info <- as.numeric(cluster.info)
      if (b == 1){
        cluster.ind <- c(cluster.ind,as.numeric(d[1]))
        cluster.info <- c(cluster.info, 1)
      }else {
        for(i in 1:b){
          n <- nPoints[i]
          cluster.ind <- c(cluster.ind,as.numeric(d[i,n-minn+1]))
          cluster.info <- c(cluster.info,i)
        }
      }
      cluster.info <- factor(cluster.info)
    }
    return(list(cluster.ind,cluster.info))

  }else if(diff < 0){

    over.n <- length(levels(cluster.info))
    a <- floor((-diff)/over.n)
    b <- ((-diff) %% over.n)

    if(a>0){
      nPoints <- as.numeric(summary(cluster.info))
      cluster.info <- as.numeric(cluster.info)
      for(i in 1:over.n){
        if(i == 1){
          cluster.ind[1:a] <- rep(NaN,a)
          cluster.info[1:a] <- rep(NaN,a)
        }else{
          cluster.ind[(sum(nPoints[1:(i-1)])+1):(sum(nPoints[1:(i-1)])+a)] <- rep(NaN,a)
          cluster.info[(sum(nPoints[1:(i-1)])+1):(sum(nPoints[1:(i-1)])+a)] <- rep(NaN,a)
        }
      }
      cluster.ind <- cluster.ind[!is.nan(cluster.ind)]
      cluster.info <- cluster.info[!is.nan(cluster.info)]
      cluster.info <- factor(cluster.info)
    }
    if(b>0){
      nPoints <- as.numeric(summary(cluster.info))
      cluster.info <- as.numeric(cluster.info)
      for(i in 1:b){
        if(i == 1){
          cluster.ind[1] <- NaN
          cluster.info[1] <- NaN
        }else{
          cluster.ind[(sum(nPoints[1:(i-1)])+1)] <- NaN
          cluster.info[(sum(nPoints[1:(i-1)])+1)] <- NaN
        }
      }
    }
    cluster.ind <- cluster.ind[!is.nan(cluster.ind)]
    cluster.info <- cluster.info[!is.nan(cluster.info)]
    cluster.info <- factor(cluster.info)

    return(list(cluster.ind,cluster.info))

  }else{return(list(cluster.ind,cluster.info))}
}
