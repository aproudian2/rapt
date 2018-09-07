# Script for code profiling and figuring out what to make faster.

# Source files from the correct folders. Note that they need to be sourced for profiling purposes.
source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/cluster_functions.R')
source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/rcp_functions.R')
source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/envelope_functions.R')
source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/rapt-file.R')
source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/rapt-extend.R')
source('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/cluster_test.R')
library(spatstat)
library(rgl)
library(profvis)

profvis({

  # Upload RCP patterns
  r1 <- as.numeric(levels(read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/system_test_1')$V1)[2])
  r2 <- as.numeric(levels(read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/system_test_2')$V1)[2])
  rcp1 <- read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/rcp_test_1',sep=" ",col.names=c("x","y","z","type"))
  rcp1 <- scaleRCP(createSpat(rcp1[,c("x","y","z")]),newRadius = 0.5,oldRadius = r1)
  rcp2 <- read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/rcp_test_2',sep=" ",col.names=c("x","y","z","type"))
  rcp2 <- scaleRCP(createSpat(rcp2[,c("x","y","z")]),newRadius = 0.5,oldRadius = r2)


  a <- makecluster(rcp1,rcp2,0.5,0.5,type="ppc",ppc=10,pic=1,pcp=0.06,toPlot=FALSE,showOverPts=FALSE)
  b <- makecluster(rcp1,rcp2,0.5,0.5,type="cr",cr=2,pic=1,pcp=0.06,toPlot=FALSE,showOverPts=FALSE)
  c <- makecluster(rcp1,rcp2,0.5,0.5,type="dist",d=10,pic=1,pcp=0.06,toPlot=FALSE,showOverPts=FALSE)
  d <- makecluster.fast(rcp1,rcp2,0.5,0.5,type="cr",cr=2,pic=1,pcp=0.06,toPlot=TRUE,showOverPts=TRUE)

})
