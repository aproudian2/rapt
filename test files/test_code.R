# Test Script!
# Run these lines of code to check that the outputs are what they should be after changing code

# Source files from the correct folders. Note that they need to be sourced for profiling purposes.
source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/cluster_functions.R')
source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/rcp_functions.R')
source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/envelope_functions.R')
source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/rapt-file.R')
source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/rapt-extend.R')
library(spatstat)
library(rgl)

##############################
# makecluster

# Upload RCP patterns
r1 <- as.numeric(levels(read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/system_test_1')$V1)[2])
r2 <- as.numeric(levels(read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/system_test_2')$V1)[2])
rcp1 <- read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/rcp_test_1',sep=" ",col.names=c("x","y","z","type"))
rcp1 <- scaleRCP(createSpat(rcp1[,c("x","y","z")]),newRadius = 0.5,oldRadius = r1)
rcp2 <- read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/rcp_test_2',sep=" ",col.names=c("x","y","z","type"))
rcp2 <- scaleRCP(createSpat(rcp2[,c("x","y","z")]),newRadius = 0.5,oldRadius = r2)

# make a cluster of type "ppc"
ppc_cluster <- makecluster(rcp1,rcp2,0.5,0.5,type="ppc",ppc=10,pic=1,pcp=0.06,toPlot=TRUE,showOverPts=TRUE)
# check out the plots. Make sure they look like they're clustered.
# ppc_cluster[[3]] should read: 9 10 10 39
ppc_cluster[[3]]
# other tests you could run:
  # ppc = 2 returns 1 24 2 228
  # ppc = 3 returns 2 21 3 146
  # ppc = 5 returns 4 20 5 80
  # ppc = 20 returns 20 24 21 0

# make a cluster of type "cr"
cr_cluster <- makecluster(rcp1,rcp2,0.5,0.5,type="cr",cr=2,pic=1,pcp=0.06,toPlot=TRUE,showOverPts=TRUE)
# check out the plots. Make sure they look like they're clustered.
# cr_cluster[[3]] should read: 14 480
cr_cluster[[3]]
# other tests you could run:
  # cr = 1 returns 121 480
  # cr = 1.25 returns 57 480
  # cr = 1.1 returns 88 480

# make a cluster of type "dist"
dist_cluster <- makecluster(rcp1,rcp2,0.5,0.5,type="dist",d=10,pic=1,pcp=0.06,toPlot=TRUE,showOverPts=TRUE)
# check out the plots. Make sure they look like they're clustered.
# dist_cluster[[3]] should read 68 3 69 4
dist_cluster[[3]]
# other tests you could run:
  # d = 5 returns 7 56 8 11
  # d = 3 returns 1 146 2 167
  # d = 7 returns 20 24 21 0

# make a cluster of type "ppc" fail
makecluster(rcp1,rcp2,0.5,0.5,type="ppc",ppc=1,pic=1,pcp=0.06,toPlot=TRUE,showOverPts=TRUE)
# should return "Points per cluster is too small".

# make a cluster of type "cr" fail
makecluster(rcp1,rcp2,0.5,0.5,type="cr",cr=5,pic=1,pcp=0.06,toPlot=TRUE,showOverPts=TRUE)
# should return "Error in X.nn[,1] : incorrect number of dimensions"

# make a cluster of type "dist" fail
makecluster(rcp1,rcp2,0.5,0.5,type="dist",d = 2,pic=1,pcp=0.06,toPlot=TRUE,showOverPts=TRUE)
# should return "Distance between clusters is too small"
makecluster(rcp1,rcp2,0.5,0.5,type="dist",d = 20,pic=1,pcp=0.06,toPlot=TRUE,showOverPts=TRUE)
# should return "Error in if (diff > 0) { : missing value where TRUE/FALSE needed}"

# Test out the percent in clusters variable. Check that 50% (and others) in clusters works.

# make a cluster of type "ppc" with pic = 0.5
ppc_cluster <- makecluster(rcp1,rcp2,0.5,0.5,type="ppc",ppc=10,pic=0.5,pcp=0.06,toPlot=TRUE,showOverPts=TRUE)
# check out the plot. Should look like clusters with points randomly dispursed also.
# ppc_cluster[[3]] should read 10 24 11 0
ppc_cluster[[3]]
# other tests you could run:
  # ppc = 10, pic = 0.25 returns 9 10 10 3
  # ppc = 3, pic = 0.5 returns 2 9 3 74
  # ppc = 5, pic = 0.8 returns  4 21 5 60
  # ppc = 2, pic = 0.25 returns 2 60 3 0

# make a cluster of type "cr" with pic = 0.05
cr_cluster <- makecluster(rcp1,rcp2,0.5,0.5,type="cr",cr=2,pic=0.5,pcp=0.06,toPlot=TRUE,showOverPts=TRUE)
# check out the plot. Should look like clusters with points also randomly dispursed.
# cr_cluster[[3]] should read 5 240
cr_cluster[[3]]
# other tests you could run:
  # cr = 1, pic = 0.5 returns 56 240
  # cr = 1.1, pic = 0.25 returns 21 120
  # cr = 1.8, pic = 0.8 returns 15 384

# make a cluster of type "dist" with pic = 0.05
dist_cluster <- makecluster(rcp1,rcp2,0.5,0.5,type="dist",d=10,pic=0.5,pcp=0.06,toPlot=TRUE,showOverPts=TRUE)
# plots should be good.
# dist_cluster[[3]] should read 34 5 35 2
# other tests you could run:
  # d = 5, pic = 0.5 returns 3 28 4 39
  # d = 3, pic = 0.25 returns "Distance between clusters is too small"
  # d = 7, pic = 0.25 returns 5 24 6 0
  # d = 3, pic = 0.8 returns 1 242 2 71


