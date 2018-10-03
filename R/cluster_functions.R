# Functions to simulate clusers using RCP data sets

#### makecluster ####
#' Simulate point clustering in an RCP matrix
#'
#' The \code{makecluster} function simulates point clusters using two RCP point
#' clouds. The first point cloud is the "underlaying" pattern. This is the set
#' of points that will be used as actual cluster point locations in the final
#' product. The second point cloud is the "overlaying" pattern. This is the set
#' of points that will determine the positions of the clusters within the
#' underlaying pattern.
#'
#' @param under The underlaying RCP pattern. A \code{\link[spatstat]{pp3}}
#'   object containing the correctly scaled RCP pattern point locations. Should
#'   have used scaleRCP prior to putting hte object into this argument.
#' @param over The overlaying RCP pattern. A \code{\link[spatstat]{pp3}} object
#'   containing the RCP pattern point locations. Scaled equal to the underlaying
#'   pattern.
#' @param radius1 The small radius of the underlaying RCP pattern. Can be found
#'   in the system output file, or as whatever the pattern is scaled to.
#' @param radius2 The small radius of the overlaying RCP pattern. Can be found
#'   in the system output file, or as whatever the pattern is scaled to.
#' @param type How to create the clusters. Either "ppc", "cr", or "dist". See
#' below for more information on each.
#' @param ppc Number of points per cluster if type = "ppc", otherwise NULL.
#' @param cr Cluster radius if type = "cr", otherwise NULL.
#' @param speed "slow", "fast", or "superfast" for type = "cr".
#' @param d Distance between clusters if type = "dist", otherwise NULL.
#' @param pic Percent In Clusters. Percent of the points marked as cluster type
#'   that should actually be contained within the clusters. Number between 0 and
#'   1 For example, if pic = 0.5, 50 percent of the cluster type points will be
#'   in clusters and 50 percent will be randomly spread through the point cloud.
#' @param pcp Percent Cluster Points. Percent of the points in the underlaying
#' RCP pattern that should be marked as cluster type points. Number between 0
#' and 1.
#' @param den Intra cluster density of cluster points. Number between 0 and 1.
#'   Currently written for application to \code{type = "cr"} \code{speed =
#'   "fast"} and \code{"superfast"}, and to \code{type = "dist"}.
#' @param s Seed for the random parts of the cluster generation process.
#' @param toPlot Show a 3D plot of the cluster points once generation is done?
#'   TRUE or FALSE.
#' @param showOverPts If toPlot = TRUE, show all of the points in the RCP
#'   pattern, with clustered points marked in red, or just show the clustered
#'   points? TRUE or FALSE.
#'
#' @section Cluster Creation Methods: What the different arguments for
#'   \code{type} argument mean.
#'   \subsection{Points Per Cluster - "ppc"}{
#'   Specify the number of points that should be in each cluster. Needs to be
#'   paired with the \code{ppc} argument The program will scale the overlaying
#'   pattern so that there are approximately the correct number of clusters
#'   placed through the underlaying patten, mark the n closest points to these
#'   cluster centers as cluster points, where n is the value input to the ppc
#'   argument. The cluster centers are defined by the point locations of the
#'   scaled overlaying RCP pattern. At this point, there should be approximately
#'   N x \code{pcp} x \code{pip} points containined in clusters (N = total
#'   points in the underlaying RCP pattern). The method cleans up by marking or
#'   removing marks around the centers as necesary so that there are exactly N x
#'   pcp x pip points in clusters. The rest of the N x \code{pcp} x
#'   (1-\code{pip}) points are then placed randomly through the remaining
#'   non-cluster-marked points.}
#'   \subsection{Cluster Radius - "cr"}{
#'   Specify the radius of the clusters. Needs to be paried with the \code{cr}
#'   argument. The overlaying RCP pattern will be scaled so that there are
#'   approximately the correct number of clusters placed through the underlaying
#'   pattern. All points in the underlaying pattern within distance \code{cr} of
#'   the cluster centers are marked as cluster points. The cluster centers are
#'   defined by the point locations of the scaled overlaying RCP pattern. At
#'   this point, there should be approximately N x \code{pcp} x \code{pip}
#'   points containined in clusters (N = total points in the underlaying RCP
#'   pattern). The method cleans up by marking or removing marks around the
#'   centers as necesary so that there are exactly N x pcp x pip points in
#'   clusters. The rest of the N x \code{pcp} x (1-\code{pip}) points are then
#'   placed randomly through the remaining non-cluster-marked points.
#'   \code{speed = "slow"} gives the most accurate, but slowest realization.
#'   \code{speed = "fast"} gives the second fastest, but less accurate
#'   realization (fails with large data sets). Clusters may not be 100\% dense
#'   with points. \code{speed = "superfast"} gives the fastest, but least
#'   accurate realization (works well with large data sets). Clusters may not be
#'   100\% dense with points and additional random points may be scattered
#'   through the pattern.}

#'   \subsection{Distance - "dist"}{
#'   Specify the distance between clusters. Needs to be paired with the \code{d}
#'   argument. The overlaying RCP pattern will be scaled according to the
#'   \code{d} argumnet so that the small radius of the overlaying RCP pattern is
#'   \code{d}/2. Selects the closest n points around each cluster center so that
#'   there are N x \code{pcp} x \code{pip} points contained in clusters. The
#'   rest of the N x \code{pcp} x (1-\code{pip}) points are then placed randomly
#'   through the remaining non-cluster-marked points.}
#'
#' @return Returns are different based on \code{type}.
#'   \subsection{\code{type} = "ppc" or "dist"}{List of: [[1]]
#'   \code{\link[spatstat]{pp3}} object containing the cluster marked point
#'   locations. [[2]] \code{\link[spatstat]{pp3}} object containing the
#'   overlaying RCP pattern after scaling. [[3]] Numeric vector containing: [1]
#'   points per cluster 1 [2] number of points with points per cluster 1 [3]
#'   points per cluster 2 [4] number of points with points per cluster 2. }
#'   \subsection{\code{type} = "cr", speed = "slow"}{List of: [[1]]
#'   \code{\link[spatstat]{pp3}} object containing the cluster marked point
#'   locations. [[2]] \code{\link[spatstat]{pp3}} object containing the
#'   overlaying RCP pattern after scaling. [[3]] Numeric vector containing: [1]
#'   number of clusters [2] number of points in those clusters.}
#'   \subsection{\code{type} = "cr", speed = "fast"}{List of: [[1]]
#'   \code{\link[spatstat]{pp3}} object containing the cluster marked point
#'   locations. [[2]] \code{\link[spatstat]{pp3}} object containing the
#'   overlaying RCP pattern after scaling. [[3]] factor containing the number of
#'   points in each cluster.}
#'   \subsection{\code{type} = "cr", speed = "superfast"}{List of: [[1]]
#'   \code{\link[spatstat]{pp3}} object containing the cluster marked point
#'   locations. [[2]] \code{\link[spatstat]{pp3}} object containing the
#'   overlaying RCP pattern after scaling. [[3]] number of points put in or
#'   taken away at random.}

makecluster <- function(under,over,radius1,radius2,
                        type = "ppc",
                        ppc=NULL,
                        cr=NULL,speed = "superfast",
                        d=NULL,
                        pic = 1,
                        pcp = 0.06,
                        den = 1,
                        s = 100,
                        toPlot=FALSE,showOverPts=FALSE){
  ############################################################################################
  # POINTS PER CLUSTER METHOD

  if(type == "ppc"){
    #real cluster percent
    rcp <- pcp*pic

    under.r <- radius1
    over.r <- radius2
    over.rf <- under.r*(ppc/rcp)^(1/3)

    over.scaled <- scaleRCP(over,newRadius = over.rf, oldRadius = over.r,win = domain(over))
    over.scaledf <- subSample(under,over.scaled)

    ppc <- floor(npoints(under)*rcp/npoints(over.scaledf))

    if(ppc == 0){
      print("Points per cluster is too small.")
      return(NULL)
    }

    diff <- npoints(under)*rcp-ppc*npoints(over.scaledf)

    if(diff > 0){
      over.split <- splitpp3(over.scaledf,diff)

      cluster.inddf1 <- nncross(over.split[[1]],under,what="which",k=1:ppc)
      cluster.inddf2 <- nncross(over.split[[2]],under,what="which",k=1:(ppc+1))
      cluster.ind <- NULL

      for(i in 1:npoints(over.split[[1]])){
        if(ppc==1){
          cluster.ind <- c(cluster.ind,cluster.inddf1[i])
        }else{
          cluster.ind <- c(cluster.ind,as.numeric(cluster.inddf1[i,]))
        }
      }
      for(i in 1:npoints(over.split[[2]])){
        cluster.ind <- c(cluster.ind,as.numeric(cluster.inddf2[i,]))
      }
    }else{
      cluster.ind1 <- nncross(over.scaledf,under,what="which",k=1:ppc)
      cluster.ind <- NULL
      for(i in 1:npoints(over.scaledf)){
        cluster.ind <- c(cluster.ind,as.numeric(cluster.ind1[i,]))
      }
    }

    more <- npoints(under)*pcp-npoints(under)*rcp
    if(more==0){

    }else{
      cluster.ind <- randomInsert(cluster.ind,more,npoints(under),s)
    }

    cluster.xyz <- coords(under)[cluster.ind,]
    cluster <- createSpat(cluster.xyz)

    if(toPlot==TRUE){
      plot3d.pp3(cluster,col="red",size=5)
      plot3d.pp3(under,col="lightgray",add=TRUE)
      if(showOverPts==TRUE){
        plot3d.pp3(over.scaledf,size= 6,col="black",add=TRUE)
      }
    }

    return(list(cluster,over.scaledf,c(ppc,npoints(over.scaledf)-diff,ppc+1,diff)))

    ########################################################################################
    # CHOOSE RADIUS METHOD

  }else if(type == "cr"){
    if (speed == "slow"){
      #real cluster percent
      rcp <- pcp*pic

      under.r <- radius1
      over.r <- radius2
      under.vol <- volume(domain(under))

      over.rf <- under.r*cr*((4*pi*npoints(under))/(3*under.vol*rcp))^(1/3)

      over.scaled <- scaleRCP(over,newRadius = over.rf, oldRadius = over.r,win = domain(over))
      over.scaledf <- subSample(under,over.scaled)

      cluster.nnR <- nncrossR(over.scaledf,under,cr)
      cluster.ind <- cluster.nnR[[1]]
      cluster.matrix <- cluster.nnR[[2]]
      diff <- round(rcp*npoints(under)-length(cluster.ind))

      cluster.ind <- crAdjust(cluster.matrix,diff,over.scaledf,under)

      more <- npoints(under)*pcp-npoints(under)*rcp
      if(more==0){

      }else{
        cluster.ind <- randomInsert(cluster.ind,more,npoints(under),s)
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

      return(list(cluster,over.scaledf,c(npoints(over.scaledf),npoints(cluster)-more)))
    }else if (speed == "fast"){
      #real cluster percent
      rcp <- pcp*pic

      under.r <- radius1
      over.r <- radius2
      under.vol <- volume(domain(under))

      over.rf <- under.r*cr*((4*pi*npoints(under))/(3*under.vol*rcp))^(1/3)

      over.scaled <- scaleRCP(over,newRadius = over.rf, oldRadius = over.r,win = domain(over))
      over.scaledf <- subSample(under,over.scaled)

      cluster.nnR.new <- crosspairs.pp3(over.scaledf,under,cr,what="indices",twice=FALSE,distinct=TRUE,neat=TRUE)

      cluster.ind <- cluster.nnR.new[[2]]
      cluster.info <- factor(cluster.nnR.new[[1]])
      diff <- round(rcp*npoints(under)-length(cluster.ind))

      cluster.adj <- crAdjust.new(cluster.ind,cluster.info,diff,over.scaledf,under)
      cluster.ind <- cluster.adj[[1]]
      cluster.info <- as.numeric(cluster.adj[[2]])

      if(den < 1 & den >= 0){
        set.seed(s)
        cluster.ind.split <- split(cluster.ind,cluster.info,drop=FALSE)
        cluster.ind.thinned <- lapply(cluster.ind.split,function(x){
          return(sample(x,round(den*length(x)),replace=FALSE))})
        cluster.ind <- unlist(cluster.ind.thinned)
        cluster.ind.split <- unlist(cluster.ind.split)
      }else{
        cluster.ind.split <- cluster.ind
      }

      more <- npoints(under)*pcp-npoints(under)*rcp*den
      if(more==0){

      }else{
        cluster.ind <- randomInsert(cluster.ind,more,npoints(under),s,cluster.ind.split)
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
    }else if (speed == "superfast"){
      #real cluster percent
      rcp <- pcp*pic

      under.r <- radius1
      over.r <- radius2
      under.vol <- volume(domain(under))

      over.rf <- under.r*cr*((4*pi*npoints(under))/(3*under.vol*rcp))^(1/3)
      #over.rf <- over.rf * .97
      over.scaled <- scaleRCP(over,newRadius = over.rf, oldRadius = over.r,win = domain(over))
      over.scaledf <- subSample(under,over.scaled)

      cluster.nnR.new <- crosspairs.pp3(over.scaledf,under,cr,what="indices",twice=FALSE,distinct=TRUE,neat=TRUE)

      if(den < 1 & den >= 0){
        set.seed(s)
        cluster.ind.split <- split(cluster.nnR.new[[2]],cluster.nnR.new[[1]],drop=FALSE)
        cluster.ind.thinned <- lapply(cluster.ind.split,function(x){
          return(sample(x,round(den*length(x)),replace=FALSE))})
        cluster.ind <- unlist(cluster.ind.thinned)
        cluster.ind.split <- unlist(cluster.ind.split)
      }else {
        cluster.ind <- cluster.nnR.new[[2]]
        cluster.ind.split <- cluster.ind
      }

      #diff <- round(rcp*npoints(under)-length(cluster.ind))

      #if(diff > 0){
      # cluster.ind <- randomInsert(cluster.ind,diff,npoints(under),s,cluster.ind.split)
      #}else if(diff < 0){
      #  cluster.ind <- randomTakeAway(cluster.ind,-1*diff,npoints(under),s)
      #}

      more <- round(npoints(under)*pcp)-length(cluster.ind)
      if(more==0){

      }else if(more > 0){
        cluster.ind <- randomInsert(cluster.ind,more,npoints(under),s,cluster.ind.split)
      }else if(more < 0){
        cluster.ind <- randomTakeAway(cluster.ind,-1*diff,npoints(under),s)
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

      return(list(cluster,over.scaledf,more))
    }else{
      print('Instert speed for cluster radius method. slow, fast, or superfast.')
    }

    ###########################################################################################
    # CHOOSE DISTANCE BETWEEN CLUSTERS METHOD

  }else if(type == "dist"){
    #real cluster percent
    rcp <- pcp*pic

    under.r <- radius1
    over.r <- radius2
    under.vol <- volume(domain(under))

    over.rf <- d/2

    over.scaled <- scaleRCP(over,newRadius = over.rf, oldRadius = over.r, win = domain(over))
    over.scaledf <- subSample(under,over.scaled)

    ppc <- floor(npoints(under)*rcp/npoints(over.scaledf))

    if(ppc == 0){
      print("Distance between clusters is too small")
      return(NULL)
    }

    diff <- npoints(under)*rcp-ppc*npoints(over.scaledf)

    if(diff > 0){
      over.split <- splitpp3(over.scaledf,diff)

      cluster.inddf1 <- nncross(over.split[[1]],under,what="which",k=1:ppc)
      cluster.inddf2 <- nncross(over.split[[2]],under,what="which",k=1:(ppc+1))
      cluster.ind <- vector("numeric",ppc*npoints(over.split[[1]])+(ppc+1)*npoints(over.split[[2]]))
      cnt <- 1
      for(i in 1:npoints(over.split[[1]])){
        if(ppc==1){
          cluster.ind[cnt] <- cluster.inddf1[i]
          cnt <- cnt + 1
        }else{
          cluster.ind[cnt:(cnt+ppc-1)] <- as.numeric(cluster.inddf1[i,])
          cnt <- cnt + ppc
        }
      }
      for(i in 1:npoints(over.split[[2]])){
        cluster.ind[cnt:(cnt+ppc)] <- as.numeric(cluster.inddf2[i,])
        cnt <- cnt + ppc + 1
      }
    }else{
      cluster.ind1 <- nncross(over.scaledf,under,what="which",k=1:ppc)
      cluster.ind <- vector("numeric",ppc*npoints(over.scaledf))
      cnt <- 1
      for(i in 1:npoints(over.scaledf)){
        cluster.ind[cnt:(cnt+ppc-1)] <- as.numeric(cluster.ind1[i,])
        cnt <- cnt + ppc
      }
    }

    if(den < 1 & den >= 0){

      cluster.ind.info <- vector("numeric",length(cluster.ind))
      cnt <- 1
      for(i in 1:npoints(over.split[[1]])){
        cluster.ind.info[cnt:(cnt + ppc -1)] <- rep(i,ppc)
        cnt <- cnt + ppc
      }
      for(i in 1:npoints(over.split[[2]])){
        cluster.ind.info[cnt:(cnt + ppc)] <- rep(i,ppc+1)
        cnt <- cnt + ppc + 1
      }

      set.seed(s)
      cluster.ind.split <- split(cluster.ind,cluster.ind.info,drop=FALSE)
      cluster.ind.thinned <- lapply(cluster.ind.split,function(x){
        return(sample(x,round(den*length(x)),replace=FALSE))})
      cluster.ind <- unlist(cluster.ind.thinned)
      cluster.ind.split <- unlist(cluster.ind.split)
    }

    more <- round(npoints(under)*pcp)-length(cluster.ind)
    if(more==0){
    }else if(more > 0){
      cluster.ind <- randomInsert(cluster.ind,more,npoints(under),s,cluster.ind.split)
    }else if(more < 0){
      cluster.ind <- randomTakeAway(cluster.ind,-1*diff,npoints(under),s)
    }

    cluster.xyz <- coords(under)[cluster.ind,]
    cluster <- createSpat(cluster.xyz)

    if(toPlot==TRUE){
      plot3d.pp3(cluster,col="red",size=5)
      plot3d.pp3(under,col="lightgray",add=TRUE)
      if(showOverPts==TRUE){
        plot3d.pp3(over.scaledf,size= 6,col="black",add=TRUE)
      }
    }

    return(list(cluster,over.scaledf,c(ppc,npoints(over.scaledf)-diff,ppc+1,diff)))

  } else{
    print("Please input a valid type")
    return()
  }
}


#########################################
# Helper functions

#### subSample ####
#' Helper for \code{\link{makecluster}} to cut \code{\link[spatstat]{pp3}} object to size.
#'
#' Takes one \code{\link[spatstat]{pp3}} object, and cuts its volume down to the
#' size of another \code{\link[spatstat]{pp3}} object. Only keeps the points of
#' the first object that lay within the volume of the second object.
#'
#' @param underPattern The second \code{\link[spatstat]{pp3}} object. The volume
#'   that you want to cut down to.
#' @param overPattern The first \code{\link[spatstat]{pp3}} object. The object
#'   that you wish to cut down to the volume of \code{underPattern}.
#' @return A \code{\link[spatstat]{pp3}} object containing \code{overPattern}
#'   cut down to the volume of \code{underPattern}.

subSample <- function(underPattern, overPattern){

  xr <- domain(underPattern)$xrange
  yr <- domain(underPattern)$yrange
  zr <- domain(underPattern)$zrange

  a <- box3(xrange = xr, yrange = yr, zrange = zr)

  tflist <- inside.boxx(overPattern,w=a)
  overPattern.xyz <- coords(overPattern)

  newPattern <- overPattern.xyz[tflist,]
  newPattern <- createSpat(newPattern, win = a)

  return(newPattern)
}

#### splitpp3 ####
#' Helper for \code{\link{makecluster}} that splits a
#' \code{\link[spatstat]{pp3}} into two.
#'
#' Splits a \code{\link[spatstat]{pp3}} function into two smaller subset
#' \code{\link[spatstat]{pp3}} objects, given an input for how many points to
#' put in the first object.
#'
#' @param overPattern The \code{\link[spatstat]{pp3}} object to be split.
#' @param num The number of points from \code{overPattern} to be put into the
#'   first \code{\link[spatstat]{pp3}} object.
#' @return List containing the two \code{\link[spatstat]{pp3}} patterns. If
#'   these two patterns were combined, they would create the original
#'   \code{overPattern} object.

splitpp3 <- function(overPattern, num){
  pat.xyz <- coords(overPattern)

  pat1.xyz <- pat.xyz[1:num,]
  pat2.xyz <- pat.xyz[(num+1):npoints(overPattern),]

  pp3.1 <- createSpat(pat1.xyz)
  pp3.2 <- createSpat(pat2.xyz)

  return(list(pp3.2,pp3.1))
}

#### trueBox ####
#' Helper for \code{\link{makecluster}} that determines a
#' \code{\link[spatstat]{pp3}} object true dimensions.
#'
#' RCP pattern generations, when loaded into R, lose their boundary information.
#' R interprets their boundary to be the extreme point locations in each
#' direction, when really they are usually nice integer values. This function
#' simply rounds the R boundary values to the true ones.
#'
#' @param pp3file The \code{\link[spatstat]{pp3}} object to find the bounds of.
#' @return A \code{\link[spatstat]{as.box3}} object containing the true volume
#'   dimensions.

trueBox <- function(pp3file) {

  xr <- round(domain(pp3file)$xrange)
  yr <- round(domain(pp3file)$yrange)
  zr <- round(domain(pp3file)$zrange)

  a <- box3(xrange = xr, yrange = yr, zrange = zr)

  return(a)
}

#### crAdjust ####
#' Helper for \code{\link{makecluster}} for \code{type} = "cr" and \code{fast} =
#' FALSE
#'
#' Adjustment method to get correct number of points in each cluster for the
#' \code{type} = "cr" cluster generation method of \code{\link{makecluster}}
#' function.
#'
#' @param mat Matrix filled with points and associated cluster index values.
#' @param diff The difference between the number of values needed and the number
#'   in mat
#' @param X The overlaying point pattern. A \code{\link[spatstat]{pp3}} object.
#' @param Y The underlaying point pattern. A \code{\link[spatstat]{pp3}} object.
#'
#' @return Updated mat matrix containing the correct number of cluster points.

crAdjust <- function(mat, diff, X, Y){
  # mat is the matrix filled with points and associated cluster index values
  # diff is the difference between the number of values needed and the number in mat
  # X is the over point pattern
  # Y is the under point pattern

  if(diff > 0){
    over.n <- length(mat[1,])
    a <- floor(diff/over.n)
    b <- (diff %% over.n)

    for(i in 1:(a+1)){
      mat <- rbind(mat,rep(NA,over.n))
    }

    if(a>0){
      for(i in 1:over.n){
        nPoints <- length(mat[,i][!is.na(mat[,i])])
        d <- nncross(X,Y,k = (nPoints+1):(nPoints+a),what="which")
        e <- vector("numeric",a)
        for(j in 1:a){
          e[j] <- d[i,j]
        }
        mat[(nPoints+1):(nPoints+a),i] <- e
      }
    }
    if(b>0){
      for(i in 1:b){
        nPoints <- length(mat[,i][!is.na(mat[,i])])
        d <- nncross(X,Y,k = (nPoints+1),what="which")
        mat[(nPoints+1),i] <- d[i]
      }
    }

    cluster.ind <- NULL

    for(i in 1:over.n){
      nPoints <- length(mat[,i][!is.na(mat[,i])])
      cluster.ind <- c(cluster.ind,mat[1:nPoints,i])
    }

    return(cluster.ind)

  }else if(diff < 0){

    over.n <- length(mat[1,])
    a <- floor((-diff)/over.n)
    b <- ((-diff) %% over.n)

    if(a>0){
      for(i in 1:over.n){
        nPoints <- length(mat[,i][!is.na(mat[,i])])
        mat[(nPoints-a+1):nPoints,i] <- rep(NA,a)
      }
    }
    if(b>0){
      for(i in 1:b){
        nPoints <- length(mat[,i][!is.na(mat[,i])])
        mat[nPoints,i] <- NA
      }
    }

    cluster.ind <- NULL

    for(i in 1:over.n){
      nPoints <- length(mat[,i][!is.na(mat[,i])])
      cluster.ind <- c(cluster.ind,mat[1:nPoints,i])
    }

    return(cluster.ind)

  }else{return(mat)}
}

#### crAdjust.new ####
#' Helper for \code{\link{makecluster}} for \code{type} = "cr" and \code{fast} =
#' TRUE
#'
#' Adjustment method to get correct number of points in each cluster for the
#' \code{type} = "cr" cluster generation method of \code{\link{makecluster}}
#' function. This is a version of \code{\link{crAdjust}} meant to work with the
#' faster, updated verios of \code{type}="cr" \code{\link{makecluster}}.
#'
#' @param cluster.ind Vector of indices of the cluster points.
#' @param cluster.info Factor containing the number of points in each cluster.
#' @param diff The difference between the number of values needed and the number
#'   in cluster.ind.
#' @param X The overlaying point pattern. A \code{\link[spatstat]{pp3}} object.
#' @param Y The underlaying point pattern. A \code{\link[spatstat]{pp3}} object.
#'
#' @return A list of [[1]] The updated cluster.ind vector and [[2]] the updated
#'   cluster.info factor.

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


#### randomInsert ####
#' Helper for \code{\link{makecluster}} to insert random cluster points
#'
#' When \code{pip} argument of \code{\link{makecluster}} is not equal to 1, or
#' \code{superfast == TRUE}, there may be random points that need to be marked
#' as cluster type placed within the underlaying pattern. This function does
#' just that.
#'
#' @param cluster.indices A vector containing the indices of the current cluster
#'   points
#' @param n The number of points that need to be placed randomly
#' @param N The number of points in the entire underlaying pattern.
#' @param s The seed for the random inserted points.
#' @param points.avoid A vector containing the indices of the points you don't
#'   want the function to insert points to.
#' @return New indices vector containing new cluster points randomly placed.

# function to randomly place points within the under data set, if not 100% of the cluster points are set to be in the clusters
randomInsert <- function(cluster.indices,n,N,s,points.avoid = cluster.indices){
  #cluster.Indices is a vector containig the indices of the current cluster points
  #n is the number of points that need to be placed randomly
  #N is the number of points in the entire underlying pattern

  full <- 1:N
  nonclust.ind <- full[!(full %in% points.avoid)]

  set.seed(s)
  inds <- sample(nonclust.ind,n)

  all.ind <- c(cluster.indices,inds)

  return(all.ind)
}

#### randomTakeAway ####
#' Helper for \code{\link{makecluster}} to take away random cluster points
#'
#' When \code{superfast} argument of \code{\link{makecluster}} is \code{TRUE},
#' there may be random points that need to be removed from the existing clusters
#' of the the underlaying pattern. This function does just that.
#'
#' @param cluster.indices A vector containing the indices of the current cluster
#'   points
#' @param n The number of points that need to be removed randomly
#' @param N The number of points in the entire underlaying pattern.
#' @param s The random seed for the random selection of points to take away.
#' @return New indices vector containing new cluster points randomly placed.

# function to randomly place points within the under data set, if not 100% of the cluster points are set to be in the clusters
randomTakeAway <- function(cluster.indices,n,N,s){
  #cluster.Indices is a vector containig the indices of the current cluster points
  #n is the number of points that need to be placed randomly
  #N is the number of points in the entire underlying pattern

  set.seed(s)
  inds <- sample(1:length(cluster.indices),n)
  cluster.indices[inds] <- NaN

  all.ind <- cluster.indices[!is.nan(cluster.indices)]

  return(all.ind)
}
