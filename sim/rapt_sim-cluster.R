# Functions to simulate clusters using RCP data sets

#' Simulate Clusters
#'
#' \code{makeCluster} uses two point patterns to generate clusters. One is the
#' underlying point pattern, from which the actual simulated data will be taken.
#' The other is to randomly place the cluster locations through the underlying
#' pattern.
#'
#' @param under A pp3. The underlying point pattern
#' @param over A pp3. The pattern of clusters
#' @param radius1 A numeric. The small radius of the rcpUnder pattern
#' @param radius2 A numeric. The small radius of the rcpOver pattern
#' @param type A string. One of "ppc" (points per cluster), "cr" (cluster
#' radius), or "dist" (distance between points). See Details.
#' @param ppc A numeric. The points per cluster that you want for "ppc." This is
#' a SUGGESTION; it will get close.
#' @param cr A numeric. The radius of cluster that you want for "cr."
#' @param d A numeric. The distance between clusters that you want for "dist."
#' @param pic A numeric. The % of the points designated as cluster type points
#' that you would like to actually be in clusters.
#' @param pcp A numeric. The % of points out of the underlying pattern that you
#' want to be included as points of the clustering type.
#'
makeCluster <- function(under, over, radius1, radius2, type = "ppc",
                        ppc = NULL, cr = NULL, d = NULL,
                        pic = 1, pcp = 0.06,
                        toPlot = FALSE, showOverPts = FALSE) {

  #### POINTS PER CLUSTER (ppc) METHOD ####
  if(type == "ppc") {
    #real cluster percent
    rcp <- pcp * pic

    under.r <- radius1
    over.r <- radius2
    over.rf <- under.r*(ppc/rcp)^(1/3)

    over.scaled <- scale(over, newRadius = over.rf, oldRadius = over.r)
    over.scaledf <- subSample(under, over.scaled)

    ppc <- floor(npoints(under) * rcp / npoints(over.scaledf))

    if(ppc == 0) {
      print("Points per cluster is too small.")
      return(NULL)
    }

    diff <- npoints(under) * rcp-ppc * npoints(over.scaledf)

    if(diff > 0){
      over.split <- splitpp3(over.scaledf, diff)

      cluster.inddf1 <- nncross(over.split[[1]], under, what = "which",
                                k = 1:ppc)
      cluster.inddf2 <- nncross(over.split[[2]], under, what = "which",
                                k = 1:(ppc+1))
      cluster.ind <- NULL

      for(i in 1:npoints(over.split[[1]])) { #clean for
        if(ppc == 1) {
          cluster.ind <- c(cluster.ind, cluster.inddf1[i])
        } else {
          cluster.ind <- c(cluster.ind, as.numeric(cluster.inddf1[i,]))
        }
      }
      for(i in 1:npoints(over.split[[2]])){
        cluster.ind <- c(cluster.ind, as.numeric(cluster.inddf2[i,]))
      }
    } else {
      cluster.ind1 <- nncross(over.scaledf, under, what = "which", k = 1:ppc)
      cluster.ind <- NULL
      for(i in 1:npoints(over.scaledf)) {
        cluster.ind <- c(cluster.ind, as.numeric(cluster.ind1[i,]))
      }
    }

    more <- npoints(under) * pcp-npoints(under) * rcp
    if(more == 0) {

    } else {
      cluster.ind <- randomInsert(cluster.ind, more, npoints(under))
    }

    cluster.xyz <- coords(under)[cluster.ind,]
    cluster <- createSpat(cluster.xyz)

    if(toPlot == TRUE) {
      plot3d.pp3(cluster, col = "red", size = 5)
      plot3d.pp3(under, col = "lightgray", add = TRUE)
      if(showOverPts == TRUE) {
        plot3d.pp3(over.scaledf, size = 6, col = "black", add = TRUE)
      }
    }

    return(list(cluster, over.scaledf,
                c(ppc, npoints(over.scaledf) - diff, ppc + 1, diff)
                ))

  #### CHOOSE RADIUS METHOD ####

  } else if(type == "cr") {
    #real cluster percent
    rcp <- pcp * pic

    under.r <- radius1
    over.r <- radius2
    under.vol <- volume(trueBox(under))

    over.rf <- under.r * cr * ((4 * pi * npoints(under)) /
                                 (3 * under.vol * rcp))^(1/3)

    over.scaled <- scale(over, newRadius = over.rf, oldRadius = over.r)
    over.scaledf <- subSample(under, over.scaled)

    cluster.nnR <- nncrossR(over.scaledf, under,cr)
    cluster.ind <- cluster.nnR[[1]]
    cluster.matrix <- cluster.nnR[[2]]
    diff <- round(rcp * npoints(under) - length(cluster.ind))

    cluster.ind <- crAdjust(cluster.matrix, diff, over.scaledf, under)

    more <- npoints(under) * pcp - npoints(under) * rcp
    if(more == 0) {

    } else{
      cluster.ind <- randomInsert(cluster.ind, more, npoints(under))
    }

    cluster.xyz <- coords(under)[cluster.ind,]
    cluster.xyz <- na.omit(cluster.xyz)
    cluster <- createSpat(cluster.xyz)

    if(toPlot == TRUE) {
      plot3d.pp3(cluster, col = "red", size = 5)
      plot3d.pp3(under, col = "lightgray", add = TRUE)
      if(showOverPts == TRUE) {
        plot3d.pp3(over.scaledf, size = 6,col = "black", add = TRUE)
      }
    }

    return(list(cluster, over.scaledf))

  #### CHOOSE DISTANCE BETWEEN CLUSTERS METHOD ####
  } else if(type == "dist") {
    #real cluster percent
    rcp <- pcp * pic

    under.r <- radius1
    over.r <- radius2
    under.vol <- volume(trueBox(under))

    over.rf <- d / 2

    over.scaled <- scale(over, newRadius = over.rf, oldRadius = over.r)
    over.scaledf <- subSample(under, over.scaled)

    ppc <- floor(npoints(under) * rcp / npoints(over.scaledf))

    if(ppc == 0) {
      print("Distance between clusters is too small")
      return(NULL)
    }

    diff <- npoints(under) * rcp - ppc * npoints(over.scaledf)

    if(diff > 0) {
      over.split <- splitpp3(over.scaledf, diff)

      cluster.inddf1 <- nncross(over.split[[1]], under, what = "which",
                                k = 1:ppc)
      cluster.inddf2 <- nncross(over.split[[2]], under, what = "which",
                                k = 1:(ppc+1))
      cluster.ind <- NULL

      for(i in 1:npoints(over.split[[1]])) {
        if(ppc == 1) {
          cluster.ind <- c(cluster.ind, cluster.inddf1[i])
        } else {
          cluster.ind <- c(cluster.ind, as.numeric(cluster.inddf1[i,]))
        }
      }
      for(i in 1:npoints(over.split[[2]])) {
        cluster.ind <- c(cluster.ind, as.numeric(cluster.inddf2[i,]))
      }
    } else {
      cluster.ind1 <- nncross(over.scaledf, under, what = "which", k = 1:ppc)
      cluster.ind <- NULL
      for(i in 1:npoints(over.scaledf)) {
        cluster.ind <- c(cluster.ind, as.numeric(cluster.ind1[i,]))
      }
    }

    more <- npoints(under) * pcp - npoints(under) * rcp
    if(more ==0) {

    } else {
      cluster.ind <- randomInsert(cluster.ind, more, npoints(under))
    }

    cluster.xyz <- coords(under)[cluster.ind,]
    cluster <- createSpat(cluster.xyz)

    if(toPlot == TRUE) {
      plot3d.pp3(cluster, col = "red", size = 5)
      plot3d.pp3(under, col = "lightgray", add = TRUE)
      if(showOverPts == TRUE) {
        plot3d.pp3(over.scaledf, size = 6, col = "black", add = TRUE)
      }
    }

    return(list(cluster, over.scaledf,
                c(ppc, npoints(over.scaledf) - diff, ppc+1, diff)
                ))

  } else {
    print("Please input a valid type")
    return()
  }
}

#### Helper functions ####

# This is duplicate. Remove.
# Helper function to subset the points of one pattern which lay within the
# window of another pattern
subSample <- function(underPattern, overPattern) {

  xr <- round(domain(underPattern)$xrange)
  yr <- round(domain(underPattern)$yrange)
  zr <- round(domain(underPattern)$zrange)

  a <- box3(xrange = xr, yrange = yr, zrange = zr)

  tflist <- inside.boxx(overPattern, w = a)
  overPattern.xyz <- coords(overPattern)

  newPattern <- overPattern.xyz[tflist,]
  newPattern <- createSpat(newPattern, win = a)

  return(newPattern)
}

# Examine.
# function to split a pp3 into two smaller pp3s, depending on how many clusters
# need an extra point for the "ppc" method
splitpp3 <- function(overPattern, num) {
  pat.xyz <- coords(overPattern)

  pat1.xyz <- pat.xyz[1:num,]
  pat2.xyz <- pat.xyz[(num+1):npoints(overPattern),]

  pp3.1 <- createSpat(pat1.xyz)
  pp3.2 <- createSpat(pat2.xyz)

  return(list(pp3.2, pp3.1))
}

# Duplicate. Remove.
#function to get the true box dimensions of an rcp generation from the original file
trueBox <- function(pp3file) {

  xr <- round(domain(pp3file)$xrange)
  yr <- round(domain(pp3file)$yrange)
  zr <- round(domain(pp3file)$zrange)

  a <- box3(xrange = xr, yrange = yr, zrange = zr)

  return(a)
}

#' Correct points per cluster
#' \code{crAdjust} adjusts the number of points in each cluster defined by the
#' "cr" method to get the correct number of points in each cluster.
#' @param mat A matrix. The matrix filled with points and associated cluster
#' index values.
#' @param diff A numeric. The difference between the number of values needed
#' and the number in \code{mat}.
#' @param X A \code{pp3}. The over point pattern.
#' @param Y A \code{pp3}. The under point pattern.
crAdjust <- function(mat, diff, X, Y) {

  if(diff > 0) {
    over.n <- length(mat[1,])
    a <- floor(diff / over.n)
    b <- (diff %% over.n)

    for(i in 1:(a+1)) {
      mat <- rbind(mat, rep(NA, over.n))
    }

    if(a > 0) {
      for(i in 1:over.n) {
        nPoints <- length(mat[,i][!is.na(mat[,i])])
        d <- nncross(X, Y, k = (nPoints+1):(nPoints+a), what = "which")
        e <- vector("numeric", a)
        for(j in 1:a) {
          e[j] <- d[i,j]
        }
        mat[(nPoints+1):(nPoints+a),i] <- e
      }
    }
    if(b > 0) {
      for(i in 1:b) {
        nPoints <- length(mat[,i][!is.na(mat[,i])])
        d <- nncross(X, Y, k = (nPoints+1), what = "which")
        mat[(nPoints+1),i] <- d[i]
      }
    }

    cluster.ind <- NULL

    for(i in 1:over.n) {
      nPoints <- length(mat[,i][!is.na(mat[,i])])
      cluster.ind <- c(cluster.ind, mat[1:nPoints,i])
    }

    return(cluster.ind)

  } else if(diff < 0) {

    over.n <- length(mat[1,])
    a <- floor((-diff)/over.n)
    b <- ((-diff) %% over.n)

    if(a>0) {
      for(i in 1:over.n) {
        nPoints <- length(mat[,i][!is.na(mat[,i])])
        mat[(nPoints-a+1):nPoints,i] <- rep(NA, a)
      }
    }
    if(b > 0) {
      for(i in 1:b) {
        nPoints <- length(mat[,i][!is.na(mat[,i])])
        mat[nPoints,i] <- NA
      }
    }

    cluster.ind <- NULL

    for(i in 1:over.n) {
      nPoints <- length(mat[,i][!is.na(mat[,i])])
      cluster.ind <- c(cluster.ind, mat[1:nPoints,i])
    }

    return(cluster.ind)

  } else{return(mat)}
}

#' Create cluster random background
#'
#' \code{randomInsert} randomly places points within the under data set, if not
#' 100% of the cluster points are set to be in the clusters
#'
#' @param cluster.Indices A numeric vector. A vector containing the indices of
#' the current cluster points.
#' @param n A numeric. The number of points that need to be placed randomly.
#' @param N A numeric. The number of points in the entire underlying pattern.
#' @return A numeric vector containing the indicies of the background points.
randomInsert <- function(cluster.Indices, n, N) {
  full <- 1:N
  nonclust.ind <- full[!full %in% cluster.Indices]
  inds <- sample(nonclust.ind, n)
  all.ind <- c(cluster.Indices, inds)
  return(all.ind)
}