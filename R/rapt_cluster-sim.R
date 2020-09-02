#
# This file contains functions pertaining to the simulation of clustered point
# patterns.
#

#### clustersim ####
#' Simulate clusters of marks on an RCP background point pattern.
#'
#' The \code{clustersim} function simulates point clusters using two RCP point
#' clouds. The first point cloud is the "underlaying" pattern. This is the set
#' of points that will be used as actual cluster point locations in the final
#' product. The second point cloud is the "overlaying" pattern. This is the set
#' of points that will be scaled to determine the positions of the cluster
#' centroids within the underlying point pattern.
#'
#' @param under The underlying RCP point pattern (UPP). A
#'   \code{\link[spatstat]{pp3}} object containing points at the centers of RCP
#'   spheres.
#' @param over The overlying RCP pattern. A \code{\link[spatstat]{pp3}} object
#'   containing points at the centers of RCP spheres. Scaled equal to the
#'   underlying point pattern.
#' @param rcp_rad The radius of the spheres within the underlying and overlying
#'   RCP patterns.
#' @param pcp The fraction of type-A (clustering) points within the UPP. A
#'   decimal value between 0 and 1.
#' @param cr Mean cluster radius. Any value larger than zero.
#' @param rho1 Intra-cluster type-A point concentration. A decimal value between
#'   0 and 1.
#' @param rho2 Background type-A point concentration. A decimal value between 0
#'   and the value of \code{pcp}.
#' @param rb Radius blur. A decimal value between 0 and 1.
#' @param pb Position blur. A decimal value between 0 and 1.
#' @param tol Tolerance value for \code{pcp}. The tru fraction of type-A points
#'   in the pattern will be within this tolerance of the value specified by the
#'   \code{pcp} parameter. If not, function will return a null value.
#' @param s Random seed for the simulation.
#' @param toplot Show a 3D plot of the cluster points once generation is done?
#'   \code{TRUE} or \code{FALSE}.
#'
#' @return List of: [[1]] A \code{\link[spatstat]{pp3}} object containing the
#'   final locations of only type-A points within the final marked point
#'   pattern. [[2]] A \code{\link[spatstat]{pp3}} object containing the full
#'   marked underlying point pattern with points marked as either type A (cluser
#'   type point in a cluster), B (cluster type point not in a cluster) or C
#'   (non-cluster type point). [[3]] A vector containing the simulated radius of
#'   each cluster in the final point pattern. [[4]] A single value containing
#'   the true fraction of type-A points in the final simulated pattern.
#' @export
clustersim <- function(under, over, rcp_rad,
                       pcp = 0.1,
                       cr,
                       rho1, rho2,
                       rb = 0.1, # Radius blur (AS A PERCENT OF CR)
                       pb = 0.1, # Position blur (AS A PERCENT OF AVERAGE CLUSTER SEPARATION)
                       tol = 0.005,
                       s = 103,
                       toplot = F){
  set.seed(s)
  # Total volume
  sidelength <- diff(domain(under)$xrange)
  vt <- sidelength^3

  # Calculate volume needed in clusters
  if(rho2 >= pcp){
    print('Background density (rho2) cannot be larger than or equal to',
    'cluster point concentration (pcp).')
    return(-1)
  }
  vc <- vt*(pcp - rho2)/(rho1-rho2)

  # Calculate sigma (radius blur sd)
  sigma <- cr*rb

  # Calculate guess RCP concentration
  alpha <- 1
  rcp.conc <- vc/(4/3 * pi * cr*(cr^2 + 3*sigma^2) * vt * alpha)
  # Scale over RCP pattern to match this concentration
  over.vol <- volume(domain(over)) # Original volume
  over.vol.new <- npoints(over)/rcp.conc # New volume
  over.sf <- over.vol.new^(1/3)/over.vol^(1/3) # Scale factor to get to new volume
  over.scaled <- pp3_scale(over, over.sf)

  # Generate and add position blur
  over.sub <- over_cut(over.scaled, rep(sidelength, 3), cr)
  nnd <- nndist.pp3(over.sub)
  avg.sep <- mean(nnd)
  pb.sig <- avg.sep*pb

  pb.shifts <- pb_gen(npoints(over.scaled), mean = 0, sd = pb.sig)
  over.coo <- coords(over.scaled)
  over.coo.pb <- over.coo + pb.shifts
  over.scaled.pb <- pp3(x = over.coo.pb$x,
                        y = over.coo.pb$y,
                        z = over.coo.pb$z,
                        domain(over.scaled))

  # Re-size and shift the over point pattern so that it lines up with
  # the under point pattern with buffers
  over.scaledcut <- over_cut(over.scaled.pb, rep(sidelength, 3), 2.5*cr)
  over.scaledcut.coo <- coords(over.scaledcut)

  # Generate normal distributed radii for cluster centers
  cr.rand <- rnorm(npoints(over.scaledcut), mean = cr, sd = sigma)
  cr.rand[cr.rand < 0] <- 0
  marks(over.scaledcut) <- cr.rand

  # Check volume in clusters - optimize to target volume
  sf.optim <- optim(par = 1, fn = check_vol, gr = NULL,
                    over.scaledcut.coo, vc, domain(under), cr.rand,
                    method = "L-BFGS-B", lower = 0.5, upper = 3)
  over.scaled2 <- pp3_scale(over.scaledcut, sf.optim$par)

  # shift points to remove overlaps
  over.nolap <- overlap_fix(over.scaled2, cr.rand)
  if(is.numeric(over.nolap)){return(-1)}
  over.nolap.coo <- coords(over.nolap)

  # Re-check once for volume correctness
  sf.optim2 <- optim(par = 1, fn = check_vol, gr = NULL,
                     over.nolap.coo, vc, domain(under), cr.rand,
                     method = "L-BFGS-B", lower = 0.5, upper = 3)
  over.scaled3 <- pp3_scale(over.nolap, sf.optim2$par)
  over.final <- over_cut(over.scaled3, rep(sidelength, 3), cr.rand)

  # Select points that fall within the correct radius of each cluster center
  cr.rand.final <- marks(over.final)
  cluster.inds.all <- list()

  if(!(is.numeric(max(cr.rand.final)) && length(max(cr.rand.final)) == 1L &&
       max(cr.rand.final) >= 0)){
    print('Error with cr.rand.final')
    return(-1)
  }

  nnR <- crosspairs.pp3(over.final, under, rmax = max(cr.rand.final),
                        what = 'ijd', neat = TRUE, distinct = TRUE,
                        twice = FALSE)
  if(is.empty(nnR$i)){
    print('No cluster centers in domain.')
    return(-1)
  }
  nnR.split <- list()
  nnR.split$d <- split(nnR$d, nnR$i, drop=FALSE)
  nnR.split$j <- split(nnR$j, nnR$i, drop=FALSE)
  nnR.split$i <- as.numeric(attr(nnR.split$d, 'name'))

  cluster.inds.all <- lapply(1:length(nnR.split$i), function(k){
    nnR.split$j[[k]][nnR.split$d[[k]] < cr.rand.final[nnR.split$i[[k]]]]
  })

  cluster.inds.thinned <- lapply(cluster.inds.all, function(x){
    sample(x, round(rho1*length(x)))
  })

  cluster.inds <- unlist(cluster.inds.thinned)
  attr(cluster.inds, 'names') <- NULL

  # Select background points
  inds.all <- 1:npoints(under)
  bgnd.inds.all <- inds.all[!(inds.all %in% unlist(cluster.inds.all))]
  bgnd.inds <-  sample(bgnd.inds.all, round(rho2*length(bgnd.inds.all)))

  # Create final pp3
  c.marks <- rep(NA, npoints(under))
  c.marks[cluster.inds] <- 'A'
  c.marks[bgnd.inds] <- 'B'
  c.marks[-(c(cluster.inds, bgnd.inds))] <- 'C'

  marks(under) <- c.marks
  just.cluster.points <- under[marks(under) == 'A' | marks(under) == 'B']

  if(toplot == TRUE){
    plot3d.pp3(just.cluster.points)
  }

  pcp.real <- (length(bgnd.inds) + length(cluster.inds))/npoints(under)
  if(pcp.real < (pcp - tol) | pcp.real > (pcp + tol)){
    print('Pcp outside of tolerance range')
    return(-1)
  }
  #print(pcp.real)
  dat <- list(clusters = just.cluster.points,
              all = under,
              radii = cr.rand.final,
              perc = pcp.real,
              centers = over.final)
  return(dat)
}


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
#' Hello friends. A quick update as of 11/8/18; I have added quite a few
#' features and upgrades to this function, but mostly just to the cr superfast
#' type. It should be fairly easy to generalize these upgrades to the other
#' types, but I haven't gotten around to it yet. Some time soon perhaps. Let me,
#' Galen Vincent, know if you need some upgrades sooner. Thanks!
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
#' Currently written for application to \code{type = "cr"} \code{speed = "fast"}
#' and \code{"superfast"}, and to \code{type = "dist"}.
#' @param gb \code{TRUE} or \code{FALSE}. Whether or not to apply a gaussian
#'   blur to the cluster center positions. See \code{\link{rgblur}} for more
#' information on the blur function.
#' @param gbp Parameters for the cluster center gaussian blurs. \code{gbp =
#'   c(mean, sd)}. Default is mean = 0, sd = 1.
#' @param gbmethod See \code{method} argument in \code{\link{rgblur}}.
#' @param rb \code{TRUE} or \code{FALSE}. Whether or not to apply a gaussian
#'   blur to the cluster radius.
#' @param rbp Parameters for the cluster radius blurs. If \code{rbmethod = 1},
#'   then \code{rbp = sd} for the normal distribution. If \code{rbmethod = 2},
#'   then \code{rbp = c(p, r1, r2))}, where p is the percent of r1, and (1-p) is
#'   then the percent of r2. r1 < r2.
#' @param rbmethod Metod for distributing radius blur. \code{rbmethod = 1} means
#'   gaussian distributed, \code{rbmethod = 2} mean split between two radii. Can
#'   ignore \code{cr} if \code{rbmethod = 2}.
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
#'   \subsection{\code{type} = "cr", speed = "superfast", rb = FALSE}{List of: [[1]]
#'   \code{\link[spatstat]{pp3}} object containing the cluster marked point
#'   locations. [[2]] \code{\link[spatstat]{pp3}} object containing the
#'   overlaying RCP pattern after scaling. [[3]] number of points put in or
#'   taken away at random. [[4]] Vector of nearest neighbor distance (cluster
#'   center to center) from each cluster}
#'   \subsection{\code{type} = "cr", speed = "superfast", rb = TRUE}{List of:
#'   [[1]] \code{\link[spatstat]{pp3}} object containing the cluster marked
#'   point locations. [[2]] \code{\link[spatstat]{pp3}} object containing the
#'   overlaying RCP pattern after scaling. [[3]] number of points put in or
#'   taken away at random. [[4]] Vector of nearest neighbor distance (cluster
#'   center to center) from each cluster. [[5]] A vector contining the radius of
#'   each cluster.}

makecluster <- function(under,over,radius1,radius2,
                        type = "cr",
                        ppc=NULL,
                        cr=NULL,speed = "superfast",
                        d=NULL,
                        pic = 1,
                        pcp = 0.06,
                        den = 1,
                        gb = FALSE,
                        gbp = c(0,1),
                        gbmethod = 1,
                        rb = FALSE,
                        rbp = 1,
                        rbmethod = 1,
                        s = 100,
                        toPlot=FALSE,showOverPts=FALSE){
  #############################################################################
  # POINTS PER CLUSTER METHOD
  if(type == "ppc"){
    #real cluster percent
    if(gb == TRUE){
      print("ERROR - Gaussian blur only implemented for type = cr, speed = superfast.")
      return()
    }
    rcp <- pcp*pic

    under.r <- radius1
    over.r <- radius2
    over.rf <- under.r*(ppc/rcp)^(1/3)
    over.scaled <- scaleRCP(over, newRadius = over.rf, oldRadius = over.r,
                            win = domain(over))
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
      cluster.ind1 <- nncross(over.scaledf, under, what = "which", k = 1:ppc)
      cluster.ind <- NULL
      for(i in 1:npoints(over.scaledf)){
        cluster.ind <- c(cluster.ind,as.numeric(cluster.ind1[i,]))
      }
    }

    more <- npoints(under)*pcp-npoints(under)*rcp
    if(more==0){

    }else{
      cluster.ind <- randomInsert(cluster.ind, more, npoints(under), s)
    }

    cluster.xyz <- coords(under)[cluster.ind,]
    cluster <- createSpat(cluster.xyz)

    if(toPlot==TRUE){
      plot3d.pp3(cluster,col="red",size=5)
      #plot3d.pp3(under,col="lightgray",add=TRUE)
      if(showOverPts==TRUE){
        plot3d.pp3(over.scaledf,size= 6,col="black",add=TRUE)
      }
    }

    return(list(cluster,over.scaledf,c(ppc,npoints(over.scaledf)-diff,ppc+1,diff)))

  }

  ########################################################################################
  # CHOOSE RADIUS METHOD
  else if(type == "cr"){
    if (speed == "slow"){
      if(gb == TRUE){
        print("ERROR - Gaussian blur only implemented for type = cr, speed = superfast.")
        return()
      }
      #real cluster percent
      rcp <- pcp*pic

      under.r <- radius1
      over.r <- radius2
      under.vol <- volume(domain(under))

      over.rf <- under.r*cr*((4*pi*npoints(under))/(3*under.vol*rcp))^(1/3)

      over.scaled <- scaleRCP(over, newRadius = over.rf, oldRadius = over.r,
                              win = domain(over))
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
        plot3d.pp3(cluster,col="black",size=5)
        #plot3d.pp3(under,col="lightgray",add=TRUE)
        if(showOverPts==TRUE){
          plot3d.pp3(over.scaledf,size= 6,col="red",add=TRUE)
        }
      }

      return(list(cluster,over.scaledf,c(npoints(over.scaledf),npoints(cluster)-more)))
    }
    else if (speed == "fast"){
      if(gb == TRUE){
        print("ERROR - Gaussian blur only implemented for type = cr, speed = superfast.")
        return()
      }
      #real cluster percent
      rcp <- pcp*pic

      under.r <- radius1
      over.r <- radius2
      under.vol <- volume(domain(under))

      over.rf <- under.r*cr*((4*pi*npoints(under))/(3*under.vol*rcp))^(1/3)

      over.scaled <- scaleRCP(over, newRadius = over.rf, oldRadius = over.r,
                              win = domain(over))
      over.scaledf <- subSample(under,over.scaled)

      cluster.nnR.new <- crosspairs.pp3(over.scaledf, under, cr, what="indices",
                                        twice=FALSE, distinct=TRUE, neat=TRUE)

      cluster.ind <- cluster.nnR.new[[2]]
      cluster.info <- factor(cluster.nnR.new[[1]])
      diff <- round(rcp*npoints(under)-length(cluster.ind))

      cluster.adj <- crAdjust.new(cluster.ind, cluster.info, diff, over.scaledf,
                                  under)
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
        cluster.ind <- randomInsert(cluster.ind, more, npoints(under), s,
                                    cluster.ind.split)
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
    else if (speed == "superfast"){
      #browser()
      # If rb is true, we need to set a new cr to deal with spacing - Added 2/4/19
      if(rb == TRUE){
        if(rbmethod == 1){
          #scale cr to average *volume*
          cr <- (cr * (cr^2 + 3*rbp^2))^(1/3)
          #print(cr)
        }else if(rbmethod == 2){
          #scale cr to average *volume*
          cr <- (rbp[1]*rbp[2]^3 + (1-rbp[1])*rbp[3]^3)^(1/3)
        }
      }

      #real cluster percent
      set.seed(s)
      rcp <- pcp*pic

      under.r <- radius1
      over.r <- radius2
      over.vol <- volume(domain(over))
      under.vol <- volume(domain(under))

      # Calculate scaling factor for over point pattern
      # Added 2/4/19
      ####
      # There's some hairy math here that took me a while to figure out. Come see me if you need it explained,
      # or see my reseatch notebook.
      if(rb == TRUE){
        z63 <- (((4/3)*pi*over.r^3)*npoints(over)*0.75 +
                  ((4/3)*pi*(over.r*1.20)^3)*npoints(over)*0.25) / over.vol
        under.xdim <- domain(under)$xrange[2]
        under.ydim <- domain(under)$yrange[2]
        under.zdim <- domain(under)$zrange[2]

        innervol <- (under.xdim - 2 * cr)*(under.ydim - 2 * cr)*(under.zdim - 2 * cr)
        outervol <- (under.xdim + 2 * cr)*(under.ydim + 2 * cr)*(under.zdim + 2 * cr)
        middlevol <- outervol - innervol

        volfactor <- ((innervol/outervol) + (middlevol/outervol)*(5/12))

        over.rf <- cr * ((under.xdim + 2*cr) * (under.ydim + 2*cr) *
                           (under.zdim + 2*cr) * z63 * volfactor /
                           (under.vol * rcp * 1.182))^(1/3)

        if(rbmethod == 1){
          sdfactor <- (rbp/cr)*0.37
          over.rf <- over.rf + over.rf*sdfactor
        }
      } else {
        over.rf <- under.r*cr*((4*pi*npoints(under))/(3*under.vol*rcp))^(1/3)
      }
      #####
      over.sep <- over.rf*2
      over.scaled <- scaleRCP(over,newRadius = over.rf, oldRadius = over.r,win = domain(over))
      #browser()
      if(gb == TRUE){
        n <- npoints(over.scaled)
        gbval <- rgblur(n,gbp[1],gbp[2],coords = "rec", method = gbmethod)
        over.xyz <- coords(over.scaled)
        over.xyz.new <- over.xyz + gbval
        over.scaled.new <- createSpat(over.xyz.new, win = domain(over.scaled))
      }else{
        over.scaled.new <- over.scaled
      }

      # new addition as of 11/8/2018 - use cluster centers that can fall outside of the under pattern domain. This reduces
      # number of cluster points error significantly
      under.coo <- coords(under)
      under.new <- createSpat(under.coo + cr)
      over.scaled.domain <- c(domain(under)$xrange[2]+2*cr,domain(under)$yrange[2]+2*cr,domain(under)$zrange[2]+2*cr)
      over.scaledf <- subSquare(over.scaled.new, over.scaled.domain)

      if(rb == TRUE){
        n <- npoints(over.scaledf)
        if(rbmethod == 1){
          crrand <- rnorm(n,mean = cr, sd = rbp)
          crrand[crrand < 0] <- 0
        } else if(rbmethod == 2){
          n1 <- round(n*rbp[1])
          n2 <- n - n1
          a <- c(rep(rbp[2], n1), rep(rbp[3], n2))
          crrand <- sample(a, replace = FALSE)
        }
      }

      # Deal with overlapping clusters here
      # Check for overlap, seperate to no overlap if so
      if(gb == TRUE){
        nnd <- nndist.pp3(over.scaledf)
        if(rb == TRUE){
          scrrand <- sort(crrand,decreasing = TRUE)
          comp <- scrrand[1] + scrrand[2]
        }else{
          comp <- 2*cr
        }

        if(any(nnd < comp)){
          #browser()
          check <- which(nnd < comp)
          lc <- 0
          t1 <- Sys.time()
          while(!is.empty(check)) {
            t2 <- Sys.time()
            nnw <- nnwhich.pp3(over.scaledf)
            direction <- (coords(over.scaledf)[nnw[check[1]],]-coords(over.scaledf)[check[1],])/nnd[check[1]]
            coords(over.scaledf)[check[1],] <- coords(over.scaledf)[check[1],]+(nnd[check[1]]-(comp+0.00001))*direction
            nnd <- nndist.pp3(over.scaledf)
            check <- which(nnd < comp)
            if((as.numeric(t2) - as.numeric(t1)) > 15){
              return(-1)
            }

          }
          #browser()
        }
      }

      #rb
      if(rb == TRUE){
        cluster.nnR.ind1 <- list()
        cluster.nnR.ind2 <- list()
        #t1 <- Sys.time()
        cluster.nnR.full <- crosspairs.pp3(over.scaledf, under.new, max(crrand),
                                           what = "all", twice = FALSE,
                                           distinct = TRUE, neat = TRUE)
        cluster.nnR.split <- list()
        cluster.nnR.split$d <- split(cluster.nnR.full$d, cluster.nnR.full$i,
                                     drop = FALSE)
        cluster.nnR.split$j <- split(cluster.nnR.full$j, cluster.nnR.full$i,
                                     drop = FALSE)
        split.vals <- sort(unique(cluster.nnR.full$i))
        for(i in 1:length(split.vals)){
          cluster.nnR.ind2[[i]] <- cluster.nnR.split$j[[i]][cluster.nnR.split$d[[i]] < crrand[split.vals[i]]]
          cluster.nnR.ind1[[i]] <- rep(split.vals[i],length(cluster.nnR.ind2[[i]]))
        }
        #t2 <- Sys.time()
        #print(t2 - t1)
        cluster.nnR.new <- list()
        cluster.nnR.new[[1]] <- unlist(cluster.nnR.ind1)
        cluster.nnR.new[[2]] <- unlist(cluster.nnR.ind2)
      }else{
        cluster.nnR.new <- crosspairs.pp3(over.scaledf, under.new, cr,
                                          what = "indices", twice = FALSE,
                                          distinct = TRUE, neat = TRUE)
      }


      if(den < 1 & den >= 0){
        cluster.ind.split <- split(cluster.nnR.new[[2]], cluster.nnR.new[[1]],
                                   drop = FALSE)
        cluster.ind.thinned <- lapply(cluster.ind.split,function(x){
          return(sample(x,round(den*length(x)),replace=FALSE))})
        cluster.ind <- unlist(cluster.ind.thinned)
        cluster.ind.split <- unlist(cluster.ind.split)
      }else {
        cluster.ind <- cluster.nnR.new[[2]]
        cluster.ind.split <- cluster.ind
      }

      more <- round(npoints(under.new)*pcp)-length(cluster.ind)
      if(more==0){

      }else if(more > 0){
        cluster.ind <- randomInsert(cluster.ind,more,npoints(under.new),s,cluster.ind.split)
      }else if(more < 0){
        cluster.ind <- randomTakeAway(cluster.ind,-1*more,npoints(under.new),s)
      }

      cluster.xyz <- coords(under)[cluster.ind,]
      cluster.xyz <- na.omit(cluster.xyz)
      cluster <- createSpat(cluster.xyz, win = domain(under))

      c.marks <- rep(NA, npoints(under))
      c.marks[cluster.ind] <- 'A'
      c.marks[-cluster.ind] <- 'B'

      marks(under) <- c.marks

      over.scaledf.coo <- coords(over.scaledf)
      over.scaledf <- createSpat(over.scaledf.coo - cr, win = domain(under))

      cluster.nn <- nndist.pp3(over.scaledf, k = 1)

      if(toPlot==TRUE){
        plot3d.pp3(cluster,col="red",size=5)
        #plot3d.pp3(under,col="lightgray",add=TRUE)
        if(showOverPts==TRUE){
          plot3d.pp3(over.scaledf,size= 6,col="black",add=TRUE)
        }
      }

      if(rb == TRUE){
        return(list(cluster,over.scaledf,more,cluster.nn,crrand,under))
      }else{
        return(list(cluster,over.scaledf,more,cluster.nn,under))
      }
    }
    else{
      print('Insert speed for cluster radius method. slow, fast, or superfast.')
    }
  }

  ###########################################################################################
  # CHOOSE DISTANCE BETWEEN CLUSTERS METHOD
  else if(type == "dist"){
    if(gb == TRUE){
      print("ERROR - Gaussian blur only implemented for type = cr, speed = superfast.")
      return()
    }
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
      cluster.ind <- randomTakeAway(cluster.ind,-1*more,npoints(under),s)
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

  }
  else{
    print("Please input a valid type")
    return()
  }
}
################################################################################

#### hcpcluster ####
#' Generate spherical clusters with hexagonal close packed spacing.
#'
#' Generates perfectly spherical clusters with structred order of hexagonal
#' close packing (HCP). Can define cluster separation, cluster radius, inter
#' cluster density, background density, and window size.
#'
#' @param csep_r Cluster separation radius. That is, the radius of the spheres
#'   in the HCP structure.
#' @param R Cluster radius.
#' @param sigma1 Inter cluster density. Value between 0 and 1.
#' @param sigma2 Background density. Value between 0 and 1.
#' @param win A \code{\link[spatstat]{box3}} object containing the window of the
#'   cluster set you want to make.
#' @param background Either \code{'poisson'} or \code{'rcp'}. Whether to have
#'   Poission distributed points or RCP points for the points in the clusters
#'   and in the background.
#' @param filepath Needed if \code{background = 'rcp'}. Vector with the filepath
#'   to [1] the FinalConfig file of the RCP pattern desired, [2] the system file
#'   of the RCP pattern desired
#'
#' @return A list of [[1]] A \code{\link[spatstat]{pp3}} object containing the
#'   cluster points, [[2]] The overall intensity of the point pattern; Total
#'   number of points/total volume.
#' @export

hcpcluster <- function(csep_r, R, sigma1, sigma2, win, background, filepath){
  hcp.c <- hcp.gen(csep_r,win)
  hcp.pp3 <- createSpat(hcp.c,win)
  tr <- inside.boxx(hcp.pp3,w=win)
  hcp.pp3 <- hcp.pp3[tr]

  if(background == 'poisson'){
    # Generate poission clusters around the hcp lattice
    inside.pp3 <- rpoispp3(sigma1, domain = win)
    outside.pp3 <- rpoispp3(sigma2, domain = win)
  }else if(background == 'rcp'){
    # OR upload rcp pattern from filepath given
    inside.pp31 <- read.rcp(filepath[1],filepath[2],scaleUp = TRUE,newRadius = 0.50475)
    inside.pp32 <- stitch.size(inside.pp31,domain(inside.pp31),c(60,60,60))
    inside.pp3 <- percentSelect(sigma1, inside.pp32)
    outside.pp3 <- percentSelect(sigma2, inside.pp32)
  }

  nnr <- crosspairs.pp3(hcp.pp3, inside.pp3, rmax = R, what ="indices")
  cluster.pts <- inside.pp3[nnr[[2]]]

  #Generate background points ouside the clusters
  nnr2 <- crosspairs.pp3(hcp.pp3, outside.pp3, rmax = R, what ="indices")
  bg.pts <- outside.pp3[-nnr2[[2]]]

  full.coo <- rbind(coords(cluster.pts),coords(bg.pts))
  full.pp3 <- createSpat(full.coo, win)

  sigma3 <- npoints(full.pp3)/volume(win)

  return(list(full.pp3,sigma3))
}

#### planes ####
#' Simulate lamellar morphology for mixed guest-host system.
#'
#' @param lambda Intensity for the background Poisson process.
#' @param frac Fraction of points to select as guest points. Between zero and
#'   one.
#' @param plane.norm Vector you wish to define normal of the planes with.
#'   Doesn't need to be normalized. Positive values only
#' @param plane.den Linear plane density.
#' @param point.den Point density within the plane volumes.
#' @param toplot \code{TRUE} or \code{FALSE}. Plot results or not.
#' @param win The simulation window.
#'
#' @return A list of: [[1]] a \code{\link[spatstat]{pp3}} object contining the
#'   guest points. [[2]] A \code{pp3} object containing the entire background
#'   underlying point pattern.
#' @export

morph_lamellar <- function(lambda,
                           frac,
                           plane.norm = c(1, 0, 0),
                           plane.den,
                           point.den = 1,
                           toplot = FALSE,
                           win = box3(c(0,1), c(0,1), c(0,1))){

  bgnd <- rpoispp3(lambda = lambda, domain = win)
  coo <- coords(bgnd)
  p.norm <- plane.norm
  p.den <- plane.den

  #plane density (planes per unit length)
  p.spacing <- 1/p.den

  #plane orientation (vector holding normal direction for planes, doesn't have to be normalized)
  p.direc <- p.norm/sqrt(sum(p.norm^2)) * p.spacing

  #points on planes:
  p.points <- list()
  p.points[[1]] <- c(-1, -1, -1)
  nxt <- p.points[[1]] + p.direc
  i <- 1
  while(nxt[1] <= win$xrange[2]*1.5 &
        nxt[2] <= win$yrange[2]*1.5 &
        nxt[3] <= win$zrange[2]*1.5){
    p.points[[i+1]] <- nxt
    nxt <- p.points[[i+1]] + p.direc
    i <- i + 1
  }

  offset <- runif(1, 0, 1)
  p.points <- lapply(p.points, function(x){x + p.norm/sqrt(sum(p.norm^2)) * offset
    print(x + p.norm/sqrt(sum(p.norm^2)) * offset)})

  #select points within x distance of planes to be type-A
  distmat <- matrix(FALSE, npoints(bgnd), length(p.points))
  p.norm.normed <- p.norm/sqrt(sum(p.norm^2))

  for(i in 1:length(p.points)){
    vecs <- t(apply(coo, 1, function(x){x - p.points[[i]]}))
    distmat[,i] <- abs(vecs %*% p.norm.normed)
  }

  nkeep <- round(frac*npoints(bgnd))
  mindistlist <- apply(distmat, 1, min)
  gbwhich <- order(mindistlist)[1:nkeep]

  if(point.den != 1){
    nin <- round(point.den*length(gbwhich))
    nout <- length(gbwhich) - nin
    points.final <- sample(gbwhich, nin)

    nonplane.ind <- (1:npoints(bgnd))[-gbwhich]
    noise <- sample(nonplane.ind, nout)
    total.inds <- c(points.final, noise)
  }else{
    total.inds <- gbwhich
  }

  p.selected <- bgnd[total.inds]

  if(toplot == TRUE){
    plot3d.pp3(p.selected, col = 'red', xlim = bgnd$domain$xrange, ylim = bgnd$domain$yrange, zlim = bgnd$domain$zrange)
  }

  return(p.selected)
}

#### rods ####
#' Simulate morphology of guest rods in a two-phase system
#'
#' @param lambda Intensity for the background Poisson process.
#' @param frac Fraction of points to select as guest points. Between zero and
#'   one.
#' @param rod.norm Either 'x', 'y', or 'z'; defines which axis the rods should
#'   be alligned with.
#' @param rod.den Rod density; 2D vector containing a linear density for both
#'   directions.
#' @param rod.spacing One of: "grid", "rcp", "hexagonal", or "random". The
#'   method used to space the rods.
#' @param rcp.path If \code{rod.spacing} is set to 'rcp', then set the file path
#'   to the folder contining the rcp FinalConfig and system files
#' @param point.den Point density within the volume of the rods.
#' @param toplot \code{TRUE} or \code{FALSE}; Plot results or not.
#' @param win The simulation window.
#'
#' @return A list of: [[1]] a \code{\link[spatstat]{pp3}} object contining the
#'   guest points. [[2]] A \code{pp3} object containing the entire background
#'   underlying point pattern.
#' @export

morph_rods <- function(lambda,
                       frac,
                       rod.norm = 'z',
                       rod.den,
                       rod.spacing = 'grid',
                       rcp.path = 'C:/Users/galen/Documents/Research/point_patterns/2D/Final',
                       point.den = 1,
                       toplot = FALSE,
                       win = box3(c(0,1), c(0,1), c(0,1))){
  # generate points
  bgnd <- rpoispp3(lambda = lambda, domain = win)
  coo <- coords(bgnd)

  # grid of rods using the appropriate method and axis
  if(rod.norm == 'z') {
    comp.points <- data.frame(coo$x, coo$y)
    xdist <- diff(win$xrange)
    ydist <- diff(win$yrange)
  }
  else if (rod.norm == 'y') {
    comp.points <- data.frame(coo$x, coo$z)
    xdist <- diff(win$xrange)
    ydist <- diff(win$zrange)
  }
  else if (rod.norm == 'x') {
    comp.points <- data.frame(coo$y, coo$z)
    xdist <- diff(win$yrange)
    ydist <- diff(win$zrange)
  }
  else {
    print('Please input one of \'x\', \'y\', or \'z\' for rod.norm')
  }

  # create rod spacing
  buffer <- max(comp.points)*0.1

  if(rod.spacing == 'grid'){
    rod.x <- seq(min(comp.points[,1]) - 1 - buffer, max(comp.points[,1]) + 1 + buffer, by = 1/rod.den[1])
    rod.y <- seq(min(comp.points[,2]) - 1 - buffer, max(comp.points[,2]) + 1 + buffer, by = 1/rod.den[2])
    offset <- c(runif(1,-1, 1), runif(1,-1,1))
    rod.x <- rod.x + offset[1]
    rod.y <- rod.y + offset[2]
    rod.x <- rod.x[rod.x > (min(comp.points[,1]) - buffer) & rod.x < (max(comp.points[,1]) + buffer)]
    rod.y <- rod.y[rod.y > (min(comp.points[,2]) - buffer) & rod.y < (max(comp.points[,2]) + buffer)]
    rod.xy <- expand.grid(rod.x, rod.y)
  }
  else if (rod.spacing == 'rcp') {
    rcp.upload <- fread(paste(rcp.path, '/FinalConfig', sep = ''))
    rcp <- ppp(rcp.upload$V1, rcp.upload$V2)

    rcp.den <- rcp$n/area(rcp$window)
    scaling.factor <- sqrt(rcp.den/(rod.den[1]*rod.den[2]))

    rcp.upload.scaled <- coords(rcp)*scaling.factor
    rcp.scaled <- ppp(rcp.upload.scaled$x, rcp.upload.scaled$y,
                      window = owin(c(0,max(rcp.upload.scaled$x)),
                                    c(0,max(rcp.upload.scaled$y))))

    xmu <- mean(rcp.scaled$window$xrange)
    ymu <- mean(rcp.scaled$window$xrange)

    win.select <- owin(c(xmu - xdist/2 - buffer, xmu + xdist/2 + buffer),
                       c(ymu - ydist/2 - buffer, ymu + ydist/2 + buffer))
    rod.xy <- coords(rcp.scaled[inside.owin(rcp.scaled$x, rcp.scaled$y, win.select)])
    rod.xy$x <- rod.xy$x - xmu + xdist/2
    rod.xy$y <- rod.xy$y - ymu + ydist/2
  }
  else if (rod.spacing == 'hexagonal') {
    hex <- hexgrid(owin(c(0,1), c(0,1)), 0.01)

    hex.den <- hex$n/area(hex$window)
    scaling.factor <- sqrt(hex.den/(rod.den[1]*rod.den[2]))

    hex.scaled.coo <- coords(hex)*scaling.factor
    offset <- c(runif(1, 0, 1), runif(1, 0, 1))

    hex.scaled <- ppp(hex.scaled.coo$x + offset[1], hex.scaled.coo$y + offset[2],
                      window = owin(c(0,max(hex.scaled.coo$x + offset[1])),
                                    c(0,max(hex.scaled.coo$y + offset[2]))))

    xmu <- mean(hex.scaled$window$xrange)
    ymu <- mean(hex.scaled$window$xrange)

    win.select <- owin(c(xmu - xdist/2 - buffer, xmu + xdist/2 + buffer),
                       c(ymu - ydist/2 - buffer, ymu + ydist/2 + buffer))
    rod.xy <- coords(hex.scaled[inside.owin(hex.scaled$x, hex.scaled$y, win.select)])
    rod.xy$x <- rod.xy$x - xmu + xdist/2
    rod.xy$y <- rod.xy$y - ymu + ydist/2
  }
  else if (rod.spacing == 'random') {
    xrng <- win$xrange*10
    yrng <- win$yrange*10
    mat <- rMaternII((rod.den[1]*rod.den[2])*1.25, buffer*2, owin(xrng, yrng))

    xmu <- mean(mat$window$xrange)
    ymu <- mean(mat$window$xrange)

    win.select <- owin(c(xmu - xdist/2 - buffer, xmu + xdist/2 + buffer),
                       c(ymu - ydist/2 - buffer, ymu + ydist/2 + buffer))
    rod.xy <- coords(mat[inside.owin(mat$x, mat$y, win.select)])
    rod.xy$x <- rod.xy$x - xmu + xdist/2
    rod.xy$y <- rod.xy$y - ymu + ydist/2
  }

  distmat <- matrix(FALSE, npoints(bgnd), nrow(rod.xy))
  for(i in 1:nrow(rod.xy)){
    distmat[,i] <- sqrt((comp.points[,1] - rod.xy[i,1])^2 + (comp.points[,2] - rod.xy[i,2])^2)
  }

  nkeep <- round(frac*npoints(bgnd))
  mindistlist <- apply(distmat, 1, min)
  gbwhich <- order(mindistlist)[1:nkeep]

  if(point.den != 1){
    nin <- round(point.den*length(gbwhich))
    nout <- length(gbwhich) - nin
    points.final <- sample(gbwhich, nin)

    nonplane.ind <- (1:npoints(bgnd))[-gbwhich]
    noise <- sample(nonplane.ind, nout)
    total.inds <- c(points.final, noise)
  }else{
    total.inds <- gbwhich
  }

  p.selected <- bgnd[total.inds]

  if(toplot == TRUE){
    plot3d.pp3(p.selected, col = 'red')
  }

  return(p.selected)
}

#### Gyroid ####
#' Simulate morphology of guest molecules in gyroid pattern
#'
#' @param lambda Intensity for the background Poisson process.
#' @param frac Fraction of points to select as guest points. Between zero and
#'   one.
#' @param gyroid.scale Scale factor for the gyroid pattern. Larger = tighter
#'   spacing of the gyroid.
#' @param point.den Point density within the volume of the guest volume.
#' @param toplot \code{TRUE} or \code{FALSE}; Plot results or not.
#' @param win The simulation window.
#'
#' @return A list of: [[1]] a \code{\link[spatstat]{pp3}} object contining the
#'   guest points. [[2]] A \code{pp3} object containing the entire background
#'   underlying point pattern.
#' @export

morph_gyroid <- function(lambda,
                         frac,
                         gyroid.scale,
                         point.den = 1,
                         toplot = FALSE,
                         win = box3(c(0,1), c(0,1), c(0,1))){

  shift <- c(runif(1, 0, 1000), runif(1, 0, 1000), runif(1, 0, 1000))

  gyr <- function(x, y, z, sc){
    return(sin(8*sc*x + shift[1])*cos(8*sc*y + shift[2]) +
             sin(8*sc*y + shift[2])*cos(8*sc*z + shift[3]) +
             sin(8*sc*z + shift[3])*cos(8*sc*x + shift[1]))
  }

  bgnd <- rpoispp3(lambda = lambda, domain = win)
  coo <- coords(bgnd)

  coo.gyr <- gyr(coo$x, coo$y, coo$z, gyroid.scale)
  nkeep <- round(frac*npoints(bgnd))

  gbwhich <- tail(order(coo.gyr), n = nkeep)

  if(point.den != 1){
    nin <- round(point.den*length(gbwhich))
    nout <- length(gbwhich) - nin
    points.final <- sample(gbwhich, nin)

    nonplane.ind <- (1:npoints(bgnd))[-gbwhich]
    noise <- sample(nonplane.ind, nout)
    total.inds <- c(points.final, noise)
  }else{
    total.inds <- gbwhich
  }

  p.selected <- bgnd[total.inds]

  if(toplot == TRUE){
    plot3d.pp3(p.selected, col = 'red')
  }

  return(p.selected)
}

#### Grain Boundary ####
#' Simulate morphology of guest molecules on a grain boundary
#'
#' @param lambda Intensity for the background Poisson process.
#' @param frac Fraction of points to select as guest points. Between zero and
#'   one.
#' @param point.den Point density within the volume of the guest volume.
#' @param rcp.rad Radius to scale the RCP pattern to.
#' @param rcp.path File path to folder containing the 3D RCP patterns which
#'   contain both the 'FinalConfigx' and 'systemx' files.
#' @param rcp.number Which RCP file to use.
#' @param toplot \code{TRUE} or \code{FALSE}; Plot results or not.
#' @param win The simulation window.
#'
#' @return A list of: [[1]] a \code{\link[spatstat]{pp3}} object contining the
#'   guest points. [[2]] A \code{pp3} object containing the entire background
#'   underlying point pattern.
#' @export

morph_gb <- function(lambda,
                     frac,
                     point.den = 1,
                     rcp.rad,
                     rcp.path = NULL,
                     rcp.number = 1,
                     toplot = FALSE,
                     win = box3(c(0,1), c(0,1), c(0,1))){

  if(rcp.number == 'rand'){
    rcp.number <- sample(1, 1:523)
  }
  bgnd <- rpoispp3(lambda, domain = win)
  rcp <- read.rcp(paste(rcp.path, '/FinalConfig', toString(rcp.number), sep = ''),
                  paste(rcp.path, '/system', toString(rcp.number), sep = ''),
                  scaleUp = TRUE, newRadius = rcp.rad)
  xmu <- mean(rcp$domain$xrange)
  ymu <- mean(rcp$domain$yrange)
  zmu <- mean(rcp$domain$zrange)

  xdist <- diff(win$xrange)
  ydist <- diff(win$yrange)
  zdist <- diff(win$zrange)
  offset <- c(runif(1, -1, 1), runif(1, -1, 1), runif(1, -1, 1))

  win.select <- box3(c(xmu - xdist/2 - rcp.rad - 1 + offset[1],
                       xmu + xdist/2 + rcp.rad + 1 + offset[1]),
                     c(ymu - ydist/2 - rcp.rad - 1 + offset[2],
                       ymu + ydist/2 + rcp.rad + 1 + offset[2]),
                     c(zmu - zdist/2 - rcp.rad - 1 + offset[3],
                       zmu + zdist/2 + rcp.rad + 1 + offset[3]))
  rcp.xyz <- coords(rcp[inside.boxx(rcp, w = win.select)])
  rcp.xyz$x <- rcp.xyz$x - xmu + xdist/2 - offset[1]
  rcp.xyz$y <- rcp.xyz$y - ymu + ydist/2 - offset[2]
  rcp.xyz$z <- rcp.xyz$z - zmu + zdist/2 - offset[3]
  rcp.cut <- pp3(rcp.xyz$x, rcp.xyz$y, rcp.xyz$z,
                 c(min(rcp.xyz$x), max(rcp.xyz$x)),
                 c(min(rcp.xyz$y), max(rcp.xyz$y)),
                 c(min(rcp.xyz$z), max(rcp.xyz$z)))

  nnc <- nncross(bgnd, rcp.cut, k = 1:2)

  diff <- abs(nnc$dist.1 - nnc$dist.2)
  nkeep <- round(frac*npoints(bgnd))
  gbwhich <- order(diff)[1:nkeep]

  if(point.den != 1){
    nin <- round(point.den*length(gbwhich))
    nout <- length(gbwhich) - nin
    points.final <- sample(gbwhich, nin)

    nonplane.ind <- (1:npoints(bgnd))[-gbwhich]
    noise <- sample(nonplane.ind, nout)
    total.inds <- c(points.final, noise)
  }else{
    total.inds <- gbwhich
  }

  p.selected <- bgnd[total.inds]

  if(toplot == TRUE){
    plot3d.pp3(p.selected, col = 'red')
  }

  return(p.selected)
}

#########################################
# Helper functions
#### Pre 1/21/2020 helpers: ####
#### bcc.gen ####
#' Generate body centerd cubic (BCC) lattice of points in 3D.
#'
#' Generates a BCC point pattern with a specified window size and number of
#' points desired in the pattern.
#'
#' @param npoint Approximate number of points to have in the pattern.
#' @param win A \code{\link[spatstat]{box3}} object containing the window for
#'   the final pattern.
#'
#' @return A \code{\link[spatstat]{pp3}} object with the lattice points.
#'
#' @name bcc.gen-deprecated
#' @seealso \code{\link{rapt-deprecated}}
#' @keywords internal
NULL
#' @rdname rapt-deprecated
#' @section \code{bcc.gen}:
#'   For \code{bcc.gen}, use \code{\link{lattice}}
#'
#' @export
bcc.gen <- function(npoint, win){
  vol <- volume(win)
  n.units <- npoint/2
  vol.per.unit <- vol/n.units
  a <- (vol.per.unit)^(1/3)

  x <- seq(win$xrange[1], win$xrange[2], by = a)
  y <- seq(win$yrange[1], win$yrange[2], by = a)
  z <- seq(win$zrange[1], win$zrange[2], by = a)

  grid <- expand.grid(x, y, z)

  x.c <- seq(win$xrange[1] + a/2, win$xrange[2] - a/2, by = a)
  y.c <- seq(win$yrange[1] + a/2, win$yrange[2] - a/2, by = a)
  z.c <- seq(win$zrange[1] + a/2, win$zrange[2] - a/2, by = a)

  grid.c <- expand.grid(x.c, y.c, z.c)

  grid.full <- rbind(grid, grid.c)
  names(grid.full) <- c('x','y','z')

  lat.pp3 <- pp3(grid.full$x, grid.full$y, grid.full$z, win)

  return(lat.pp3)
}

#### hcp.gen ####
#' Helper for \code{\link{hcpcluster}} which generates a HCP lattice with the correct spacing.
#'
#' Generates a hexagonal close packed (HCP) structure based on spheres with
#' specified radius within a window of specified size.
#'
#' @param r The radius of the spheres for the HCP structure.
#' @param win A \code{\link[spatstat]{box3}} object defining the window for the
#'   generation.
#'
#' @return A data frame object containing xyz coordinates of the HCP lattice
#'   centers.
#' @seealso \code{\link{hcpcluster}}

hcp.gen <- function(r, win){
  xyz <- list()
  i <- 0
  j <- 0
  k <- 0
  count <- 1
  xyz[[count]] <- r*c(2*i + ((j+k)%%2), sqrt(3)*(j + (1/3)*(k%%2)), (2*sqrt(6)/3)*k)
  count <- count + 1
  i <- -1
  while(xyz[[count-1]][1] < win$xrange[2]){
    j <- 0
    xyz[[count]] <- r*c(2*i + ((j+k)%%2), sqrt(3)*(j + (1/3)*(k%%2)), (2*sqrt(6)/3)*k)
    count <- count + 1
    while(xyz[[count-1]][2] < win$yrange[2]){
      k <- 0
      xyz[[count]] <- r*c(2*i + ((j+k)%%2), sqrt(3)*(j + (1/3)*(k%%2)), (2*sqrt(6)/3)*k)
      count <- count + 1
      while(xyz[[count-1]][3] < win$zrange[2]){
        k <- k + 1
        xyz[[count]] <- r*c(2*i + ((j+k)%%2), sqrt(3)*(j + (1/3)*(k%%2)), (2*sqrt(6)/3)*k)
        count <- count + 1
      }
      j <- j +1
    }
    i <- i+1
  }
  xyzdf <- do.call(rbind.data.frame, xyz)
  colnames(xyzdf) <- c("x","y","z")
  return(xyzdf)
}

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

randomInsert <- function(cluster.indices, n, N,s, points.avoid = cluster.indices){
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

randomTakeAway <- function(cluster.indices, n, N, s){
  #cluster.Indices is a vector containig the indices of the current cluster points
  #n is the number of points that need to be placed randomly
  #N is the number of points in the entire underlying pattern

  set.seed(s)
  inds <- sample(1:length(cluster.indices),n)
  cluster.indices[inds] <- NaN

  all.ind <- cluster.indices[!is.nan(cluster.indices)]

  return(all.ind)
}

#### rgblur ####
#' Helper for \code{\link{makecluster}} to 3D gaussian blur cluster positions.
#'
#' Function which returns a list of coordinates of radial gaussian blurred
#' points with uniform probability to be sent in any direction. This function
#' was written to apply a blur in positions of RCP spaced clusters. The radial
#' distribution is actually "folded normal" due to the fact that r > 0. See
#' \url{https://en.wikipedia.org/wiki/Folded_normal_distribution}.
#'
#' @param n The number of points you want returned
#' @param mean Mean of the regular normal distribution for radius.
#' @param sd Standard deviation for the normal distribution for radius.
#' @param coords Return comand. Either "rec" or "sph". See below for more.
#' @param method 1 for uniformly distributing direction and normally
#' distributing r, 2 for normally distributing x y and z
#'
#' @return If \code{coords = "rec"}, returns a vector of cartesian coordinates.
#'   If \code{coords = "sph"}, returns a vector of spherical coordinates.

rgblur <- function(n = 1,mean = 0,sd = 1, coords = "rec", method = 1){

  if(method == 1){
    r <- abs(rnorm(n,mean,sd))
    theta <- runif(n,0,2*pi)
    phi <- acos(runif(n,-1,1))

    if(coords == "sph"){
      return(as.data.frame(cbind(r,theta,phi)))
    }else{
      x <- r*sin(phi)*cos(theta)
      y <- r*sin(phi)*sin(theta)
      z <- r*cos(phi)
      return(as.data.frame(cbind(x,y,z)))
    }
  }else if(method == 2){
    x <- rnorm(n, mean, sd)
    y <- rnorm(n, mean, sd)
    z <- rnorm(n, mean, sd)
    return(as.data.frame(cbind(x,y,z)))
  }
}


#### Post 1/21/2020 helpers: ####
#### pp3_scale ####
# Scale pp3 pattern according to scale factor
pp3_scale <- function(X, sf){
  mks <- marks(X)
  X.d <- domain(X)
  X.d.new <- box3(X.d$xrange*sf, X.d$yrange*sf, X.d$zrange*sf)
  X.coo.new <- coords(X)*sf
  X.new <- pp3(X.coo.new$x, X.coo.new$y, X.coo.new$z, X.d.new)
  marks(X.new) <- mks
  return(X.new)
}

#### over_cut ####
# Get correctly sized (size win) and positioned over rcp
# win is the dimensions of the under point pattern
over_cut <- function(X, win, cr.rand){
  win.full <- win + 2*max(cr.rand)
  add <- max(cr.rand)
  X.cut <- subSquare(X, win.full)
  mks <- marks(X.cut)
  X.shift.coo <- coords(X.cut) - add
  X.shift.win <- box3(c(-add, win[1] + add), c(-add, win[2] + add), c(-add, win[3] + add))
  X.fin <- pp3(X.shift.coo$x, X.shift.coo$y, X.shift.coo$z, X.shift.win)
  marks(X.fin) <- mks
  return(X.fin)
}

#### bdist.complex ####
# Return distance to closest edge of window, and whether it is inside (+) or
# outside (-) that window
bdist.complex <- function(X.df, win){

  x <- X.df$x
  y <- X.df$y
  z <- X.df$z

  xmin <- min(win$xrange)
  xmax <- max(win$xrange)
  ymin <- min(win$yrange)
  ymax <- max(win$yrange)
  zmin <- min(win$zrange)
  zmax <- max(win$zrange)

  result <- data.frame(x = pmin.int(x - xmin, xmax - x),
                       y = pmin.int(y - ymin, ymax - y),
                       z = pmin.int(z - zmin, zmax - z))

  return(result)
}

#### pb_gen ####
# Generate random position shifts
pb_gen <- function(n, mean, sd){
  r <- abs(rnorm(n,mean,sd))
  theta <- runif(n,0,2*pi)
  phi <- acos(runif(n,-1,1))
  x <- r*sin(phi)*cos(theta)
  y <- r*sin(phi)*sin(theta)
  z <- r*cos(phi)
  return(as.data.frame(cbind(x,y,z)))
}

#### overlap_fix ####
# Shift points so none of the clusters overlap after a position blur
overlap_fix <- function(X, cr.rand){
  nnd <- nndist(X)
  nnw <- nnwhich(X)

  check <- which(nnd < (cr.rand + cr.rand[nnw]))
  #print(check)

  t1 <- Sys.time()
  while(!is.empty(check)){
    ind <- check[1]
    minsep <- cr.rand[ind] + cr.rand[nnw[ind]]
    direction <- (coords(X)[nnw[ind],]-coords(X)[ind,])/nnd[ind]
    coords(X)[ind,] <- coords(X)[ind,]+(nnd[ind]-(minsep+0.00001))*direction
    nnd <- nndist.pp3(X)
    nnw <- nnwhich.pp3(X)
    check <- which(nnd < (cr.rand + cr.rand[nnw]))
    t2 <- Sys.time()
    if((as.numeric(t2) - as.numeric(t1)) > 15){
      print('Impossible to find no-overlap solution.')
      return(-1)
    }
  }
  return(X)
}

#### check_vol ####
# Function to check volume in clusters and optimize
check_vol <- function(sf, coo, vol.target, under.domain, cr.rand){
  coo.scaled <- coo*sf

  bdist <- bdist.complex(coo.scaled, under.domain)
  io <- as.numeric(inside.boxx(coo.scaled$x, coo.scaled$y, coo.scaled$z,
                               w = under.domain))

  fullin <- apply(bdist > cr.rand, 1, all) & io == 1

  partin.face <- !apply(bdist > cr.rand, 1, all) &
    apply(bdist < cr.rand , 1, sum) == 1 & io == 1
  partin.edge <- !apply(bdist > cr.rand, 1, all) &
    apply(bdist < cr.rand , 1, sum) == 2 & io == 1
  partin.corner <- !apply(bdist > cr.rand, 1, all) &
    apply(bdist < cr.rand , 1, sum) == 3 & io == 1
  partin.face.h <- cr.rand[partin.face] - apply(bdist[partin.face,], 1, min)
  partin.edge.h <- cr.rand[partin.edge] - t(apply(bdist[partin.edge,], 1, function(x){sort(x)[1:2]}))
  partin.corner.h <- cr.rand[partin.corner] - t(apply(bdist[partin.corner,], 1, function(x){sort(x)}))

  fullout <- apply(bdist < -cr.rand, 1, any) & io == 0

  partout.face <- !apply(bdist < -cr.rand, 1, any) &
    apply(bdist < 0 , 1, sum) == 1 & io == 0
  partout.edge <- !apply(bdist < -cr.rand, 1, any) &
    apply(bdist < 0 , 1, sum) == 2 & io == 0
  partout.corner <- !apply(bdist < -cr.rand, 1, any) &
    apply(bdist < 0 , 1, sum) == 3 & io == 0
  partout.face.h <- cr.rand[partout.face] + apply(bdist[partout.face,], 1, min)
  partout.edge.h <- cr.rand[partout.edge] + apply(bdist[partout.edge,], 1, min)
  partout.corner.h <- cr.rand[partout.corner] + apply(bdist[partout.corner,], 1, min)

  #fullin volume
  vol.fullin <- 4/3*pi*cr.rand[fullin]^3

  #partin volume
  vol.partin.face <- 4/3*pi*cr.rand[partin.face]^3 - 1/3*pi*partin.face.h^2*(3*cr.rand[partin.face] - partin.face.h)
  if(sum(partin.edge) > 0){
    vol.partin.edge <- 4/3*pi*cr.rand[partin.edge]^3 -
      1/3*pi*partin.edge.h[,1]^2*(3*cr.rand[partin.edge] - partin.edge.h[,1]) -
      1/3*pi*partin.edge.h[,2]^2*(3*cr.rand[partin.edge] - partin.edge.h[,2])
  }else{
    vol.partin.edge <- c(0)
  }
  if(sum(partin.corner) > 0){
    vol.partin.corner <- 4/3*pi*cr.rand[partin.corner]^3 -
      1/3*pi*partin.corner.h[,1]^2*(3*cr.rand[partin.corner] - partin.corner.h[,1]) -
      1/3*pi*partin.corner.h[,2]^2*(3*cr.rand[partin.corner] - partin.corner.h[,2]) -
      1/3*pi*partin.corner.h[,3]^2*(3*cr.rand[partin.corner] - partin.corner.h[,3])
  }else{
    vol.partin.corner <- c(0)
  }

  #partout volume
  vol.partout.face <- 1/3*pi*partout.face.h^2*(3*cr.rand[partout.face] - partout.face.h)
  vol.partout.edge <- 0.5*(1/3*pi*partout.edge.h^2*(3*cr.rand[partout.edge] - partout.edge.h))
  vol.partout.corner <- 0.25*(1/3*pi*partout.corner.h^2*(3*cr.rand[partout.corner] - partout.corner.h))

  #fullout volume = zero

  # add to get total volume
  vols.total <- sum(vol.fullin) +
    sum(vol.partin.face) + sum(vol.partin.edge) + sum(vol.partin.corner) +
    sum(vol.partout.face) + sum(vol.partout.edge) + sum(vol.partout.corner)

  return(abs(vol.target - vols.total))
}
