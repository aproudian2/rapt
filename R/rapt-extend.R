#### marktable.pp3 ####
#' Extends \code{\link[spatstat]{marktable}} to \code{\link[spatstat]{pp3}}.
#'
#' \code{marktable.pp3}
#' @seealso \code{\link[spatstat]{marktable}}
#' @export
marktable.pp3 <- function (X, R, N, exclude = TRUE, collapse = FALSE) {
  verifyclass(X, "pp3")
  if (!is.marked(X, dfok = FALSE))
    stop("point pattern has no marks")
  gotR <- !missing(R) && !is.null(R)
  gotN <- !missing(N) && !is.null(N)
  if (gotN == gotR)
    stop("Exactly one of the arguments N and R should be given")
  stopifnot(is.logical(exclude) && length(exclude) == 1)
  m <- marks(X)
  if (!is.factor(m))
    stop("marks must be a factor")
  if (gotR) {
    stopifnot(is.numeric(R) && length(R) == 1 && R > 0)
    p <- closepairs.pp3(X, R, what = "indices")
    pi <- p$i
    pj <- p$j
    if (!exclude) {
      n <- npoints(X)
      pi <- c(pi, 1:n)
      pj <- c(pj, 1:n)
    }
  }
  else {
    stopifnot(is.numeric(N) && length(N) == 1)
    ii <- seq_len(npoints(X))
    nn <- nnwhich(X, k = 1:N)
    if (N == 1)
      nn <- matrix(nn, ncol = 1)
    if (!exclude)
      nn <- cbind(ii, nn)
    pi <- as.vector(row(nn))
    pj <- as.vector(nn)
  }
  if (!collapse) {
    i <- factor(pi, levels = seq_len(npoints(X)))
    mj <- m[pj]
    mat <- table(point = i, mark = mj)
  }
  else {
    mi <- m[pi]
    mj <- m[pj]
    mat <- table(point = mi, neighbour = mj)
  }
  return(mat)
}

#### superimpose.pp3 ####
#' Extends \code{\link[spatstat]{superimpose}} to \code{\link[spatstat]{pp3}}.
#'
#' \code{superimpose.pp3}
#' @seealso \code{\link[spatstat]{superimpose}}
#' @export
superimpose.pp3 <- function(..., W = NULL, check = F) {
  # Add ability to superimpose with marks
  input.list <- list(...)
  df.list <- lapply(input.list, as.data.frame)
  df.comb <- Reduce(rbind, df.list)
  out.pp3 <-  createSpat(df.comb, win = W)
  return(out.pp3)
}

#### shift.pp3 ####
#' Extends \code{\link[spatstat]{shift}} to \code{\link[spatstat]{pp3}}.
#'
#' @seealso \code{\link[spatstat]{shift}}
#' @export
shift.pp3 <- function (X, vec = c(0, 0, 0), ..., origin = NULL)
{
  verifyclass(X, "pp3")
  if (!is.null(origin)) {
    if (!missing(vec))
      warning("argument vec ignored; overruled by argument origin")
    if (is.numeric(origin)) {
      locn <- origin
    }
    else if (is.character(origin)) {
      origin <- pickoption("origin", origin,
                           c(midpoint = "midpoint", bottomleft = "bottomleft"))
      W <- X$domain
      locn <- switch(origin, midpoint = {
        c(mean(W$domain$xrange), mean(W$domain$yrange), mean(W$domain$zrange))
      }, bottomleft = {
        c(W$domain$xrange[1], W$domain$yrange[1], W$domain$zrange[1])
      })
    }
    else stop("origin must be a character string or a numeric vector")
    vec <- -locn
  }
  Y <- pp3(X$data$x + vec[1], X$data$y + vec[2], X$data$z + vec[3],
           xrange = X$domain$xrange + vec[1],
           yrange = X$domain$yrange + vec[2],
           zrange = X$domain$zrange + vec[3])
  attr(Y, "lastshift") <- vec
  return(Y)
}
#### inside.pp3 ####
#' Extends \code{\link[spatstat]{inside}} to \code{\link[spatstat]pp3}}.
#'
#' @seealso \code{\link[spatstat]{inside}}
#' @export
inside.pp3 <- function(points, domain = NULL) {
  if(is.null(domain)) {
    domain <- points$domain
  }
  if (length(points) == 0)
    return(logical(0))
  xr <- domain$xrange
  yr <- domain$yrange
  zr <- domain$zrange
  x <- points$data$x
  y <- points$data$y
  z <- points$data$z
  eps <- sqrt(.Machine$double.eps)
  frameok <- (x >= xr[1] - eps) & (x <= xr[2] + eps) &
    (y >= yr[1] - eps) & (y <= yr[2] + eps) &
    (z >= zr[1] - eps) & (z <= zr[2] + eps)
  return(frameok)
}
#### sample.ppp ####
#' Extends \code{\link[base]{sample}} to handle \code{\link[spatstat]{ppp}}.
#'
#' @param X A \code{ppp}. The point pattern from which to sample.
#' @param size A numeric. The number of points to sample.
#' @return A \code{ppp}. The sampled point pattern.
#'
#' @seealso \code{\link[base]{sample}}
#' @export
sample.ppp <- function(X, size) {
  sam.n <- npoints(X)
  sam.pts <- sample(1:sam.n, size)
  sam.dat <- X[sam.pts]
  return(sam.dat)
}
#### sample.pp3 ####
#' Extends \code{\link[base]{sample}} to handle \code{\link[spatstat]{pp3}}.
#'
#' @param X A \code{pp3}. The point pattern from which to sample.
#' @param size A numeric. The number of points to sample.
#' @return A \code{pp3}. The sampled point pattern.
#'
#' @seealso \code{\link[base]{sample}}
#' @export
sample.pp3 <- function(X, size) {
  sam.lab <- rownames(as.data.frame(X$data))
  sam.pts <- sample(sam.lab, size)
  sam.dat <- X[sam.pts]
  return(sam.dat)
}
#### findClusters.pp3 ####
#' Finds clusters by NN adjacency marks.
findClusters.pp3 <- function(X, mark, k = 1) {
  if(!(mark %in% marks(X)))
    stop('The specified mark does not exist in the pattern')
  else if(length(mark) != 1)
    stop('Only one mark may be tested')
  spl <- split(X)[[mark]]
  # required because pp3 functions don't work on marked pp3 by default
  class(spl) <- c('pp3', 'ppx')
  ind <- rownames.pp3(spl)
  fCl <- lapply(ind, function(x) {
    clus <- x
    rem <- rownames.pp3(X)
    nn <- 0
    while(length(nn) > 0) {
      rem <- setdiff(rem, nn)
      nn <- nncross(spl[clus], X[rem], k = 1:k, what = 'which')
      nn <- unlist(nn)
      nn <- rownames.pp3(X[rem])[nn]
      nn <- nn[marks(X[nn]) == mark]
      clus <- append(clus, nn)
    }
    clus <- sort(unique(clus))
    return(clus)
  })
  return(fCl)
}
#### intensity.pp3 ####
#' Extends \code{\link[spatstat]{intensity}} to \code{\link[spatstat]{pp3}}.
#'
#' @seealso \code{\link[spatstat]{intensity}}
#' @export
intensity.pp3 <- function(X, weights = NULL) {
  n <- npoints(X)
  a <- volume(domain(X))
  if (is.null(weights)) {
    if (is.multitype(X)) {
      mks <- marks(X)
      answer <- as.vector(table(mks))/a
      names(answer) <- levels(mks)
    }
    else answer <- n/a
    return(answer)
  } else
    warning('Weights is not yet implemented for pp3')
}

#### rownames.pp3 ####
#' Extends \code{\link[base]{rownames}} to \code{\link[spatstat]{pp3}}.
#'
#' @param pat A \code{pp3}. The point pattern from which to extract rownames.
#' @return A string vector. The rownames of the point pattern.
#' @seealso \code{\link[base]{rownames}}
rownames.pp3 <- function(pat) {
  dat <- rownames(as.data.frame(pat))
  return(dat)
}
#### plot3d.pp3 ####
#' Plot a \code{\link[spatstat]{pp3}} in a manipulatable 3D plot.
#'
#' (requires the rgl library)
#' @param X A \code{pp3}. The point pattern to visualize
#' @param ... Other arguments to pass to \code{plot3d} from the \code{rgl}
#' library.
#' @seealso \code{\link[rgl]{plot3d}}
plot3d.pp3 <- function(X, ...) {
  rgl::plot3d(as.data.frame(X$data), ...)
}

<<<<<<< HEAD
#### localK3est ####
#' Extends \code{\link[spatstat]{localK}} to 3D
#'
#' Outputs local K function information for a \code{\link[spatstat]{pp3}}
#' pattern. Currently only implemented for the translation edge correction.
#'
#' @param X A \code{\link[spatstat]{pp3}} object to test.
#' @param rmax See \code{\link[spatstat]{K3est}}.
#' @param nrval See \code{\link[spatstat]{K3est}}.
#' @param correction Currently only "translation" implemented.
#'
#' @return Data frame with columns for the local K function around each point in
#'   the pattern.
localK3est <- function(X, rmax=NULL, nrval=128, correction="translation") {

  stopifnot(inherits(X, "pp3"))
  correction <- pickoption("correction", correction,
                           c(translation="translation",
                             trans="translation",
                             isotropic="isotropic",
                             iso="isotropic",
                             best="isotropic"),
                           multi=TRUE)

  if(correction != "translation"){
    print("Local K3est function is only implemented for the translation edge correction.")
    return()
  }

  B <- X$domain
  if(is.null(rmax))
    rmax <- diameter(B)/2
  r <- seq(from=0, to=rmax, length.out=nrval)
  np <- npoints(X)

  # extract the x,y,z ranges as a vector of length 6
  flatbox <- unlist(B[1:3])

  # extract coordinates
  coo <- coords(X)

  u <- localk3engine(coo$x, coo$y, coo$z, flatbox, rmax=rmax, nrval=nrval)
  um <- matrix(u, nrow = nrval, ncol = np)
  kavg <- apply(um, 1, mean)

  rseq <- seq(0, rmax, len = nrval)
  rtheo <- (4/3)*pi*rseq^3
  K <- matrix(c(rseq,rtheo,kavg,u), nrow = nrval, ncol = (np+3), byrow = FALSE)
  K <- as.data.frame(K)

  names <- vector("character", np+3)
  names[1] <- "r"
  names[2] <- "theo"
  names[3] <- "kavg"
  for(i in 1:np){
    names[i+3] <- paste("k",toString(i), sep = "")
  }

  names(K) <- names

  return(K)
}

#### localk3engine ####
#' Wrapper function for the C code under \code{\link{localK3est}}.
localk3engine <- function(x, y, z, box=c(0,1,0,1,0,1), rmax=1, nrval=100){
  res <- .C("RcallK3local",
            as.double(x), as.double(y), as.double(z),
            as.integer(length(x)),
            as.double(box[1L]), as.double(box[2L]),
            as.double(box[3L]), as.double(box[4L]),
            as.double(box[5L]), as.double(box[6L]),
            as.double(0), as.double(rmax),
            as.integer(nrval),
            f = as.double(numeric(nrval)),
            num = as.double(numeric(nrval)),
            denom = as.double(numeric(nrval)),
            full = as.double(numeric(nrval*length(x))),
            PACKAGE = "rapt")
  return(res$full)
}

#### quadratcount.pp3 ####
#' Extension of \code{\link[spatstat]{quadrat}} to \code{\link[spatstat]{pp3}} objects.
#'
#' Divides volume into quadrats and counts the number of points in each quadrat.
#'
#' @param X The \code{\link[spatstat]{pp3}} object to split up.
#' @param nx,ny,nz Number of ractangular quadrats in the x, y, and z directions.
#'
#' @return A \code{data.frame} object containing the number of counts in each quadrat.
quadratcount.pp3 <- function(X, nx = 5, ny = 5, nz = 5){
  verifyclass(X, "pp3")
  w <- domain(X)

  # create box3objects for each quadrat
  xlim <- w$xrange
  ylim <- w$yrange
  zlim <- w$zrange

  xbreaks <- seq(xlim[1],xlim[2],length.out = (nx+1))
  ybreaks <- seq(ylim[1],ylim[2],length.out = (ny+1))
  zbreaks <- seq(zlim[1],zlim[2],length.out = (nz+1))

  ntot <- nx*ny*nz
  gridvals <- list()
  cnt <- 1

  for(i in 1:nx){
    for(j in 1:ny){
      for(k in 1:nz){
        gridvals[[cnt]] <- box3(xrange = xbreaks[i:(i+1)], yrange = ybreaks[j:(j+1)], zrange = zbreaks[k:(k+1)])
        cnt <- cnt + 1
      }
    }
  }

  inside.tf <- lapply(gridvals, function(x){inside.boxx(X, w = x)})
  counts <- lapply(inside.tf, function(x){sum(x)})
  counts <- unlist(counts)
  return(data.frame(quad.no = seq(1,ntot), count = counts))
}

#### nndensity.pp3 ####
#' Extension of \code{\link[spatstat]{nndensity}} to handle pp3 objects.
#'
#' Calculates the 3D nearest-neighbor intensity estimate of a point process at
#' either a grid of points or at the point locations in the data set. Utilizes
#' the volume weighted edge correction. See Statistics for Spatial Data by
#' Cressie pg. 654 for more info.
#'
#' @param X The point pattern to estimate the intensity of.
#' @param k Vector containing the nearest-neighbor #s that the estimate should
#'   be calculated for.
#' @param nx,ny,nz If estimating on a grid, the number of grid points in x, y,
#'   and z.
#' @param dz The z spacing for numeric integration for the edge correction
#'   (suggested ~ 0.1-0.5)
#' @param at.points \code{TRUE} or \code{FALSE}. Whether or not to estimate
#'   intensity at points in pattern. If \code{TRUE}, nx, ny, and nz are not
#'   used.
#' @param par \code{TRUE} or \code{FALSE}: whether or not to calculate in
#'   parallel
#' @param cores If \code{par = TRUE}, this is the number of cores to use for the
#'   parallel calculation.
#'
#' @return List containing: [[1]] A data frame of the intensity estimates for
#'   each nearest neighbor value. [[2]] The coordinates of the estimates. [[3]]
#'   The coordinates of the original points from the data set.
nndensity.pp3 <- function(X, k, nx, ny, nz, dz, at.points = FALSE, par = TRUE, cores = 7){

  if(at.points == FALSE){
    # set up grid of points and find nearest neighbors from grid to data set
    d <- domain(X)

    xsep <- (d$xrange[2]-d$xrange[1])/nx
    xstart <- d$xrange[1] + xsep/2
    xend <- d$xrange[2] - xsep/2

    ysep <- (d$yrange[2]-d$yrange[1])/ny
    ystart <- d$yrange[1] + ysep/2
    yend <- d$yrange[2] - ysep/2

    zsep <- (d$zrange[2]-d$zrange[1])/nz
    zstart <- d$zrange[1] + zsep/2
    zend <- d$zrange[2] - zsep/2

    x <- seq(xstart,xend, xsep)
    y <- seq(ystart,yend, ysep)
    z <- seq(zstart,zend, zsep)

    coo <- expand.grid(x,y,z)
    names(coo) <- c('x','y','z')
    grid <- pp3(coo$x, coo$y, coo$z, domain(X))

    nnk <- nncross(grid, X, what = "dist", k = k)
    bdist <- bdist.points3.multi(grid) # distances to nearest boundaries
    est.points <- coo

  }else{
    #find nearest neighbors from data set to itself
    nnk <- nndist(X, k = k)
    bdist <- bdist.points3.multi(X) # distances to nearest boundaries
    est.points <- coords(X)
  }

  lambda.est <- local.den.engine(bdist, nnk, k, dz, par, cores)

  res <- list(lambda.est = lambda.est, estimate.coords = est.points, x = coords(X))

  return(res)
}

#### nncrossden.pp3 ####
#' Similar to \code{\link{nndensity.pp3}}, but specifically for marked patterns
#' with a global intensity inhomogeneity.
#'
#' Calculates the 3D nearest-neighbor intensity estimate of a point process at
#' either a grid of points or at the point locations in the data set. Requires
#' two point patterns: One (\code{X}) that contains just points of the type you
#' want to estimate the intesnity of. A second (\code{Y}) that contains all
#' points from the sample (note that this includes the points from pattern
#' \code{X}) for removing global intensity. Note that this function runs faster
#' on linux and mac machines than windows, but will work on both.
#'
#' @param X The point pattern to estimate the intensity of.
#' @param Y The full point pattern to estimate global intensity.
#' @param k Vector containing the nearest-neighbor #s that the estimate should
#'   be calculated for.
#' @param nx,ny,nz If estimating on a grid, the number of grid points in x, y,
#'   and z.
#' @param at.points \code{TRUE} or \code{FALSE}. Whether or not to estimate
#'   intensity at points in pattern. If \code{TRUE}, nx, ny, and nz are not
#'   used.
#' @param nsplit In this function, large \code{pp3} objects are split into
#'   multiple smaller data sets so that the memory is not overloaded while doing
#'   computations. This parameter is the number of points per split set.
#' @param cores Number of cores to use for parallelization. Set to one for
#'   serial calculation.
#' @param os Either 'windows', 'mac', or 'linux'. Changes the parallelization
#'   method used.
#' @return List containing: [[1]] A data frame of the intensity estimates for
#'   each nearest neighbor value. [[2]] The coordinates of the estimates. [[3]]
#'   The coordinates of the original points from the data set. [[4]] A vector
#'   containing all k values tested.
nncrossden.pp3 <- function(X, Y, k, nx, ny, nz, at.points = FALSE, nsplit = 1000, cores = 8, os = "linux"){
  # yes this is complicated as hell... but it works I promise -GV

  t1 <- Sys.time()
  if(at.points == FALSE){
    # set up grid of points
    d <- domain(X)

    xsep <- (d$xrange[2]-d$xrange[1])/nx
    xstart <- d$xrange[1] + xsep/2
    xend <- d$xrange[2] - xsep/2

    ysep <- (d$yrange[2]-d$yrange[1])/ny
    ystart <- d$yrange[1] + ysep/2
    yend <- d$yrange[2] - ysep/2

    zsep <- (d$zrange[2]-d$zrange[1])/nz
    zstart <- d$zrange[1] + zsep/2
    zend <- d$zrange[2] - zsep/2

    x <- seq(xstart,xend, xsep)
    y <- seq(ystart,yend, ysep)
    z <- seq(zstart,zend, zsep)

    coo <- expand.grid(x,y,z)
    names(coo) <- c('x','y','z')

    grid.n <- nrow(coo)
    if(grid.n > nsplit){
      coo.split.ind <- split(1:nrow(coo), ceiling(1:nrow(coo)/nsplit))
      coo.split <- lapply(coo.split.ind, function(x){coo[x,]})
      grid.split <- lapply(coo.split, function(x){pp3(x$x, x$y, x$z, domain(X))})
    }else{
      grid.split <- list(pp3(coo$x, coo$y, coo$z, domain(X)))
      coo.split.ind <- list(1:grid.n)
    }

    est.points <- coo

  }else{
    grid.n <- npoints(X)
    if(grid.n > nsplit){
      coo <- coords(X)
      coo.split.ind <- split(1:nrow(coo), ceiling(1:nrow(coo)/nsplit))
      coo.split <- lapply(coo.split.ind, function(x){coo[x,]})
      grid.split <- lapply(coo.split, function(x){pp3(x$x, x$y, x$z, domain(X))})
    }else{
      grid.split <- list(X)
      coo.split.ind <- list(1:grid.n)
    }

    est.points <- coords(X)
  }

  lambda.global.Y <- npoints(Y)/volume(domain(Y))

  if(length(k) > 1){
    lambda.est <- matrix(NA, nrow = grid.n, ncol = length(k))

    if(os == "windows"){
      cl <- makePSOCKcluster(cores)
      clusterExport(cl, "nncross")
      clusterExport(cl, c("X", "k"), envir = environment())
    }else{
      cl <- makeForkCluster(cores)
    }

    nnk.X.split <- parLapply(cl, grid.split, function(x){nncross(x, X, what = "dist", k = k)})

    if(os == "windows"){
      clusterExport(cl,c("lambda.global.Y","nnk.X.split"), envir = environment())
    }

    tot <- length(k) * length(grid.split)
    cnt <- 0

    for(i in 1:length(k)){
      for(j in 1:length(grid.split)){
        cp.Y <- crosspairs(grid.split[[j]], Y, rmax = max(nnk.X.split[[j]]), what = "ijd")
        cp.Y <- data.frame(i = cp.Y[[1]], j = cp.Y[[2]], dist = cp.Y[[3]])
        cp.Y.list <- split(cp.Y, factor(cp.Y$i))
        xi <- 1:nrow(nnk.X.split[[j]])
        lambda.est[coo.split.ind[[j]],i] <- parSapply(cl, xi, function(xq,i,cp.Y.list,j){(k[i]/sum(cp.Y.list[[xq]]$dist < nnk.X.split[[j]][[xq,i]]))*lambda.global.Y},i,cp.Y.list,j)
        cnt <- cnt + 1
        print(paste(toString(round(100*cnt/tot, 1)), "%", sep = ""))
      }
    }
    stopCluster(cl)

  }else{
    lambda.est <- matrix(NA, nrow = grid.n, ncol = 1)

    if(os == "windows"){
      cl <- makePSOCKcluster(cores)
      clusterExport(cl, "nncross")
      clusterExport(cl, c("X", "k"), envir = environment())
    }else{
      cl <- makeForkCluster(cores)
    }

    nnk.X.split <- parLapply(cl, grid.split, function(x){nncross(x, X, what = "dist", k = k)})

    if(os == "windows"){
      clusterExport(cl,c("lambda.global.Y","nnk.X.split"), envir = environment())
    }

    for(i in 1:length(grid.split)){
      cp.Y <- crosspairs(grid.split[[i]], Y, rmax = max(nnk.X.split[[i]]), what = "ijd")
      cp.Y <- data.frame(i = cp.Y[[1]], j = cp.Y[[2]], dist = cp.Y[[3]])
      cp.Y.list <- split(cp.Y, factor(cp.Y$i))
      xi <- 1:length(nnk.X.split[[i]])

      lambda.est[coo.split.ind[[i]],1] <- parSapply(cl, xi, function(xq, i, cp.Y.list){(k/sum(cp.Y.list[[xq]]$dist < nnk.X.split[[i]][[xq]]))*lambda.global.Y}, i, cp.Y.list)

      print(paste(toString(round(100*i/length(grid.split), 1)), "%", sep = ""))
    }
    stopCluster(cl)
  }

  res <- list(lambda.est = lambda.est, estimate.coords = est.points, x = coords(X))

  res <- as.data.frame(lambda.est)
  names <- sapply(k,function(x){return(paste("nn",toString(x),sep = ""))})
  colnames(res) <- names

  t2 <- Sys.time()
  print("Time to complete:")
  print(t2-t1)
  return(list(lambda.est = res, estimate.coords = est.points, x = coords(X), k = k))
}

=======
>>>>>>> a54e6d28facbb2490101091bb7211c58558492b0
#### K3cross ####
K3multi <- function(X, I, J, r, breaks,
              correction = c("none", "isotropic", "translation"),
              ..., ratio = FALSE) {
  verifyclass(X, "pp3")
  npts <- npoints(X)
  W <- X$domain
  volW <- volume(W)
  I <- ppsubset(X, I)
  J <- ppsubset(X, J)
  if (is.null(I) || is.null(J))
    stop("I and J must be valid subset indices")
  if (!any(I))
    stop("no points belong to subset I")
  if (!any(J))
    stop("no points belong to subset J")
  nI <- sum(I)
  nJ <- sum(J)
  lambdaI <- nI/volW
  lambdaJ <- nJ/volW
  rmax <- max(r)
  alim <- c(0, rmax)
  K <- data.frame(r = r, theo = (4/3) * pi * r^3)
  desc <- c("distance argument r", "theoretical Poisson %s")
  K <- fv(K, "r", quote(K[IJ](r)), "theo", , alim, c("r", "{%s[%s]^{pois}}(r)"),
          desc, fname = c("K", "list(I,J)"), yexp = quote(K[list(I,
                                                                 J)](r)))
  if (ratio) {
    denom <- lambdaI * lambdaJ * volW
    numK <- eval.fv(denom * K)
    denK <- eval.fv(denom + K * 0)
    attributes(numK) <- attributes(denK) <- attributes(K)
    attr(numK, "desc")[2] <- "numerator for theoretical Poisson %s"
    attr(denK, "desc")[2] <- "denominator for theoretical Poisson %s"
  }
  XI <- X[I]
  XJ <- X[J]
  close <- crosspairs(XI, XJ, max(r))
  orig <- seq_len(npts)
  imap <- orig[I]
  jmap <- orig[J]
  iX <- imap[close$i]
  jX <- jmap[close$j]
  if (any(I & J)) {
    ok <- (iX != jX)
    if (!all(ok)) {
      close$i <- close$i[ok]
      close$j <- close$j[ok]
      close$d <- close$d[ok]
    }
  }
  dcloseIJ <- close$d
  icloseI <- close$i
  jcloseJ <- close$j
  if (any(correction == "none")) {
    wh <- whist(dcloseIJ, breaks)
    numKun <- cumsum(wh)
    denKun <- lambdaI * lambdaJ * volW
    Kun <- numKun/denKun
    K <- bind.fv(K, data.frame(un = Kun), "{hat(%s)[%s]^{un}}(r)",
                 "uncorrected estimate of %s", "un")
    if (ratio) {
      numK <- bind.fv(numK, data.frame(un = numKun), "{hat(%s)[%s]^{un}}(r)",
                      "numerator of uncorrected estimate of %s", "un")
      denK <- bind.fv(denK, data.frame(un = denKun), "{hat(%s)[%s]^{un}}(r)",
                      "denominator of uncorrected estimate of %s",
                      "un")
    }
  }
  formula(K) <- . ~ r
  unitname(K) <- unitname(X)
  if (ratio) {
    formula(numK) <- formula(denK) <- . ~ r
    unitname(numK) <- unitname(denK) <- unitname(K)
    K <- rat(K, numK, denK, check = FALSE)
  }
  return(K)
}
