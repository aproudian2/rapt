#
# This file contains spatial processing functions for APT data
#
rownames.pp3 <- function(pat) {
  dat <- rownames(as.data.frame(pat))
  return(dat)
}

## Shape Subsets ##
cylinder <- function(x, ...) {
  UseMethod('cylinder', x)
}
cylinder.pp3 <- function(pat, ind, r = 1) {
  cylX <- pat$data$x[ind]
  cylY <- pat$data$y[ind]
  dat <- subset(pat, subset = r > sqrt((x-cyl)^2 + (y-cylY)^2))
  return(dat)
}
cone <- function(x, ...) {
  UseMethod('cone', x)
}
cone.pp3 <- function(pat, r = 1, phi = pi/4) {
  coneX <- as.numeric(pat$data$x[1])
  coneY <- as.numeric(pat$data$y[1])
  coneZ <- as.numeric(pat$data$z[1])
  a <- tan(phi)
  k <- a * r
  dat <- subset(pat, subset =
    z >=  a * sqrt((x-coneX)^2 + (y-coneY)^2) + coneZ - k |
    z <= -a * sqrt((x-coneX)^2 + (y-coneY)^2) + coneZ + k
  )
  return(dat)
}
cone.ppp <- function(pat, r = 1, phi = pi/4) {
  coneX <- pat$x[1]
  coneY <- pat$y[1]
  a <- tan(phi)
  k <- a * r
  dat <- subset(pat, subset =
    y >  a * sqrt((x-coneX)^2) + coneY - k |
    y < -a * sqrt((x-coneX)^2) + coneY + k
  )
  return(dat)
}

pointDist <- function(pat) {
  orig <- as.data.frame(pat[1])
  pat <- as.data.frame(pat)
  dat <- apply(pat, 1, function(point) {
     sqrt(sum((point-orig)^2))
   })
  dat <- as.numeric(dat)
  dat <- dat[dat > 0]
  return(dat)
}
coneRDF <- function(x, ...) {
  UseMethod('coneRDF', x)
}
coneRDF.ppp <- function(pat, cl, r = 1, phi = pi/4, nrval = 128) {

}
coneRDF.pp3 <- function(pat, cl, r = 1, phi = pi/4, nrval = 128) {
  rdf.n <- npoints(pat)
  rdf.ind <- 1:rdf.n
  rdf.dist <- parLapply(cl, rdf.ind, function(ind){
    subpat <- pat[ind:rdf.n]
    pointDist(cone(subpat, r = r, phi = phi))
  })
  rdf.dist <- unlist(rdf.dist)
  rdf.breaks <- seq(0, max(rdf.dist), length.out = nrval+1)
  rdf.hist <- hist(rdf.dist, breaks = rdf.breaks, plot = F)
  rdf.r <- rdf.hist$mids
  rdf.dr <- diff(rdf.breaks)[1]
  rdf.a <- tan(phi)
  rdf.density <- intensity.ppx(pat) # Improve accuracy
  rdf.dvol <- 4*pi * (1 - 1/sqrt(1 + (rdf.a + rdf.a^2)^2)) * rdf.dr
  rdf.g <- rdf.hist$counts / (rdf.density * rdf.dvol * rdf.r^2)
  rdf.dat <- data.frame(r = rdf.r, g = rdf.g)
  return(rdf.dat)
}

# Extends marktable (from 'Spatstat') to pp3
marktable.pp3 <- function (X, R, N, exclude = TRUE, collapse = FALSE)
{
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

# Extends the superimpose function from "SpatStat" to handle pp3
## Add ability to superimpose with marks
superimpose.pp3 <- function(..., W = NULL, check = F) {
  input.list <- list(...)
  df.list <- lapply(input.list, as.data.frame)
  df.comb <- Reduce(rbind, df.list)
  out.pp3 <-  createSpat(df.comb, win = W)
  return(out.pp3)
}
# Extends the shift function from "SpatStat" to handle pp3
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
# Extends the inside method to the pp3 class
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
# Extends the sample function from "base" to handle pp3
sample.pp3 <- function(X, size) {
  sam.lab <- rownames(as.data.frame(X$data))
  sam.pts <- sample(sam.lab, size)
  sam.dat <- X[sam.pts]
  return(sam.dat)
}
# Finds clusters by NN adjacency marks
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
# Extends intensity (from 'Spatstat') to pp3
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
