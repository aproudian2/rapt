#
# This file contains extensions to the spatstat package.
#

#### marktable ####
#' Tabulate Marks in Neighbourhood of Every Point in a Point Pattern
#'
#' @description This is an S3 generic that extends the use of
#'   \code{\link[spatstat]{marktable}} beyond \code{ppp} objects.
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat]{marktable}}, \code{\link{marktable.ppp}},
#'   \code{\link{marktable.pp3}}
#' @export
marktable <- function(X, ...) UseMethod("marktable")

### marktable.ppp ###
#' Tabulate Marks in Neighbourhood of Every Point in a Point Pattern
#'
#' @seealso \code{\link[spatstat]{marktable}}
#' @export
marktable.ppp <- spatstat::marktable

### marktable.pp3 ###
#' Tabulate Marks in Neighbourhood of Every Point in a Point Pattern
#'
#' @description Visit each point in a point pattern, find the neighbouring
#'   points, and compile a frequency table of the marks of these neighbour
#'   points.
#'
#' @family spatstat extensions
#'
#' @export
marktable.pp3 <- function (X, R, N, exclude = TRUE, collapse = FALSE) {
  spatstat::verifyclass(X, "pp3")
  if (!spatstat::is.marked(X, dfok = FALSE))
    stop("point pattern has no marks")
  gotR <- !missing(R) && !is.null(R)
  gotN <- !missing(N) && !is.null(N)
  if (gotN == gotR)
    stop("Exactly one of the arguments N and R should be given")
  stopifnot(is.logical(exclude) && length(exclude) == 1)
  m <- spatstat::marks(X)
  if (!is.factor(m))
    stop("marks must be a factor")
  if (gotR) {
    stopifnot(is.numeric(R) && length(R) == 1 && R > 0)
    p <- spatstat::closepairs(X, R, what = "indices")
    pi <- p$i
    pj <- p$j
    if (!exclude) {
      n <- spatstat::npoints(X)
      pi <- c(pi, 1:n)
      pj <- c(pj, 1:n)
    }
  }
  else {
    stopifnot(is.numeric(N) && length(N) == 1)
    ii <- seq_len(spatstat::npoints(X))
    nn <- spatstat::nnwhich(X, k = 1:N)
    if (N == 1)
      nn <- matrix(nn, ncol = 1)
    if (!exclude)
      nn <- cbind(ii, nn)
    pi <- as.vector(row(nn))
    pj <- as.vector(nn)
  }
  if (!collapse) {
    i <- factor(pi, levels = seq_len(spatstat::npoints(X)))
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

#### rjitter ####
#' Random Perturbation of a Point Pattern
#'
#' @description This is an S3 generic that extends the use of
#'   \code{\link[spatstat]{rjitter}} beyond \code{ppp} objects.
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat]{rjitter}}, \code{\link{rjitter.ppp}},
#'   \code{\link{rjitter.pp3}}
#'
#' @export
rjitter <- function(X, ...) UseMethod("rjitter")

### rjitter.ppp ###
#' Random Perturbation of a Point Pattern
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat]{rjitter}}, \code{\link{rjitter.pp3}}
#' @export
rjitter.ppp <- spatstat::rjitter

### rjitter.pp3 ###
#' Random Perturbation of a Point Pattern
#'
#' Applies independent random displacements to each point in a point pattern.
#' Extends \code{\link[spatstat]{rjitter}} to \code{\link[spatstat]{pp3}}.
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat]{rjitter}}, \code{\link{rjitter.ppp}}
#'
#' @export
rjitter.pp3 <- function(X, domain = box3()) {
  spatstat::verifyclass(X, "pp3")
  nX <- spatstat::npoints(X)
  if (nX == 0)
    return(X)
  W <- X$domain
  D <- spatstat::runifpoint3(nX, domain = domain)
  xnew <- X$data$x + D$data$x
  ynew <- X$data$y + D$data$y
  znew <- X$data$z + D$data$z
  new <- pp3(xnew, ynew, znew, W)
  ok <- subset(new, subset =
                 (x > W$xrange[1] & x < W$xrange[2]) &
                 (y > W$yrange[1] & y < W$yrange[2]) &
                 (z > W$zrange[1] & z < W$zrange[2])
  )
  return(ok)
}

#### superimpose.pp3 ####
#' Superimpose Several Geometric Patterns
#'
#' \code{superimpose.pp3} superimposes any number of 3D point patterns
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat]{superimpose}}
#'
#' @export
superimpose.pp3 <- function(..., W = NULL, check = FALSE) {
  input.list <- list(...)
  df.list <- lapply(input.list, as.data.frame)
  df.comb <- do.call(rbind, df.list)
  out.pp3 <-  createSpat(df.comb, win = W)
  if(!is.null(df.comb$marks)) {
    spatstat::marks(out.pp3) <- df.comb$marks
  }
  return(out.pp3)
}

#### shift.pp3 ####
#' Apply Vector Translation
#'
#' Applies a vector shift to a 3D point pattern
#'
#' @param X Point pattern(object of class"pp3")
#' @param vec Numeric. A vector of length 3 specifying the translation
#' @param origin Location that will be shifted to the origin. Either a vector of
#'   length 3 specifying the new origin or one of the character strings
#'   "midpoint" or "bottomleft" (will be extended)
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat]{shift}}
#'
#' @export
shift.pp3 <- function (X, vec = c(0, 0, 0), ..., origin = NULL)
{
  spatstat::verifyclass(X, "pp3")
  if (!is.null(origin)) {
    if (!missing(vec))
      warning("argument vec ignored; overruled by argument origin")
    if (is.numeric(origin) & length(origin) == 3) {
      locn <- origin
    } else if (is.character(origin)) {
      origin <- spatstat::pickoption("origin", origin,
                           c(midpoint = "midpoint", bottomleft = "bottomleft"))
      W <- X$domain
      locn <- switch(origin, midpoint = {
        c(mean(W$domain$xrange), mean(W$domain$yrange), mean(W$domain$zrange))
      }, bottomleft = {
        c(W$domain$xrange[1], W$domain$yrange[1], W$domain$zrange[1])
      })
    } else {
      stop("origin must be a character string or a numeric vector")
    }
    vec <- -locn
  }
  Y <- spatstat::pp3(X$data$x + vec[1], X$data$y + vec[2], X$data$z + vec[3],
           xrange = X$domain$xrange + vec[1],
           yrange = X$domain$yrange + vec[2],
           zrange = X$domain$zrange + vec[3],
           marks = marks(X))
  attr(Y, "lastshift") <- vec
  return(Y)
}

#### sample ####
### sample.ppp ###
#' Sample a Planar Point Pattern
#'
#' @param X A \code{ppp}. The point pattern from which to sample.
#' @param size A numeric. The number of points to sample.
#' @return A \code{ppp}. The sampled point pattern.
#'
#' @family spatstat extensions
#' @seealso \code{\link{sample.pp3}}, \code{\link[base]{sample}}
#'
#' @export
sample.ppp <- function(X, size) {
  sam.n <- seq_len(spatstat::npoints(X))
  sam.pts <- sample(sam.n, size, replace = FALSE)
  sam.dat <- X[sam.pts]
  return(sam.dat)
}

### sample.pp3 ###
#' Sample A 3D Point Pattern
#'
#' @param X A \code{pp3}. The point pattern from which to sample.
#' @param size A numeric. The number of points to sample.
#' @return A \code{pp3}. The sampled point pattern.
#'
#' @family spatstat extensions
#' @seealso \code{\link{sample.ppp}}, \code{\link[base]{sample}}
#'
#' @export
sample.pp3 <- function(X, size) {
  sam.n <- seq_len(spatstat::npoints(X))
  sam.pts <- sample(sam.n, size, replace = FALSE)
  sam.dat <- X[sam.pts]
  return(sam.dat)
}

#### findClusters.pp3 ####
#' Find Clusters by NN Adjacency Marks
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
#' @name intensity.pp3-deprecated
#' @seealso \code{\link{rapt-deprecated}}
#' @keywords internal
NULL
#' @rdname rapt-deprecated
#' @section \code{intensity.pp3}:
#'   For \code{intensity.pp3}, use \code{\link[spatstat]{intensity.ppx}}
#'
#' @export
intensity.pp3 <- function(X, weights = NULL) {
  n <- spatstat::npoints(X)
  a <- spatstat::volume(domain(X))
  if (is.null(weights)) {
    if (spatstat::is.multitype(X)) {
      mks <- spatstat::marks(X)
      answer <- as.vector(table(mks))/a
      names(answer) <- levels(mks)
    }
    else answer <- n/a
    return(answer)
  } else
    warning('Weights is not yet implemented for pp3')
}

#### rownames.pp3 ####
#' Extends \code{\link[base:row+colnames]{rownames}} to
#' \code{\link[spatstat]{pp3}}.
#'
#' @param pat \code{pp3}. The point pattern from which to extract rownames.
#' @return Character vector. The rownames of the point pattern.
#' @seealso \code{\link[base:row+colnames]{rownames}}
rownames.pp3 <- function(pat) {
  dat <- rownames(as.data.frame(pat))
  return(dat)
}
#### rownames.pp3<- ####
#' Extends \code{\link[base:row+colnames]{rownames}} to
#' \code{\link[spatstat]{pp3}}.
#'
#' @param pat \code{pp3}. The point pattern to assign rownames.
#' @param lab character. The new rownames.
#' @seealso \code{\link[base:row+colnames]{rownames}}
"rownames.pp3<-" <- function(pat,lab) {
  rownames(as.data.frame(pat)) <- lab
  return(pat)
}

#### plot3d.pp3 ####
#' Plot a \code{\link[spatstat]{pp3}} in a manipulatable 3D plot.
#'
#' (requires the rgl library)
#' @param X \code{pp3}. The point pattern to visualize
#' @param ... Other arguments to pass to \code{plot3d} from the \code{rgl}
#' library.
#'
#' @family visualization functions
#' @seealso \code{\link[rgl]{plot3d}}
#'
#' @export
plot3d.pp3 <- function(X, ...) {
  rgl::plot3d(coords(X), ...)
}

#### quadratcount.pp3 ####
#' Count Points in Sub-Volumes of a 3D Point Pattern
#'
#' Divides volume into quadrats and counts the number of points in each quadrat.
#'
#' @param X The \code{\link[spatstat]{pp3}} object to split up.
#' @param nx,ny,nz Number of ractangular quadrats in the x, y, and z directions.
#'
#' @return A \code{data.frame} object containing the number of counts in each
#'   quadrat.
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat]{quadratcount}}
#'
#' @export
quadratcount.pp3 <- function(X, nx = 5, ny = 5, nz = 5){
  spatstat::verifyclass(X, "pp3")
  w <- domain(X)

  # create box3objects for each quadrat
  xlim <- w$xrange
  ylim <- w$yrange
  zlim <- w$zrange

  xbreaks <- seq(xlim[1], xlim[2],length.out = (nx+1))
  ybreaks <- seq(ylim[1], ylim[2],length.out = (ny+1))
  zbreaks <- seq(zlim[1], zlim[2],length.out = (nz+1))

  ntot <- nx*ny*nz
  gridvals <- list()
  cnt <- 1

  for(i in 1:nx){
    for(j in 1:ny){
      for(k in 1:nz){
        gridvals[[cnt]] <- box3(xrange = xbreaks[i:(i+1)],
                                yrange = ybreaks[j:(j+1)],
                                zrange = zbreaks[k:(k+1)])
        cnt <- cnt + 1
      }
    }
  }

  inside.tf <- lapply(gridvals, function(x){inside.boxx(X, w = x)})
  counts <- lapply(inside.tf, function(x){sum(x)})
  counts <- unlist(counts)
  return(data.frame(quad.no = seq(1,ntot), count = counts))
}

#### quadrats.pp3 ####
#' Divide a 3D Point Pattern into Sub-Volumes
#'
#' Divides volume into quadrats and returns them.
#'
#' @param X The \code{\link[spatstat]{pp3}} object to split up.
#' @param nx,ny,nz Number of rectangular quadrats in the x, y, and z directions,
#'   if you wish to split up your point patthern by number of boxes.
#' @param box.dims Vector containing the dimensions of the subsetted 3D boxxes,
#'   if you wish to define the individual boxx size. Use either \code{nx, ny,
#'   nz} or \code{box.dims}, but not both.
#'
#' @return A list containing the split up \code{pp3} objects.
#'
#' @family spatstat extensions
#' @seealso \code{\link[spatstat]{quadrats}}
#'
#' @export
quadrats.pp3 <- function(X, nx, ny, nz, box.dims = NULL){
  # Rewrite to use n = c(nx,ny,nz) with a default of c(1,1,1). Update with
  # description of how the function handles a box.dims that's not commensurate
  # with the domain size.
  spatstat::verifyclass(X, "pp3")
  w <- domain(X)
  xlim <- w$xrange
  ylim <- w$yrange
  zlim <- w$zrange

  # create box3objects for each quadrat
  if(is.null(box.dims)){
    xbreaks <- seq(xlim[1],xlim[2],length.out = (nx+1))
    ybreaks <- seq(ylim[1],ylim[2],length.out = (ny+1))
    zbreaks <- seq(zlim[1],zlim[2],length.out = (nz+1))

    ntot <- nx*ny*nz
  } else {
    xbreaks <- seq(xlim[1],xlim[2],by = box.dims[1])
    ybreaks <- seq(ylim[1],ylim[2],by = box.dims[2])
    zbreaks <- seq(zlim[1],zlim[2],by = box.dims[3])
  }


  gridvals <- list()
  cnt <- 1

  for(i in 1:(length(xbreaks)-1)){
    for(j in 1:(length(ybreaks)-1)){
      for(k in 1:(length(zbreaks)-1)){
        gridvals[[cnt]] <- box3(xrange = xbreaks[i:(i+1)],
                                yrange = ybreaks[j:(j+1)],
                                zrange = zbreaks[k:(k+1)])
        cnt <- cnt + 1
      }
    }
  }

  #browser()
  coo <- coords(X)
  if(!is.null(marks(X))){
    marks <- marks(X)
  }

  boxes <- lapply(gridvals, function(x){
    res.coo <- coo[inside.boxx(X, w = x),]
    if(!is.null(marks(X))) {
      res.marks <- marks[inside.boxx(X, w = x)]
      return(pp3(res.coo$x, res.coo$y, res.coo$z, x, marks = res.marks))
    } else {
      return(pp3(res.coo$x, res.coo$y, res.coo$z, x))
    }
  })

  return(boxes)
}

#### K3cross ####
# barely works... Needs corrections and inferface streamlining
K3multi <- function(X, I, J, r, breaks,
              correction = c("none", "isotropic", "translation"),
              ..., ratio = FALSE) {
  spatstat::verifyclass(X, "pp3")
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
          desc, fname = c("K", "list(I,J)"), yexp = quote(K[list(I,J)](r)))
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

#### studpermu.test.pp3 ####
#' Extends \code{\link[spatstat]{studpermu.test}} to pp3
#'
#' This function is not yet functional - many subfunctions must also be extended
#' to pp3.
#'
#' @param use.Tbar Logical value indicating choice of test statistic. If TRUE,
#'   use the alternative test statistic, which is appropriate for summary
#'   functions with roughly constant variance, such as K(r)/r or L(r). Defaults
#'   to FALSE
#' @param rinterval Numeric of length 2. Experimental
#'
#' @family spatstat extensions
#'
#' @seealso \code{\link[spatstat]{studpermu.test}}
# Add ability to supply a summary function directly...
studpermu.test.pp3 <- function (X, formula,
                                summaryfunction = K3est, ...,
                                nperm = 999, use.Tbar = FALSE, rinterval = NULL,
                                minpoints = 20, rmax = NULL, nrval = 128,
                                arguments.in.data = FALSE) {
  if (arguments.in.data & !is.hyperframe(X))
    stop(paste("X needs to be a hyperframe",
               "if arguments for summary function are to be retrieved"),
         call. = FALSE)
  stopifnot(is.function(summaryfunction))
  if (is.hyperframe(X)) {
    if (dim(X)[2] < 2)
      stop(paste("Hyperframe X needs to contain at least 2 columns,",
                 "one for patterns, one indicating groups"), call. = FALSE)
    data <- X
    Xclass <- unclass(X)$vclass
    factorcandidate <- Xclass %in% c("integer", "numeric",
                                     "character", "factor")
    ppcandidate <- Xclass == "pp3"
    Xnames <- names(X)
    names(factorcandidate) <- names(ppcandidate) <- names(Xclass) <- Xnames
    if (all(!factorcandidate) || all(!ppcandidate))
      stop(paste("Hyperframe X needs to contain at least a column",
                 "with point patterns, and one indicating groups"),
           call. = FALSE)
    if (!missing(formula)) {
      if (!inherits(formula, "formula"))
        stop(paste("Argument", dQuote("formula"), "should be a formula"))
      if (length(formula) < 3)
        stop(paste("Argument", sQuote("formula"), "must have a left hand side"))
      rhs <- spatstat.utils::rhs.of.formula(formula)
      ppname <- formula[[2]]
      if (!is.name(ppname))
        stop("Left hand side of formula should be a single name")
      ppname <- paste(ppname)
      if (!ppcandidate[ppname])
        stop(paste("Left hand side of formula",
                   "should be the name of a column of point patterns"),
             call. = FALSE)
      groupvars <- all.vars(as.expression(rhs))
      if (!all(groupvars %in% Xnames) || any(!factorcandidate[groupvars]))
        stop(paste("Not all variables on right hand side of formula",
                   "can be interpreted as factors"), call. = FALSE)
      group <- interaction(lapply(
        as.data.frame(data[,groupvars, drop = FALSE]), factor))
      newnames <- Xnames
      newnames[Xnames == ppname] <- "pp"
      names(data) <- newnames
      data$group <- group
    }
    else {
      thepp <- which.max(ppcandidate)
      thegroup <- which.max(factorcandidate)
      formula <- as.formula(paste(Xnames[thepp], "~", Xnames[thegroup]))
      newnames <- Xnames
      newnames[thepp] <- "pp"
      newnames[thegroup] <- "group"
      names(data) <- newnames
      data$group <- as.factor(data$group)
    }
  }
  else {
    if (!is.list(X))
      stop("X should be a hyperframe or a list of lists of point patterns")
    if (!is.list(X[[1]]) || !is.pp3(X[[1]][[1]]))
      stop("X is a list, but not a list of lists of point patterns")
    nams <- names(X)
    if (is.null(nams))
      nams <- paste("group", seq_along(X))
    pp <- list()
    group <- NULL
    for (i in seq_along(X)) {
      pp <- c(pp, X[[i]])
      group <- c(group, rep(nams[i], length(X[[i]])))
    }
    group <- as.factor(group)
    data <- hyperframe(pp = pp, group = group)
    ppname <- "pp"
  }
  framename <- deparse(substitute(X))
  fooname <- deparse(substitute(summaryfunction))
  OK <- sapply(data$pp, npoints) >= minpoints
  if ((nbad <- sum(!OK)) > 0)
    warning(paste(nbad, "patterns have been discarded",
                  "because they contained fewer than",
                  minpoints, "points"), call. = FALSE)
  data <- data[OK, , drop = FALSE]
  pp <- data$pp
  groupi <- as.integer(data$group)
  ngroups <- max(groupi)
  if (ngroups < 2)
    stop(paste("Sorry, after discarding patterns with fewer than",
               minpoints, "points,", if (ngroups < 1)
                 "nothing"
               else "only one group", "is left over.",
                "\n- nothing to compare, take a break!"),
         call. = FALSE)
  lev <- 1:ngroups
  m <- as.vector(table(groupi))
  if (any(m < 3))
    stop(paste("Data groups need to contain at least two patterns;",
               "\nafter discarding those with fewer than", minpoints,
               "points, the remaining group sizes are",
               spatstat.utils::commasep(m)),
         call. = FALSE)
  npossible <- factorial(sum(m))/prod(factorial(m))/prod(factorial(table(m)))
  if (is.nan(npossible))
    npossible <- Inf
  if (npossible < max(100, nperm))
    warning("Don't expect exact results - group sizes are too small")
  if (is.null(rmax)) {
    rmax <- sapply(pp, function(P) {
      spatstat::diameter(P$domain)/2
    })
    rmax <- min(rmax)
  }
  if (is.null(rinterval))
    rinterval <- c(0, rmax)
  ranger <- diff(range(rinterval))
  rr <- seq(0, rmax, length.out = nrval)
  taker <- rr >= rinterval[1] & rr <= rinterval[2]
  needcorx <- "correction" %in% names(formals(summaryfunction))
  gavecorx <- "correction" %in% names(list(...))
  corx <- if (needcorx && !gavecorx)
    "translation"
  else NULL
  fvlist <- if (arguments.in.data) {
    if (is.null(corx)) {
      multicall(summaryfunction, pp, data, rmax = rmax, ...)
    }
    else {
      multicall(summaryfunction, pp, data, rmax = rmax, ...,
                correction = corx)
    }
  }
  else {
    if (is.null(corx)) {
      with(data, summaryfunction(pp, rmax = rmax, ...))
    }
    else {
      with(data, summaryfunction(pp, rmax = rmax, ..., correction = corx))
    }
  }
  # can skip most of the previous if summary functions are supplied...
  # need to extract rmax, nrval from functions
  fvtemplate <- fvlist[[1]]
  valu <- attr(fvtemplate, "valu")
  argu <- attr(fvtemplate, "argu")
  foar <- sapply(lapply(fvlist, "[[", valu), "[", taker)
  combs <- combn(lev, 2)
  predigested <- list(lev = lev, foar = foar, m = m, combs = combs,
                      rrr = rr[taker], ranger = ranger)
  if (use.Tbar) {
    Tobs <- T.barstat(groupi, predigested)
    Tsim <- replicate(nperm, T.barstat(sample(groupi), predigested))
  }
  else {
    Tobs <- T.stat(groupi, predigested)
    Tsim <- replicate(nperm, T.stat(sample(groupi), predigested))
  }
  names(Tobs) <- if (use.Tbar)
    "Tbar"
  else "T"
  pval <- (1 + sum(Tobs < Tsim))/(1 + nperm)
  method <- c("Studentized permutation test for grouped point patterns",
              ifelse(is.hyperframe(X),
                     spatstat.utils::pasteFormula(formula), NULL),
              spatstat.utils::choptext(ngroups, "groups:",
                                       paste(levels(data$group),
                collapse = ", ")),
              spatstat.utils::choptext("summary function:",
                  paste0(fooname, ","), "evaluated on r in",
                  spatstat.utils::prange(rinterval)),
              spatstat.utils::choptext("test statistic:",
                                       ifelse(use.Tbar, "Tbar,", "T,"),
                       nperm, "random permutations"))
  fooshort <- switch(fooname, pcf = "pair correlation ",
                     Kinhom = "inhomogeneous K-",
                     Linhom = "inhomogeneous L-",
                     Kscaled = "locally scaled K-",
                     Lscaled = "locally scaled L-",
                     paste(substr(fooname,1, 1), "-", sep = ""))
  alternative <- c(paste("not the same ", fooshort, "function",
                         sep = ""))
  testerg <- list(statistic = Tobs, p.value = pval, alternative = alternative,
                  method = method, data.name = framename)
  class(testerg) <- c("studpermutest", "htest")
  fvs <- lapply(fvlist, "[.fv", j = c(argu, valu))
  fvs <- lapply(fvs, "attr<-", which = "alim", value = rinterval)
  testerg$curves <- spatstat::hyperframe(fvs = fvs, groups = data$group)
  fvtheo <- fvlist[[1]]
  spatstat::fvnames(fvtheo, ".y") <- "theo"
  attr(fvtheo, "alim") <- rinterval
  testerg$curvtheo <- fvtheo[, c(argu, "theo")]
  grmn <- lapply(lev, splitmean, ind = groupi, f = foar)
  testerg$groupmeans <- lapply(grmn, makefv, xvals = rr[taker],
                               template = fvtheo)
  return(testerg)
}

# These helpers are undocumented, but located at
# https://rdrr.io/github/spatstat/spatstat.core/src/R/studpermutest.R

splitmean <- function(l, ind, f) {
  apply(f[ , ind == l], 1, mean)
}
splitvarn <- function(l, ind, f, m) {
  apply(f[ , ind == l], 1, var) / m[l]
}
studentstat <- function(i, grmean, grvar) {
  (grmean[, i[1]] - grmean[, i[2]])^2 / (grvar[i[1],] + grvar[i[2], ])
}

T.stat <- function (ind = groupi, predigested) {
  # predigested should be a list with entries lev, foar, m, combs, rrr
  with(predigested, {
    grmean <- sapply(lev, splitmean, ind=ind, f=foar)
    grvar <- t(sapply(lev, splitvarn, ind=ind, f=foar, m=m))
    y <- apply(combs, 2, studentstat, grmean=grmean, grvar=grvar)
    sum(apply(y, 2, trapint, x = rrr))
  })
}

intstudent <- function(i, rrr, grmean, meangrvar) {
  trapint(rrr, (grmean[, i[1]] - grmean[, i[2]])^2 /
            (meangrvar[i[1]] + meangrvar[i[2]]))
}

T.barstat <- function (ind = groupi, predigested) {
  # predigested should be a list
  # with entries lev, foar, m, combs, rrr, ranger
  with(predigested, {
    grmean <- sapply(lev, splitmean, ind=ind, f=foar)
    grvar <- t(sapply(lev, splitvarn, ind=ind, f=foar, m=m))
    meangrvar <- apply(grvar, 1, trapint, x = rrr)/ranger
    sum(apply(combs, 2, intstudent,
              rrr=rrr, grmean=grmean, meangrvar=meangrvar))
    # trapint(rr[taker], grvar[i[1],] + grvar[i[2], ]))))
  })
}

makefv <- function(yvals, xvals, template) {
  fdf <- data.frame(r = xvals, y = yvals)
  argu <- fvnames(template, ".x")
  valu <- fvnames(template, ".y")
  names(fdf) <- c(argu,valu)
  fv(fdf, argu = argu, ylab = attr(template, "ylab"), valu = valu,
     fmla = attr(template,"fmla"), alim = attr(template, "alim"))
}

# Trapezoidal rule approximation to integral
trapint <- function(x, y) {
  nonan <- !is.na(y)
  nn <- sum(nonan)
  if(nn < 2L) return(0)
  Y <- y[nonan]
  X <- x[nonan]
  0.5 * sum( (Y[-1] + Y[-nn]) * diff(X))
}

# call foo(x, further arguments) repeatedly
# further arguments are taken from hyperframe H and ...
multicall <- function(foo, x, H, ...){
  stopifnot(is.hyperframe(H))
  if (is.hyperframe(x)) {
    x <- as.list(x)[[1]]
  } else if(!is.list(x))
    stop("in multicall: x should be a hyperframe or list", call.=FALSE)

  # check if same length
  nrows <- dim(H)[1]
  if (length(x) != nrows)
    stop(paste("in multicall: x and H need to have",
               "the same number of rows or list elements"),
         call.=FALSE)
  dotargs <- list(...)
  hnames <- names(H)
  argnames <- names(formals(foo))# always assume first argument is given

  ppname <- argnames[1]
  argnames <- argnames[-1]
  dotmatch <- pmatch(names(dotargs), argnames)
  dotmatched <- dotmatch[!is.na(dotmatch)]
  dotuseargs <- dotargs[!is.na(dotmatch)]
  restargs <- if(length(dotmatched) >0) argnames[-dotmatched] else argnames
  hmatch <- pmatch(hnames, restargs)
  huse <- !is.na(hmatch)
  lapply(seq_len(nrows), function (i)
    do.call(foo, c(list(x[[i]]),
                   as.list(H[i, huse, drop=TRUE, strip=FALSE]),
                   dotargs)))
}

#### Tstat.pp3 ####
#' Extends Tstat to pp3
#'
#' Tstat.pp3 extends the third-order summary statistic
#' \code{\link[spatstat]{Tstat}} to pp3
#'
#' @param correction One of "none" or "isotropic." "translation" correction is
#'   planned but not yet implemented.
#'
#' @family spatstat extensions
#'
#' @seealso \code{\link[spatstat]{Tstat}}
#' @references Schladitz, K. and Baddeley, A. (2000) A third order point process
#' characteristic. Scandinavian Journal of Statistics 27 (2000) 657-671.
Tstat.pp3 <- function (X, ..., rmax = NULL, nrval = 128,
                       correction = "border",
                       ratio = FALSE, verbose = TRUE) {
  spatstat::verifyclass(X, "pp3")
  npts <- spatstat::npoints(X)
  W <- spatstat::domain(X)
  areaW <- spatstat::volume(W)
  lambda <- npts/areaW
  lambda2 <- (npts * (npts - 1))/(areaW^2)
  lambda3 <- (npts * (npts - 1) * (npts - 2))/(areaW^3)
  if (is.null(rmax))
    rmax <- spatstat::diameter(W)/2
  r <- seq(0, rmax, length.out = nrval)
  breaks <- spatstat::breakpts.from.r(r)
  correction.given <- !missing(correction) && !is.null(correction)
  if (!correction.given)
    correction <- c("border", "bord.modif") #, "translate") not implemented yet
  correction <- spatstat::pickoption(
    "correction", correction,
    c(none = "none",
      border = "border", bord.modif = "bord.modif", trans = "translate",
      translate = "translate", translation = "translate", best = "best"),
    multi = TRUE)
  alim <- c(0, rmax)
  TT <- data.frame(r = r, theo = 5/12 * pi^2 * r^6)
  desc <- c("distance argument r", "theoretical Poisson %s")
  TT <- spatstat::fv(TT, "r",
                     quote(T(r)), "theo", , alim, c("r", "%s[pois](r)"),
                     desc, fname = "T") # blank is in ppp version...
  if (ratio) {
    denom <- lambda2 * areaW
    numT <- spatstat::eval.fv(denom * TT)
    denT <- spatstat::eval.fv(denom + TT * 0)
    attributes(numT) <- attributes(denT) <- attributes(TT)
    attr(numT, "desc")[2] <- "numerator for theoretical Poisson %s"
    attr(denT, "desc")[2] <- "denominator for theoretical Poisson %s"
  }
  close <- spatstat::closepairs(X, rmax, what = "ijd", twice = FALSE,
                                neat = FALSE)
  I <- close$i
  J <- close$j
  DIJ <- close$d
  nI <- length(I)
  continue <- TRUE
  if (verbose) {
    nTmax <- nI * (nI - 1)/2
    esttime <- exp(1.25 * log(nTmax) - 21.5)
    message(paste("Searching", nTmax, "potential triangles;",
                  "estimated time", spatstat::codetime(esttime)))
    if(esttime > 60)
      continue <- askYesNo("Estimated time greater than one minute. Continue?")
  }
  if (!continue) return(NULL)
  tri <- spatstat::trianglediameters(I, J, DIJ, nvert = npts)
  stopifnot(identical(colnames(tri), c("i", "j", "k", "diam")))
  II <- with(tri, c(i, j, k))
  DD <- with(tri, rep.int(diam, 3))
  if (any(correction == "none")) {
    wh <- whist(DD, breaks$val)
    numTun <- cumsum(wh)
    denTun <- lambda3 * areaW
    Tun <- numTun/denTun
    TT <- spatstat::bind.fv(TT, data.frame(un = Tun), "hat(%s)[un](r)",
                            "uncorrected estimate of %s", "un")
    if (ratio) {
      numT <- spatstat::bind.fv(numT, data.frame(un = numTun), "hat(%s)[un](r)",
                                "numerator of uncorrected estimate of %s", "un")
      denT <- spatstat::bind.fv(denT, data.frame(un = denTun), "hat(%s)[un](r)",
                                "denominator of uncorrected estimate of %s",
                                "un")
    }
  }
  if (any(correction == "border" | correction == "bord.modif")) {
    b <- bdist.points3(X)
    bI <- b[II]
    RS <- spatstat::Kount(DD, bI, b, breaks)
    if (any(correction == "bord.modif")) {
      denom.area <- eroded.volumes(W, r)
      numTbm <- RS$numerator
      denTbm <- lambda3 * denom.area
      Tbm <- numTbm/denTbm
      TT <- spatstat::bind.fv(TT, data.frame(bord.modif = Tbm),
                              "hat(%s)[bordm](r)",
                              "modified border-corrected estimate of %s",
                              "bord.modif")
      if (ratio) {
        numT <- spatstat::bind.fv(numT, data.frame(bord.modif = numTbm),
                        "hat(%s)[bordm](r)",
                        "numerator of modified border-corrected estimate of %s",
                        "bord.modif")
        denT <- spatstat::bind.fv(denT, data.frame(bord.modif = denTbm),
                                  "hat(%s)[bordm](r)",
                      "denominator of modified border-corrected estimate of %s",
                      "bord.modif")
      }
    }
    if (any(correction == "border")) {
      numTb <- RS$numerator
      denTb <- lambda2 * RS$denom.count
      Tb <- numTb/denTb
      TT <- spatstat::bind.fv(TT, data.frame(border = Tb), "hat(%s)[bord](r)",
                    "border-corrected estimate of %s", "border")
      if (ratio) {
        numT <- spatstat::bind.fv(numT, data.frame(border = numTb),
                                  "hat(%s)[bord](r)",
                                "numerator of border-corrected estimate of %s",
                                "border")
        denT <- spatstat::bind.fv(denT, data.frame(border = denTb),
                                  "hat(%s)[bord](r)",
                              "denominator of border-corrected estimate of %s",
                              "border")
      }
    }
  }
  if (any(correction == "translate")) {
    stop("Translation correction not yet implemented.")
    edgewt <- edgetri.Trans(X, tri[, 1:3])
    wh <- whist(tri$diam, breaks$val, edgewt)
    numTtrans <- 3 * cumsum(wh)
    denTtrans <- lambda3 * areaW
    Ttrans <- numTtrans/denTtrans
    h <- spatstat::diameter(W)/2
    Ttrans[r >= h] <- NA
    TT <- spatstat::bind.fv(TT, data.frame(trans = Ttrans), "hat(%s)[trans](r)",
                  "translation-corrected estimate of %s", "trans")
    if (ratio) {
      numT <- spatstat::bind.fv(numT, data.frame(trans = numTtrans),
                                "hat(%s)[trans](r)",
                            "numerator of translation-corrected estimate of %s",
                            "trans")
      denT <- spatstat::bind.fv(denT, data.frame(trans = denTtrans),
                                "hat(%s)[trans](r)",
                          "denominator of translation-corrected estimate of %s",
                          "trans")
    }
  }
  formula(TT) <- . ~ r
  spatstat::unitname(TT) <- spatstat::unitname(X)
  if (ratio) {
    formula(numT) <- formula(denT) <- . ~ r
    spatstat::unitname(denT) <- spatstat::unitname(TT)
    spatstat::unitname(numT) <- spatstat::unitname(TT)
    TT <- spatstat::rat(TT, numT, denT, check = FALSE)
  }
  return(TT)
}

#### bdist.points3 ####
#' Extension of \code{\link[spatstat]{bdist.points}}. Helper function for border correction \code{\link{bK3est}}.
#'
#' Finds the smallest distance to a boundary for each point in a point pattern.
#'
#' @param X The point pattern for analysis. A \code{\link[spatstat]{pp3}}
#'   object.
#' @return An object containing the shortest distance to the boundary for each
#'   point in the pattern X.
#' @seealso \code{\link[spatstat]{bdist.points}}
bdist.points3 <- function (X) {

  spatstat::verifyclass(X, "pp3")

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
  result <- pmin.int(x - xmin, xmax - x,
                     y - ymin, ymax - y ,
                     z - zmin , zmax - z)

  return(result)
}

#### bdist.points3.multi ####
#' Returns the shortest distances to boundaries in the x, y, and z directions separately.
#'
#' @param X The point pattern for analysis. A \code{\link[spatstat]{pp3}}
#'   object.
#' @return A data.frame containing the shortest distance to the closest three
#'   boundaries for each point in the pattern X.
#' @seealso \code{\link{bdist.points3}}
bdist.points3.multi <- function (X){

  spatstat::verifyclass(X, "pp3")

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

  result <- data.frame(x = pmin.int(x - xmin, xmax - x),
                       y = pmin.int(y - ymin, ymax - y),
                       z = pmin.int(z - zmin, zmax - z))

  return(result)
}
