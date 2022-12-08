# I removed this so that rjitter.pp3 actually workds like rjitter.ppp
rjitter.pp3 <- function(X, domain = box3()) {
  spatstat.geom::verifyclass(X, "pp3")
  nX <- spatstat.geom::npoints(X)
  if (nX == 0) {
    return(X)
  }
  W <- X$domain
  D <- spatstat.random::runifpoint3(nX, domain = domain)
  xnew <- X$data$x + D$data$x
  ynew <- X$data$y + D$data$y
  znew <- X$data$z + D$data$z
  new <- pp3(xnew, ynew, znew, W)
  ok <- subset(new,
               subset =
                 (x > W$xrange[1] & x < W$xrange[2]) &
                 (y > W$yrange[1] & y < W$yrange[2]) &
                 (z > W$zrange[1] & z < W$zrange[2])
  )
  return(ok)
}

