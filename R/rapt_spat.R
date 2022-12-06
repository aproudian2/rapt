#
# This file contains spatial processing functions for APT data
#

#### cylinder ####
cylinder <- function(x, ...) {
  UseMethod("cylinder", x)
}
cylinder.pp3 <- function(pat, ind, r = 1) {
  cylX <- pat$data$x[ind]
  cylY <- pat$data$y[ind]
  dat <- subset(pat, subset = r > sqrt((x - cyl)^2 + (y - cylY)^2))
  return(dat)
}
#### cone ####
cone <- function(x, ...) {
  UseMethod("cone", x)
}
cone.pp3 <- function(pat, r = 1, phi = pi / 4) {
  coneX <- as.numeric(pat$data$x[1])
  coneY <- as.numeric(pat$data$y[1])
  coneZ <- as.numeric(pat$data$z[1])
  a <- tan(phi)
  k <- a * r
  dat <- subset(pat,
    subset =
      z >= a * sqrt((x - coneX)^2 + (y - coneY)^2) + coneZ - k |
        z <= -a * sqrt((x - coneX)^2 + (y - coneY)^2) + coneZ + k
  )
  return(dat)
}
cone.ppp <- function(pat, r = 1, phi = pi / 4) {
  coneX <- pat$x[1]
  coneY <- pat$y[1]
  a <- tan(phi)
  k <- a * r
  dat <- subset(pat,
    subset =
      y > a * sqrt((x - coneX)^2) + coneY - k |
        y < -a * sqrt((x - coneX)^2) + coneY + k
  )
  return(dat)
}
#### pointDist ####
pointDist <- function(pat) {
  orig <- as.data.frame(pat[1])
  pat <- as.data.frame(pat)
  dat <- apply(pat, 1, function(point) {
    sqrt(sum((point - orig)^2))
  })
  dat <- as.numeric(dat)
  dat <- dat[dat > 0]
  return(dat)
}
#### coneRDF ####
coneRDF <- function(x, ...) {
  UseMethod("coneRDF", x)
}
coneRDF.ppp <- function(pat, cl, r = 1, phi = pi / 4, nrval = 128) {

}
coneRDF.pp3 <- function(pat, cl, r = 1, phi = pi / 4, nrval = 128) {
  rdf.n <- npoints(pat)
  rdf.ind <- 1:rdf.n
  rdf.dist <- parLapply(cl, rdf.ind, function(ind) {
    subpat <- pat[ind:rdf.n]
    pointDist(cone(subpat, r = r, phi = phi))
  })
  rdf.dist <- unlist(rdf.dist)
  rdf.breaks <- seq(0, max(rdf.dist), length.out = nrval + 1)
  rdf.hist <- hist(rdf.dist, breaks = rdf.breaks, plot = F)
  rdf.r <- rdf.hist$mids
  rdf.dr <- diff(rdf.breaks)[1]
  rdf.a <- tan(phi)
  rdf.density <- intensity.ppx(pat) # Improve accuracy
  rdf.dvol <- 4 * pi * (1 - 1 / sqrt(1 + (rdf.a + rdf.a^2)^2)) * rdf.dr
  rdf.g <- rdf.hist$counts / (rdf.density * rdf.dvol * rdf.r^2)
  rdf.dat <- data.frame(r = rdf.r, g = rdf.g)
  return(rdf.dat)
}
