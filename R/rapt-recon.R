#
# This file contains methods for reconstructing APT data
#

paraOrthoX <- function(x, y, a, d) {
  if(x == 0)
    x <- 1e-12
  x.dat <- (6^(2/3) * a^2 * (1 + 2*a*d) * x^2 * (x^2 + y^2) -
     6^(1/3) * (-9*a^4 * x^3 * (x^2 + y^2)^2 +
        sqrt(3) * sqrt(a^6 * x^6 * (x^2 + y^2)^3 *
            (2 + 12*a*d + 16*a^3 * d^3 + 3*a^2 *
               (8*d^2 + 9*(x^2 + y^2)))))^(2/3)) /
  (6*a^2 * (x^2 + y^2) * (-9*a^4 * x^3 * (x^2 + y^2)^2 +
    sqrt(3) * sqrt(a^6 * x^6 * (x^2 + y^2)^3 *
            (2 + 12*a*d + 16*a^3 * d^3 + 3*a^2 *
               (8 * d^2 + 9*(x^2 + y^2)))))^(1/3))
  return(x.dat)
}
paraOrthoY <- function(x, y, a, d) {
  if(x == 0)
    x <- 1e-12
  y.dat <- (y * (6^(2/3) * a^2 * (1 + 2*a*d) * x^2 * (x^2 + y^2) -
    6^(1/3) * (-9*a^4 * x^3 * (x^2 + y^2)^2 +
      sqrt(3) * sqrt(a^6 * x^6 * (x^2 + y^2)^3 *
                       (2 + 12*a*d + 16*a^3 * d^3 + 3*a^2 *
                          (8*d^2 + 9*(x^2 + y^2)))))^(2/3))) /
    (6*a^2 * x * (x^2 + y^2) * (-9*a^4 * x^3 * (x^2 + y^2)^2 +
      sqrt(3) * sqrt(a^6 * x^6 * (x^2 + y^2)^3 *
        (2 + 12*a*d + 16*a^3 * d^3 + 3*a^2 *
           (8*d^2 + 9*(x^2 + y^2)))))^(1/3))
  return(y.dat)
}
paraOrthoRec <- function(ato, a, d = 90e7, icf = 1, k = 0) {
  para.x <- apply(ato[,c('dx','dy')], 1, function(ato) {
    paraOrthoX(ato[1]*1e8, ato[2]*1e8, a = a, d = d)
  })
  para.y <- apply(ato[,c('dx','dy')], 1, function(ato) {
    paraOrthoY(ato[1]*1e8, ato[2]*1e8, a = a, d = d)
  })
  para.dat <- data.frame(x = para.x, y = para.y)
  para.dat <- apply(para.dat, 1, function(pos){
    r <- sqrt(pos[1]^2 + pos[2]^2)
    scale <- icf * (a*r)^k
    new <- scale * pos
    return(new)
  })
  para.dat <- as.data.frame(t(para.dat))
  return(para.dat)
}

paraRec <- function(ato, a, d = 90e7, icf = 1) {
  para.ind <- rownames(ato)
  para.dr <- sqrt(ato[,'dx']^2 + ato[,'dy']^2)
  para.dr <- 1e8 * para.dr
  para.theta <- atan2(ato[,'dy'], ato[,'dx'])
  para.fun <- function(r) {
    1/sqrt(2*a) * (r^2 + (d+4*a*d^2)/(4*a))^(1/4)
  }
  para.r <- para.fun(para.dr)
  para.r <- icf * para.r
  para.x <- para.r * cos(para.theta)
  para.y <- para.r * sin(para.theta)
  para.z <- -a * para.r^2
  para.dat <- data.frame(x = para.x, y = para.y, z = para.z)
  rownames(para.dat) <- para.ind
  return(para.dat)
}
