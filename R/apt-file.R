# Read a POS file (from IVAS reconstruction) into a data frame
readPOS <- function(filepath) {
   pos.len <- file.info(filepath)['size'] / 4;
   pos.len <- as.numeric(pos.len);
   pos.raw <- readBin(filepath, what = "numeric",
                      size = 4, endian = "big", n = pos.len);
   pos.mat <- matrix(pos.raw, ncol = 4, byrow = T);
   pos.dat <- as.data.frame(pos.mat);
   names(pos.dat) <- c("x", "y", "z", "mass");
   return(pos.dat);
}
# Create a pp3 object (from package "spatstat") from a pos data frame
createPP3 <- function(pos) {
   require("spatstat");
   pp3.box <- sapply(pos[-4], range);
   pp3.dat <- pp3(pos$x, pos$y, pos$z, pp3.box);
   return(pp3.dat)
}
# Create a list of mass spec formulas from a table
readForm <- function(filepath) {
   require("enviPat");
   data("isotopes");
   form.dat <- read.csv(filepath, row.names = 1, stringsAsFactors = F);
   form.chk <- which(check_chemform(isotopes, form.dat[,1])[,"warning"]);
   if(length(form.chk) == 0) {
      form.dat <- t(form.dat);
      return(form.dat);
   }else {
      print(paste("The formulae:",form.dat[form.chk],"are not valid."))
   }
}
# Shade regions under a mass spec
# From Gavin Simpson
# www.fromthebottomoftheheap.net/2013/01/11/shading-regions-under-a-curve
polyCurve <- function(x, y, from, to, n = 50, miny,
                      col = "red", border = col) {
   drawPoly <- function(fun, from, to, n = 50, miny, col, border) {
      Sq <- seq(from = from, to = to, length = n)
      polygon(x = c(Sq[1], Sq, Sq[n]),
              y = c(miny, fun(Sq), miny),
              col = col, border = border)
   }
   lf <- length(from)
   stopifnot(identical(lf, length(to)))
   if(length(col) != lf)
      col <- rep(col, length.out = lf)
   if(length(border) != lf)
      border <- rep(border, length.out = lf)
   if(missing(miny))
      miny <- min(y)
   interp <- approxfun(x = x, y = y)
   mapply(drawPoly, from = from, to = to, col = col, border = border,
          MoreArgs = list(fun = interp, n = n, miny = miny))
   invisible()
}
