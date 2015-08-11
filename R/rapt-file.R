#
# This file contains methods for importing and conditioning APT data
#

### Load Data ###
# Read a POS file (from IVAS reconstruction) into a data frame
readPOS <- function(filepath) {
  pos.len <- file.info(filepath)['size'] / 4;
  pos.len <- as.numeric(pos.len);
  pos.raw <- readBin(
    filepath, what = "numeric",
    size = 4, endian = "big", n = pos.len
  );
  pos.mat <- matrix(pos.raw, ncol = 4, byrow = T);
  pos.dat <- as.data.frame(pos.mat);
  names(pos.dat) <- c("x", "y", "z", "mass");
  pos.name <- sub(".pos", "", basename(filepath));
  attr(pos.dat, "metaData") <- list(
    name = pos.name
  )
  return(pos.dat);
}
# Create a list of mass spec formulas from a table (class "cform")
readForm <- function(filepath) {
  data("isotopes");
  form.dat <- read.csv(filepath, stringsAsFactors = F);
  form.chk <- which(check_chemform(isotopes, form.dat[, "formula"])[, "warning"]);
  if(length(form.chk) == 0) {
    form.dat <- t(form.dat);
    class(form.dat) <- c("matrix", "cform")
    return(form.dat);
  }else {
    print(paste("The formulae:", form.dat[form.chk], "are not valid."));
  }
}
### Condition Data ###
# Create a pp3 object (from package "spatstat") from a pos data frame
createSpat <- function(pos, win = NULL) {
  pp3.box <- win;
  if(is.null(win)) {
    pp3.box <- sapply(pos[-4], range);
  }
  pp3.dat <- pp3(pos$x, pos$y, pos$z, pp3.box);
  attr(pp3.dat, "metaData") <- attr(pos, "metaData");
  return(pp3.dat);
}
# Create a MassSpectrum object (from package "MALDIquant") from a pos data frame
createSpec <- function(pos, res = 0.001) {
  ms.range <- range(pos["mass"]);
  ms.range <- ms.range + c(-res, res);
  ms.breaks <- seq(ms.range[1], ms.range[2], res)
  ms.hist <- hist(pos[,"mass"], ms.breaks, plot = F);
  ms.dat <- createMassSpectrum(
    ms.hist$mids, ms.hist$counts,
    metaData = attr(pos, "metaData")
  );
  return(ms.dat);
}
### Write Data ###
# Export methods?
