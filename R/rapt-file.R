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
# Read an ATO file (from IVAS reconstruction) into a data frame
readATO <- function(filepath) {
  ato.len <- file.info(filepath)['size'] / 4;
  ato.len <- as.numeric(ato.len);
  ato.file <- file(filepath, open = "rb");
  seek(ato.file, where = 8);
  ato.raw <- readBin(
    ato.file, what = "numeric",
    size = 4, endian = "little", n = ato.len
  );
  ato.mat <- matrix(ato.raw, ncol = 14, byrow = T);
  ato.dat <- as.data.frame(ato.mat);
  names(ato.dat) <- c("x", "y", "z", "mass", "clusID", "pIndex", "Vdc", "TOF", "dx", "dy", "Vp", "shank", "FouR", "FouI");
  ato.name <- sub(".ato", "", basename(filepath));
  attr(ato.dat, "metaData") <- list(
    name = ato.name
  )
  close(ato.file);
  return(ato.dat);
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
    pp3.box <- sapply(pos[1:3], range);
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
