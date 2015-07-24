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
  attr(pos.dat, "meta") <- list(
    name = pos.name
  )
  return(pos.dat);
}
# Create a list of mass spec formulas from a table
readForm <- function(filepath) {
  data("isotopes");
  form.dat <- read.csv(filepath, row.names = 1, stringsAsFactors = F);
  form.chk <- which(check_chemform(isotopes, form.dat[, 1])[, "warning"]);
  if(length(form.chk) == 0) {
    form.dat <- t(form.dat);
    return(form.dat);
  }else {
    print(paste("The formulae:", form.dat[form.chk], "are not valid."));
  }
}
### Condition Data ###
# Create a pp3 object (from package "spatstat") from a pos data frame
createSpat <- function(pos) {
  pp3.box <- sapply(pos[-4], range);
  pp3.dat <- pp3(pos$x, pos$y, pos$z, pp3.box);
  return(pp3.dat);
}
# Create a MassSpectrum object (from package "MALDIquant") from a pos data frame
createSpec <- function(pos, res = 0.001) {
  ms.range <- range(pos["mass"]);
  ms.breaks <- seq(ms.range[1], ms.range[2], res)
  ms.hist <- hist(pos, ms.breaks, plot = F);
  ms.dat <- createMassSpectrum(
    ms.hist$mids, ms.hist$counts,
    metaData = attr(pos,"meta")
  );
  return(ms.dat);
}
### Write Data ###
# Export methods?
