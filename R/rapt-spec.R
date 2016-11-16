#
# This file contains methods for working with the mass spectrum of APT data
#

# Define ranges in MassSpectrum
rangeSpec <- function(form, charge = 1:4, resolution = 400, thresh = 20) {
  if(is(form, "cform")) {
    cform <- as.data.frame(t(form));
    form <- form["formula", ];
  } else {
    warning("The given formulae are not of class cform")
  }
  rng.pat <- lapply(charge, function(x) {
    enviPat::isopattern(
      rapt::isotopes, chemforms = form,
      threshold = 0.01, charge = x,
      verbose = F
    )
  })
  rng.env <- lapply(
    rng.pat, enviPat::envelope,
    frac = 0.1, env = "CauchyLorentz", resolution = resolution,
    verbose = F
  )
  rng.raw <- lapply(rng.env, function(x) {
    lapply(x, function(y){
      range(subset(y, select = "m/z", subset = "abundance" > thresh))
    })
  })
  rng.dat <- lapply(charge, function(x) {
    lapply(seq_along(rng.raw[[x]]), function(y) {
      data.frame(
        formula = names(rng.raw[[x]])[y],
        charge = x,
        start = rng.raw[[x]][[y]][1],
        end = rng.raw[[x]][[y]][2]
      )
    })
  })
  rng.dat <- do.call(rbind, unlist(rng.dat, recursive = F))
  rng.dat <- rng.dat[order(rng.dat$formula, rng.dat$charge),]
  rownames(rng.dat) <- apply(rng.dat, 1, function(x){
    paste(c(x["formula"], x["charge"]), collapse = "+")
  })
  if (exists("cform")) {
    rng.dat$color <- as.character(
      cform$color[match(rng.dat$formula, cform$formula)])
  } else {
    rng.dat$color <- rep_len(NA, dim(rng.dat)[1])
  }
  return(rng.dat);
}
# Merge RNG objects
mergeRange <- function(rng1, rng2) {

}
# Remove ranges from an RNG
removeRange <- function(rng, rem) {
  rng.sub <- setdiff(rownames(rng), rem);
  rng.new <- rng[rng.sub, ];
  return(rng.new);
}
# Create a reference MassPeaks object for warping
refSpec <- function(form, charge = 1:4, resolution = 400) {
  if(is(form, "cform")) {
    cform <- as.data.frame(t(form));
    form <- form["formula", ];
  }
  ref.pat <- lapply(charge, function(x){
    enviPat::isopattern(
      rapt::isotopes, chemforms = form,
      threshold = 0.01, charge = x,
      verbose = F
    )
  })
  ref.env <- lapply(
    ref.pat, enviPat::envelope,
    frac = 0.1, env = "CauchyLorentz", resolution = resolution,
    verbose = F
  )
  ref.vdet <- lapply(
    ref.env, enviPat::vdetect,
    detect = "intensoid",
    plotit = F, verbose = F
  )
  ref.vdet <- do.call(rbind,unlist(ref.vdet, recursive = F))
  ref.dat <- createMassPeaks(ref.vdet[, "m/z"], ref.vdet[, "abundance"])
  return(ref.dat)
}
# Auto define ranges from a MassPeaks object
rangePeaks <- function(spec, peaks, cut = 0.1) {
  if(!isMassSpectrum(spec)) {
    stop('spec is not a MassSpectrum object')
  }
  if(!isMassPeaks(peaks)) {
    stop('peaks is not a MassPeaks object')
  }
  rng.ind <- match(peaks@mass, spec@mass)
  rng.len <- length(spec@intensity)
  rng.bounds <- sapply(rng.ind, function(x) {
    max.ind <- which(spec@intensity[-seq(0, x)] - cut*spec@intensity[x] < 0,
                     arr.ind = T)[1]
    max.val <- spec@mass[-seq(0, x)][max.ind]
    min.ind <- rev(
      which(spec@intensity[-seq(x, rng.len)] - cut*spec@intensity[x] < 0,
            arr.ind = T))[1]
    min.val <- spec@mass[-seq(x, rng.len)][min.ind]
    c(min.val,max.val)
  })
  rng.dat <- data.frame(formula = rep_len(NA, length(peaks)),
                        charge = rep_len(NA, length(peaks)),
                        start = rng.bounds[1,],
                        end = rng.bounds[2,],
                        color = rep_len(NA, length(peaks))
                        )
  return(rng.dat)
}
