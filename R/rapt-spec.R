# Define ranges in MassSpectrum
rangeSpec <- function(form, charge = 1:4, resolution = 600) {
  rng.pat <- lapply(charge, function(x){
    enviPat::isopattern(
      isotopes, chemforms = form,
      threshold = 0.01, charge = x,
      verbose = F
    )
  })
  rng.env <- lapply(
    rng.pat, enviPat::envelope,
    frac = 0.1, env = "CauchyLorentz", resolution = resolution,
    verbose = F
  )
  rng.dat <- lapply(rng.env, function(x){
    lapply(x, function(y){
      range(subset(y, select = "m/z", subset = "abundance" > 20))
    })
  })
  return(rng.dat)
}
# Create a reference MassPeaks object for warping
refSpec <- function(form, charge = 1:4, resolution = 600) {
  ref.pat <- lapply(charge, function(x){
    enviPat::isopattern(
      isotopes, chemforms = form,
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
  ref.dat <- createMassPeaks(ref.vdet[,"m/z"],ref.vdet[,"abundance"])
  return(ref.dat)
}
