# Define ranges in MassSpectrum
rangeSpec <- function(form, resolution = 600) {
  rng.pat <- lapply(1:4, function(x){
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
