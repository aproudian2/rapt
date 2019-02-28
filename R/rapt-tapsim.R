library(Rcpp11)
sourceCpp('tapsim_reader.cpp')

#### readTAPSim ####
#' Read a TAPSim results file
#'
#' Given a filepath, read in a TAPSim file. Currently only implemented for
#' results.
#'
#' @param fp The filepath
#' @param type The file type to be read. Currently, only "results" is
#' implemented.
#' @return A data.frame
#' @seealso \url{
#'   http://www.uni-stuttgart.de/imw/mp/forschung/atom_probe_RD_center/howto.pdf
#' }
#' @export
readTAPSim <- function(fp, type = 'results') {
  df <- read_TAPSim(fp, type)
  return(df)
}
