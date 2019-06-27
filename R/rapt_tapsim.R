#### readResults ####
#' Read a TAPSim results file
#'
#' Given a filepath, read in a TAPSim results file. Currently only implemented
#' for ASCII results files.
#'
#' @param fp The filepath
#' @param type The file type to be read. Currently, only ASCII is implemented.
#' @return A data.frame
#' @seealso \url{
#'   http://www.uni-stuttgart.de/imw/mp/forschung/atom_probe_RD_center/howto.pdf
#' }
#' @seealso \url{https://d-nb.info/1138282715/34}
#' @export
readResult <- function(fp, type = 'ASCII') {
    txt <- readLines(fp, n = 100)  # Not ideal to hard code the max...
    skip <- grep(type, txt)
    dat <- read.delim(fp, header = FALSE, skip = skip)
    names(dat) <- c('index','id','number','voltage','startX','startY','startZ',
                    'stopX','stopY','stopZ','tof','probability',
                    'potentialBefore',
                    'fieldBeforeX','fieldBeforeY','fieldBeforeZ',
                    'potentialAfter','fieldAfterX','fieldAfterY','fieldAfterZ',
                    'normalX','normalY','normalZ','apexX','apexY','apexZ')
    return(dat)
}
