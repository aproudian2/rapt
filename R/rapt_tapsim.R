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
    n <- ncol(dat)
    if(n == 26) {
        names(dat) <- c('index','id','number','voltage',
                        'startX','startY','startZ',
                        'stopX','stopY','stopZ','tof','probability',
                        'potentialBefore',
                        'fieldBeforeX','fieldBeforeY','fieldBeforeZ',
                        'potentialAfter',
                        'fieldAfterX','fieldAfterY','fieldAfterZ',
                        'normalX','normalY','normalZ','apexX','apexY','apexZ')
    } else if(n == 18) {
        names(dat) <- c('index','id','number','voltage',
                        'startX','startY','startZ',
                        'stopX','stopY','stopZ','tof','probability',
                        'normalX','normalY','normalZ','apexX','apexY','apexZ')
    } else {}
    return(dat)
}

#### resultToPOS ####
#' Convert a TAPSim results file to a POS
#'
#' @param res The results data.frame
#' @param clip.radius The detector radius at which to clip the points
#' @return A data.frame containing the (x,y,z) position and "mass' (id) of the
#'   point in the stame structure of a POS as created by \code{\link{readPOS}}.
#' @seealso \code{\link{readPOS}}, \code{\link{readResult}}
#' @export
resultToPOS <- function(res, clip.radius = NULL) {
    pos <- with(res, data.frame(x = startX, y = startY, z = startZ, mass = id))
    if (!is.null(clip.radius)) {
        r <- with(res, sqrt(stopX^2 + stopY^2))
        pos <- pos[r <= clip.radius,]
    }
    return(pos)
}

#### resultToDet ####
#' Convert a TAPSim results file to a detector point pattern
#'
#' @param res The results data.frame
#' @param clip.radius The detector radius at which to clip the points
#' @return A ppp containing the detector positions marked by their id
#' @seealso \code{\link{createDet}}, \code{\link{readResult}}
#' @export
resultToDet <- function(res, clip.radius = NULL) {
    det.df <- with(res, data.frame(dx = stopX, dy = stopY, mass = id))
    det <- createDet(det.df)
    marks(det) = det.df$mass
    if (!is.null(clip.radius)) {
        Window(det) <- disc(radius = clip.radius)
    }
    return(det)
}
