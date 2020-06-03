#
# This file contains functions pertaining to TAPSim
#

#### Read Data ####

### readResult ###
#' Read TAPSim Results
#'
#' Given a filepath, read in a TAPSim results file.
#' \emph{Currently only implemented for ASCII results files.}
#'
#' @param fp Character. The filepath
#' @param type Character. The file type to be read. Currently, only ASCII is
#'   implemented.
#' @return A data.frame with the structure of the TAPSim output.
#'
#' @family TAPSim functions
#' @seealso \url{
#'   http://www.uni-stuttgart.de/imw/mp/forschung/atom_probe_RD_center/howto.pdf
#' }
#' @seealso \url{https://d-nb.info/1138282715/34}
#'
#' @export
readResult <- function(fp, type = 'ASCII') {
    txt <- readLines(fp, n = 100)  # Not ideal to hard code the max...
    skip <- grep(type, txt)
    dat <- read.delim(fp, header = FALSE, skip = skip)
    n <- ncol(dat)
    if(n == 26) {
        names(dat) <- c("index","id","number","voltage",
                        "startX","startY","startZ",
                        "stopX","stopY","stopZ","tof","probability",
                        "potentialBefore",
                        "fieldBeforeX","fieldBeforeY","fieldBeforeZ",
                        "potentialAfter",
                        "fieldAfterX","fieldAfterY","fieldAfterZ",
                        "normalX","normalY","normalZ","apexX","apexY","apexZ")
    } else if(n == 18) {
        names(dat) <- c("index","id","number","voltage",
                        "startX","startY","startZ",
                        "stopX","stopY","stopZ","tof","probability",
                        "normalX","normalY","normalZ","apexX","apexY","apexZ")
    } else {
        warning("Unknown file structure. Returning unnamed columns.")
    }
    return(dat)
}

#### Condition Data ####

### resultToPOS ###
#' Convert TAPSim Results to POS
#'
#' @param res data.frame. The results \code{data.frame} as returned by
#'   \code{\link{readResult}}
#' @param clip.radius Numeric. The detector radius at which to clip the points
#' @return A \code{data.frame} containing the (x,y,z) position and "mass" (id)
#'   of the point in the stame structure of a POS as created by
#'   \code{\link{readPOS}}.
#'
#' @family TAPSim functions
#' @seealso \code{\link{readPOS}}
#'
#' @export
resultToPOS <- function(res, clip.radius = NULL) {
    pos <- data.frame(x = res$startX*1e9,
                      y = res$startY*1e9,
                      z = res$startZ*1e9,
                      mass = res$id)
    if (!is.null(clip.radius)) {
        r <- sqrt(res$stopX^2 + res$stopY^2)
        pos <- pos[r <= clip.radius,]
    }
    return(pos)
}

### resultToDet ###
#' Convert TAPSim Result to Detector \code{\link[spatstat]{ppp}}
#'
#' @param res data.frame. The results \code{data.frame} as returned by
#'   \code{\link{readResult}}
#' @param clip.radius Numeric. The detector radius at which to clip the points.
#' @return A \code{\link[spatstat]{ppp}} containing the detector positions
#'   marked by their id.
#'
#' @family TAPSim functions
#' @seealso \code{\link{createDet}}
#'
#' @export
resultToDet <- function(res, clip.radius = NULL) {
    det.df <- data.frame(dx = res$stopX, dy = res$stopY, mass = res$id)
    det <- createDet(det.df)
    marks(det) = det.df$mass
    if (!is.null(clip.radius)) {
        Window(det) <- disc(radius = clip.radius)
    }
    return(det)
}

#### Write Data ####

### writeNode ###
#' Write a Node File
#'
#' \code{writeNode} writes a node file to the TAPSim format for creating a mesh
#' file for TAPSim.
#'
#' @param node data.frame A \code{data.frame} with columns correponding to the
#' x, y, z, and id of each point in order
#' @param fp Character. The file to which to write the data.
#' @return None.
#'
#' @seealso \url{
#'   http://www.uni-stuttgart.de/imw/mp/forschung/atom_probe_RD_center/howto.pdf
#' }
#' @seealso \url{https://d-nb.info/1138282715/34}
#'
#' @export
writeNode <- function(node, fp) {
    cat('ASCII', nrow(node), 0, 0, fill = TRUE,
        file = fp)
    write.table(format(node, digits = 6), file = fp,
                quote = FALSE, sep = '\t', eol = '\n', append = TRUE,
                row.names = FALSE, col.names = FALSE)
}

