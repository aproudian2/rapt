#' Feature Extraction
#' @description function for extracting features for ML algorithm
#'
#' @param i an integer.  This will be the iteration number.  It determines the seed number
#' in the \code{\link[rapt]{clustersim}} function.
#' @param pattern a 3D point pattern (object of type \code{pp3}) to use as the underlying point pattern (UPP) onto
#' which generate clusters
#' @param rrl output of  \code{\link[rapt]{average_relabelings}}.  Used as expected values for summary functions
#' @param pcp a numeric.  Guest concentration
#' @param tol a numeric. 'tol' or tolerance in \code{\link[rapt]{clustersim}}
#' @param maxKr a numeric.  Maximum value to calculate K to in \code{\link[spatstat]{K3est}}
#' @param nKr a numeric.  Number of intervals at which to calculate K in  \code{\link[spatstat]{K3est}}
#' @param maxGr a numeric.  Maximum value to calculate G to in \code{\link[spatstat]{G3est}}
#' @param nGr a numeric.  Number of intervals at which to calculate G in  \code{\link[spatstat]{G3est}}
#' @param maxGXGHr a numeric.  Maximum value to calculate GXGH to in \code{\link[rapt]{G3cross}}
#' @param maxGXHGr a numeric.  Maximum value to calculate GXHG to in \code{\link[rapt]{G3cross}}
#' @param nGXr a numeric.  Number of intervals at which to calculate G in  \code{\link[rapt]{G3cross}}
#' @param vside a numeric. voxel side length for \code{\link[spatstat]{F3est}}
#' @param jitter a numeric.  Scale of perturbations. 'radius' value in  \code{\link[rapt]{rjitter.pp3}}
#' @param feats a vector.  Names of all features to be kept
#' @param total_time_stat a numeric. Leave as default to let it algorithm calculate iteration time
#'
#' @details Take an point pattern, simulate clustering, calculate summary functions on this new point
#' pattern, extract features from summary functions.
#' @export
create_training_data <- function(i, pattern, rrl, nper = 1, params,
                                 pcp = 0.5, tol = 0.005,
                                 maxKr = 10, nKr = 200,
                                 maxGr = 5, nGr = 1000,
                                 maxGXGHr = 3, maxGXHGr = 8, nGXr = 1000,
                                 vside = 0.3,
                                 jitter = 0,
                                 feats = c(
                                   "G_max_diff", "G_max_diff_r", "G_min_diff",
                                   "G_zero_diff_r", "F_min_diff", "F_min_diff_F",
                                   "Tm", "Rm", "Rdm", "Rddm", "Tdm",
                                   "GXGH_min_diff", "GXGH_95diff_r", "GXGH_FWHM"
                                 ),
                                 total_time_start = as.numeric(Sys.time())) {
  # total_time_start = Sys.time()
  cycle_time_start <- as.numeric(Sys.time())
  set.seed(i)
  print(paste("Interation:", i, ",  PCP:", pcp, ",  tol:", tol))
  cluster <- clustersim(pattern, pattern, 0.5,
    pcp = pcp,
    cr = params$cr,
    rho1 = params$rhoc,
    rho2 = params$rhob,
    rb = params$rb,
    pb = params$pb,
    tol = tol,
    s = i
  )
  names(jitter) <- "jitter"
  if (is.numeric(cluster)) {
    out <- unlist(c(params, jitter, rep(NA, length(feats) + 5)))
  } else {
    if (jitter) {
      pattern <- rjitter(pattern, radius = jitter)
    }
    rm(pattern)

    # marks must be factor for G3cross to work
    # convert all B to A in cluster[[2]]
    cluster[[2]]$data$marks[cluster[[2]]$data$marks == "B"] <- "A"
    marks(cluster[[2]]) <- as.factor(marks(cluster[[2]]))
    # times = c()
    features <- c()

    # Extract Gmet
    G_feats_all <- c("G_max_diff", "G_max_diff_r", "G_min_diff", "G_zero_diff_r")
    G_feats_to_use <- G_feats_all[G_feats_all %in% feats]
    G_start <- as.numeric(Sys.time())
    if (length(G_feats_to_use)) {
      G_result <- G3est(cluster[[1]], rmax = maxGr, nrval = nGr, correction = "km")
      G_features <- g3features(rrl$rrl_G$r, G_result$km, rrl$rrl_G$mmean)
      features <- c(features, G_features)
      # G_stop <- as.numeric(Sys.time())
      # names(G_time) = "G_time"
      # times = c(times, G_time)
    }
    G_stop <- as.numeric(Sys.time())
    G_time <- G_stop - G_start


    # Extract Fmet
    F_feats_all <- c("F_min_diff", "F_min_diff_F")
    F_feats_to_use <- F_feats_all[F_feats_all %in% feats]
    F_start <- as.numeric(Sys.time())
    if (length(F_feats_to_use)) {
      F_result <- F3est(cluster[[1]], rmax = maxGr, nrval = nGr, correction = "km")
      F_features <- f3features(rrl$rrl_F$r, F_result$km, rrl$rrl_F$mmean)
      features <- c(features, F_features)
    }

    F_stop <- as.numeric(Sys.time())
    F_time <- F_stop - F_start


    # Measure T(r) on the pattern
    ##
    K_feats_all <- c("Tm", "Rm", "Rdm", "Rddm", "Tdm")
    K_feats_to_use <- K_feats_all[K_feats_all %in% feats]
    K_start <- as.numeric(Sys.time())
    if (length(K_feats_to_use)) {
      # K_start <- as.numeric(Sys.time())
      K_result <- K3est(cluster[[1]], rmax = maxKr, nrval = nKr, correction = "trans")
      T3 <- sqrt(K_result$trans) - sqrt(rrl$rrl_K$mmean)
      K_features <- k3features(K_result$r, T3, toplot = TRUE, xlim = c(0, maxKr))
      features <- c(features, K_features)
      # K_stop <- as.numeric(Sys.time())
      # K_time <- K_stop - K_start
      # names(K_time) = "K_time"
      # times = c(times, K_time)
    }
    K_stop <- as.numeric(Sys.time())
    K_time <- K_stop - K_start


    # Extract Gcross_met for  A to C
    GXGH_feats_all <- c("GXGH_min", "GXGH_95diff_r", "GXGH_FWHM", "GXGH_absmax_diff", "GXGH_max_diff")
    GXGH_feats_to_use <- GXGH_feats_all[GXGH_feats_all %in% feats]
    GXGH_start <- as.numeric(Sys.time())
    if (length(GXGH_feats_to_use)) {
      GXGH_result <- G3cross(cluster[[2]],
        i = dopant_formula, j = host_formula,
        correction = "km", rmax = maxGXGHr, nrval = nGXr
      )
      GXGH_features <- g3Xfeatures(GXGH_result$r,
        gvals_new = GXGH_result$km,
        gvals_old = rrl$rrl_GXGH$mmean, type = "GH"
      )
      features <- c(features, GXGH_features)
    }
    GXGH_stop <- as.numeric(Sys.time())
    GXGH_time <- GXGH_stop - GXGH_start

    # Extract Gcross_met for  A to C
    GXHG_feats_all <- c("GXHG_min", "GXHG_95diff_r", "GXHG_FWHM", "GXHG_absmax", "GXHG_max")
    GXHG_feats_to_use <- GXHG_feats_all[GXHG_feats_all %in% feats]
    GXHG_start <- as.numeric(Sys.time())
    if (length(GXHG_feats_to_use)) {
      GXHG_result <- G3cross(cluster[[2]],
        i = dopant_formula, j = host_formula,
        correction = "km", rmax = maxGXHGr, nrval = nGXr
      )
      GXHG_features <- g3Xfeatures(GXHG_result$r,
        gvals_new = GXHG_result$km,
        gvals_old = rrl$rrl_GXHG$mmean, type = "HG"
      )
      features <- c(features, GXHG_features)
    }
    GXHG_stop <- as.numeric(Sys.time())
    GXHG_time <- GXHG_stop - GXHG_start

    # only keep desired features
    features <- features[feats]

    met_times <- c(G_time, F_time, K_time, GXGH_time, GXHG_time)
    names(met_times) <- c("G_time", "F_time", "K_time", "GXGH_time", "GXHG_time")
    out <- unlist(c(params, jitter, features, met_times))




    rm(cluster)
    gc()
    # print(res)
  }

  cycle_time <- as.numeric(Sys.time()) - cycle_time_start
  train_gen_time <- as.numeric(Sys.time()) - as.numeric(total_time_start)
  out <- c(out, "cycle_time" = cycle_time, "train_gen_time" = train_gen_time)
  return(out)
}
