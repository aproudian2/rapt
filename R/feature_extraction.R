#' Feature Extraction
#'
#' @export
feature_extraction <- function(i, pattern, rrl,
                               pcp = 0.5, cr = 4, rho1 = 0.5, rho2 = 0.01,
                               rb = 0.1, pb = 0.1, tol = 0.005,
                               maxKr = 10, nKr = 200,
                               maxGr = 5, nGr = 1000,
                               maxGXGHr = 3, maxGXHGr = 8, nGXr = 1000,
                               vside = 0.3,
                               total_time_start = as.numeric(Sys.time())){
  #total_time_start = Sys.time()
  cycle_time_start = as.numeric(Sys.time())
  set.seed(i)
  print(paste("Interation:", i, ",  PCP:", pcp, ",  tol:", tol))
  cluster <- clustersim(pattern, pattern, 0.5,
                        pcp = pcp,
                        cr = cr,
                        rho1 = rho1,
                        rho2 = rho2,
                        rb = rb,
                        pb = pb,
                        tol = tol,
                        s = i)
  if(is.numeric(cluster)){
    out = c(params[i,], rep(NA, 18))
  }
  else {
    rrl_K = rrl$rrl_K
    rrl_G = rrl$rrl_G
    rrl_F = rrl$rrl_F
    rrl_GXGH = rrl$rrl_GXGH
    rm(rrl)
    rm(pattern)

    # marks must be factor for G3cross to work
    # convert all B to A in cluster[[2]]
    cluster[[2]]$data$marks[cluster[[2]]$data$marks == "B"] = dopant_formula
    marks(cluster[[2]]) = as.factor(marks(cluster[[2]]))
    # Measure T(r) on the pattern
    ##
    K_start = as.numeric(Sys.time())
    K3 = K3est(cluster[[1]], rmax = maxKr, nrval = nKr, correction = "trans")
    T3 = sqrt(K3$trans)- sqrt(rrl_K$mmean)
    features_K = k3features(K3$r, T3, toplot = TRUE, xlim = c(0, maxKr))
    features_K <- c(features_K[[1]], features_K[[2]], features_K[[3]], features_K[[4]], features_K[[5]])
    K_stop = as.numeric(Sys.time())
    K_time = K_stop-K_start
    ##

    # Extract Fmet
    F_start = as.numeric(Sys.time())
    result_F <- F3est(cluster[[1]], rmax = maxGr, nrval =nGr, correction = "km")
    features_F = f3features(rrl_F$r, result_F$km, rrl_F$mmean)
    F_stop = as.numeric(Sys.time())
    F_time = F_stop - F_start

    # Extract Gcross_met for  A to C
    GXGH_start = as.numeric(Sys.time())
    result_Gcross_GH = G3cross(cluster[[2]], i= dopant_formula, j = host_formula,
                               correction = "km", rmax = maxGXGHr, nrval = nGXr)


    features_Gcross_GH = g3Xfeatures(result_Gcross_GH$r,
                                     gvals_new = result_Gcross_GH$km,
                                     gvals_old = rrl_GXGH$mmean)
    GXGH_stop = as.numeric(Sys.time())
    GXGH_time = GXGH_stop-GXGH_start



    # Extract Gmet
    g_start = as.numeric(Sys.time())
    result_G <- G3est(cluster[[1]], rmax = maxGr, nrval =nGr, correction = "km")
    features_G = g3features(rrl_G$r, result_G$km, rrl_G$mmean)
    g_stop = as.numeric(Sys.time())
    g_time = g_stop - g_start
    p = params[i,]

    met_times = c(K_time, F_time, GXGH_time, g_time)
    out = c(p, features_G, features_F, features_K, features_Gcross_GH, met_times)


    rm(cluster, rvals, gvals)
    gc()
    #print(res)
  }

  cycle_time = as.numeric(Sys.time())- cycle_time_start
  train_gen_time = as.numeric(Sys.time()) - as.numeric(total_time_start)
  out = c(out, cycle_time, train_gen_time)
  return(out)
}

