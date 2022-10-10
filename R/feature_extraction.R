
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


test = lapply(1:nrow(params), function(n) {
  feature_extraction(i = n, pattern = pp3_pois, rrl = rrl_pois_list[[pattern_ind]],
  pcp = pcp, cr = params[n, 1], rho1 = params[n, 2], rho2 = params[n, 3],
  rb = params[n, 4], pb = params[n, 5],
  maxKr = maxKr, nKr = nKr, 
  maxGr = maxGr, nGr = nGr, 
  maxGXGHr = maxGXGHr, nGXr = nGXr)
})
head(params)
train_data_list <- parLapply(cl, 1:nrow(params), all_features_extract, nclust,
                             maxGr, maxKr, maxGXGHr, nGr, nKr, nGXr, 
                             rrl_pois_list[[pattern_ind]]$rrl_K, rrl_pois_list[[pattern_ind]]$rrl_G, 
                             rrl_pois_list[[pattern_ind]]$rrl_F, rrl_pois_list[[pattern_ind]]$rrl_GXGH,
                             pp3_pois_list[[pattern_ind]], pcp = pcps[[pattern_ind]], tol = tol, Sys.time())


####################################################################################
####################################################################################
## old ##
i =1
nper = 1
pattern  = pp3_all[[pattern_ind]]
rrl_K = rrl_exp_list[[pattern_ind]]$rrl_K
rrl_G = rrl_exp_list[[pattern_ind]]$rrl_G 
rrl_F = rrl_exp_list[[pattern_ind]]$rrl_F
rrl_GXGH = rrl_exp_list[[pattern_ind]]$rrl_GXGH
pcp = pcps[[pattern_ind]]
tol = tol
Sys.time()
i = 3


all_features_extract <- function(i, nper, maxGr, maxKr,maxGXGHr, nGr, nKr, nGXr,
                                 rrl_K, rrl_G, rrl_F, rrl_GXGH, pattern,
                                 pcp, tol = 0.005, total_time_start){
  #total_time_start = Sys.time()
  cycle_time_start = as.numeric(Sys.time())
  n_sets = nrow(params)
  set.seed(i)
  seeds <- round(runif(n_sets, 1, 1e8))
  print(paste("Interation:", i, ",  PCP:", pcp, ",  tol:", tol))
  cluster <- clustersim(pattern, pattern, 0.5, 
                        pcp = pcp,
                        cr = params[i,1],
                        rho1 = params[i,2],
                        rho2 = params[i,3],
                        rb = params[i,4],
                        pb = params[i,5],
                        tol = tol,
                        s = seeds[i])
  if(is.numeric(cluster)){ 
    out = c(params[i,], rep(NA, 18))
  }
  else {
    # marks must be factor for G3cross to work 
    # convert all B to A in cluster[[2]]
    cluster[[2]]$data$marks[cluster[[2]]$data$marks == "B"] = dopant_formula
    marks(cluster[[2]]) = as.factor(marks(cluster[[2]]))
    # Measure T(r) on the pattern
    ##
    K_start = as.numeric(Sys.time())
    K3 = K3est(cluster[[1]], rmax = maxKr, nrval = nKr, correction = "trans")
    T3 = sqrt(K3$trans)- sqrt(rrl_K$mmean)
    features_K = k3features(K3$r, T3, toplot = TRUE, xlim = c(0, 30))
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
    
    
    rm(cluster)
    gc()
    #print(res)
  }
  
  cycle_time = as.numeric(Sys.time())- cycle_time_start
  train_gen_time = as.numeric(Sys.time()) - as.numeric(total_time_start)
  out = c(out, cycle_time, train_gen_time)
  return(out)
}

K3 = K3est(cluster[[1]], rmax = maxKr, nrval = nKr, correction = "trans")
T3 = sqrt(K3$trans)- sqrt(rrl_K$mmean)
features_K = k3features(K3$r, T3, toplot = TRUE, xlim = c(0, 30))
features_K <- c(features_K[[1]], features_K[[2]], features_K[[3]], features_K[[4]], features_K[[5]])

rvals_new = K3$r
tvals_new = T3
plot(peak.info$y.hat)
span
peak.info

k3features <- function(rvals_new, tvals_new, toplot, ...) {
  if(any(is.infinite(tvals_new))){
    return(list(NA, NA, NA, NA, NA))
  }
  
  peak.info <- argmax(rvals_new, tvals_new, w = 3, span = 0.08)
  
  if(is.na(peak.info$x[1])){
    return(list(NA, NA, NA, NA, NA))
  }
  span <- (peak.info$x[1]/7)*(0.3)
  peak.info <- argmax(rvals_new, tvals_new, w = 3, span = span)
  peak.info$neg <- argmax(rvals_new, -1*tvals_new, w = 3, span = span)
  
  peak.info$deriv <- finite_deriv(rvals_new, peak.info$y.hat)
  
  peak.info$derivsm <- argmax(rvals_new, -1*peak.info$deriv, w = 3,
                              span = span)
  peak.info$derivsm_neg <- argmax(rvals_new, peak.info$deriv, w = 3,
                                  span = span)
  
  peak.info$dderiv <- finite_deriv(rvals_new, -1*peak.info$derivsm$y.hat)
  peak.info$dderivsm <- argmax(rvals_new, peak.info$dderiv, w = 3,
                               span = span)
  
  
  peak.info$ddderiv <- finite_deriv(rvals_new, peak.info$dderivsm$y.hat)
  peak.info$ddderivsm <- argmax(rvals_new, peak.info$ddderiv, w = 3,
                                span = span)
  
  lb <- peak.info$i[1]
  ub <- (peak.info$derivsm$i[1] + 2*peak.info$neg$i[1])/3
  if(is.na(ub)) {
    ub <- length(rvals_new)
  }
  Rddm_ind <- peak.info$ddderivsm$i[peak.info$ddderivsm$i >
                                      lb & peak.info$ddderivsm$i < ub][1]
  
  # stuff if you want to plot
  #browser()
  if(toplot == TRUE) {
    plot(rvals_new, tvals_new, type = "n", ...)
    
    lines(rvals_new, peak.info$y.hat, lwd = 2)
    points(peak.info$x[1], tvals_new[peak.info$i[1]], pch = 17, cex = 2, col= "gray60")
    #points(peak.info$neg$x, tvals_new[peak.info$neg$i], pch = 17, cex = 2, col="black")
    
    lines(rvals_new, -peak.info$derivsm$y.hat, lwd = 2, col = "red")
    points(peak.info$derivsm$x, peak.info$deriv[peak.info$derivsm$i],
           col="hotpink", cex = 2, pch = 17)
    
    lines(rvals_new, peak.info$dderivsm$y.hat, lwd = 2, col = "purple")
    #points(peak.info$dderivsm$x, peak.info$dderiv[peak.info$dderivsm$i],
    #col="green", cex = 2, pch = 17)
    
    lines(rvals_new, peak.info$ddderivsm$y.hat,lwd = 2, col = "blue")
    points(rvals_new[Rddm_ind], peak.info$ddderiv[Rddm_ind],
           col="deepskyblue", cex = 2, pch = 17)
    
    legend(max(rvals_new)*0.75, max(tvals_new)*0.9, legend = c('Tmax/Rmax', 'Tdmin/Rdmin','Rd3max'),
           col= c('gray60', 'hotpink', 'deepskyblue'), pch = 17)
  }
  
  return(list(peak.info$y.hat[peak.info$i[1]], #Km
              peak.info$x[1], #Rm
              peak.info$derivsm$x[1], #Rdm
              rvals_new[Rddm_ind[1]], #Rddm
              -peak.info$derivsm$y.hat[peak.info$derivsm$i[1]])) #Kdm
}
