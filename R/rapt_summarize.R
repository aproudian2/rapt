#' rescale_pattern Function
#'
#' @description
#' Rescales an input 3D point pattern
#'
#' @details
#' Takes an input 3D point pattern (pp3 object) and scales it to have intensity of new.intensity by multiplying
#' the coordinates of each point as well as the window.
#'
#' @param pattern 3D point pattern (pp3 object) to be scaled
#' @param new_intensity final intensity of point pattern
#'
#' @return pp3_new will preserve all properties of the input \emph{pattern}, except the intensity
#' @export
rescale_pattern <- function(pattern, new_intensity) {
  intensity <- sum(spatstat.geom::intensity(pattern))
  if (is.pp3(pattern)) {
    factor <- (intensity / new_intensity)^(1 / 3)
    coo <- coords(pattern)
    coo_scaled <- coo * factor

    win_scaled <- box3(
      (domain(pattern)$xrange * factor),
      (domain(pattern)$yrange * factor),
      (domain(pattern)$zrange * factor)
    )
    out <- pp3(coo_scaled$x, coo_scaled$y, coo_scaled$z, win_scaled)
  } else if (is.ppp(pattern)) {
    factor <- (intensity / new_intensity)^(1 / 2)
    coo <- coords(pattern)
    coo_scaled <- coo * factor
    win_scaled <- owin(
      (domain(pattern)$xrange * factor),
      (domain(pattern)$yrange * factor)
    )
    out <- ppp(coo_scaled$x, coo_scaled$y, win_scaled)
  }
  marks(out) <- marks(pattern)
  return(out)
}


#' Calculate T value
#'
#'
#' @description
#' Takes input of relabeling object and calculates the T value for either K or G function
#'
#' @details This uses the method from \emph{Baddeley et al.} to calculate the T value for an envelope for either
#' K  (\emph{K3est}) G \emph{G3est} function
#' @references
#' Baddeley, A., Diggle, P. J., Hardegen, A., Lawrence, T_, Milne, R. K., & Nair, G. (2014).
#' On tests of spatial pattern based on simulation envelopes. Ecological Monographs, 84(3), 477–489. https://doi.org/10.1890/13-2042.1
#' @export
calc_T_vals <- function(relabelings, func = "K", rmin = 0, rmax, K_cor = "border", G_cor = "km") {
  if (func == "K") {
    if (K_cor == "trans") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_K$trans
      })
    } else if (K_cor == "iso") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_K$iso
      })
    } else if (K_cor == "border") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_K$bord
      })
    } else {
      print("Incorrect K edge correction")
    }

    r <- relabelings[[1]]$rrl_K$r
    med <- apply(mmean, 1, median)
    ind <- which(r <= rmax & r >= rmin)
    interval <- rmax / (length(r) - 1)
    T <- apply(mmean, 2, function(x) {
      T_1 <- sqrt(x) - sqrt(med)
      T_1 <- T_1[ind]
      sum(T_1^2) * interval
    })
    return(T)
  }
  if (func == "G") {
    if (G_cor == "km") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_G$km
      })
    } else if (G_cor == "rs") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_G$rs
      })
    } else {
      print("Incorrect G edge correction")
    }
    r <- relabelings[[1]]$rrl_G$r
    med <- apply(mmean, 1, median)
    ind <- which(r <= rmax & r >= rmin)
    interval <- rmax / (length(r) - 1)
    T <- apply(mmean, 2, function(x) {
      T_1 <- x - med
      T_1 <- T_1[ind]
      sum(T_1^2) * interval
    })
  }
  return(T)
}

#' Functionto take all relabelings and return averages and envelope values
#' @param relabelings output of  \code{\link[rapt]{relabel_summarize}} function
#' @param envelope.value size of envelope to compute.  Should be decimal (e.g. 0.95 = 95%)
#' @export
average_relabelings <- function(relabelings, envelope.value = .95,
                                funcs = c("K", "G"),
                                K_cor = "trans", G_cor = "km", F_cor = "km",
                                GXGH_cor = "km", GXHG_cor = "km") {
  # transform envelope value to high index (0.95 envelope will be 0.025 and 0.975)
  envelope.value <- envelope.value + (1 - envelope.value) / 2

  # make the relabelings their own individual objects

  # K
  # extract K(r) values
  if ("K" %in% funcs) {
    if (K_cor == "trans") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_K$trans
      })
    } else if (K_cor == "iso") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_K$iso
      })
    } else if (K_cor == "border") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_K$bord
      })
    } else {
      print("Incorrect K edge correction")
    }

    # order K(r) values by value
    ordered <- apply(mmean, 1, sort)

    # get index of high and low envelope values and find values at each index
    hi.ind <- round((length(relabelings) + 1) * envelope.value, 0)
    lo.ind <- round((length(relabelings) + 1) * (1 - envelope.value), 0)
    if (lo.ind == 0) {
      lo.ind <- 1
    }
    lo <- ordered[lo.ind, ]
    hi <- ordered[hi.ind, ]

    # get r values
    r <- relabelings[[1]]$rrl_K$r

    # find the median at every distance
    med <- apply(mmean, 1, median)
    rrl_K <- data.frame(r = r, mmean = med, lo = lo, hi = hi)
  }

  # Repeat for G, GX, and F
  # G
  if ("G" %in% funcs) {
    if (G_cor == "km") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_G$km
      })
    } else if (G_cor == "rs") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_G$rs
      })
    } else {
      print("Incorrect G edge correction")
    }

    r <- relabelings[[1]]$rrl_G$r
    ordered <- apply(mmean, 1, sort)
    lo <- ordered[lo.ind, ]
    hi <- ordered[hi.ind, ]
    med <- apply(mmean, 1, median)
    rrl_G <- data.frame(r = r, mmean = med, lo = lo, hi = hi)
  }

  # F
  if ("F" %in% funcs) {
    if (F_cor == "km") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_F$km
      })
    } else if (F_cor == "rs") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_F$rs
      })
    } else if (F_cor == "cs") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_F$cs
      })
    } else {
      print("Incorrect F edge correction")
    }
    r <- relabelings[[1]]$rrl_F$r
    ordered <- apply(mmean, 1, sort)
    lo <- ordered[lo.ind, ]
    hi <- ordered[hi.ind, ]
    med <- apply(mmean, 1, median)
    rrl_F <- data.frame(r = r, mmean = med, lo = lo, hi = hi)
  }

  # GXGH
  if ("GXGH" %in% funcs) {
    if (GXGH_cor == "km") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXGH$km
      })
    } else if (GXGH_cor == "rs") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXGH$rs
      })
    } else if (GXGH_cor == "han") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXGH$han
      })
    } else if (GXGH_cor == "none") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXGH$raw
      })
    } else {
      print("Incorrect GXGH edge correction")
    }

    r <- relabelings[[1]]$rrl_GXGH$r
    ordered <- apply(mmean, 1, sort)
    lo <- ordered[lo.ind, ]
    hi <- ordered[hi.ind, ]
    med <- apply(mmean, 1, median)
    rrl_GXGH <- data.frame(r = r, mmean = med, lo = lo, hi = hi)
  }

  # GXGH
  if ("GXHG" %in% funcs) {
    if (GXHG_cor == "km") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXHG$km
      })
    } else if (GXHG_cor == "rs") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXHG$rs
      })
    } else if (GXHG_cor == "han") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXHG$han
      })
    } else if (GXHG_cor == "none") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXHG$raw
      })
    } else {
      print("Incorrect GXHG edge correction")
    }

    r <- relabelings[[1]]$rrl_GXHG$r
    ordered <- apply(mmean, 1, sort)
    lo <- ordered[lo.ind, ]
    hi <- ordered[hi.ind, ]
    med <- apply(mmean, 1, median)
    rrl_GXHG <- data.frame(r = r, mmean = med, lo = lo, hi = hi)
  }
  all_funcs <- c("K", "G", "F", "GXGH", "GXHG")
  all_funcs %in% funcs
  relabs <- c("rrl_K", "rrl_G", "rrl_F", "rrl_GXGH", "rrl_GXHG")
  out <- lapply(relabs[all_funcs %in% funcs], get, envir = environment())
  names(out) <- relabs[all_funcs %in% funcs]
  return(out)
}








#' Relabel input point pattern and calcluate summary functions
#' @param seed number to use in `set.seed` for reproducibility
#' @param funcs vector of summary functions to calculate
#' @param pattern point pattern of type ppp or pp3
#' @param dopant_formula mark of points from which to calculate summary functions
#' @param host_formula marks of points to not use in summary functions (except cross type functions)
#' @param K_cor edge correction(s) to be used for K function
#' @param G_cor edge correction(s) to be used for G function
#' @param F_cor edge correction(s) to be used for F function
#' @param GXGH_cor edge correction(s) to be used for GXGH function
#' @param GXHG_cor edge correction(s) to be used for GXHG function
#' @param ... maxKr, nKr, maxGr, nGr, maxGXGHr, maxGXHGr, nGXr
#'
#' @description takes an input binary point pattern, randomly relabels it while maintainig the
#' same ratio of marks, and calculates summary functions.  Can then use  \code{\link[rapt]{average_relabelings}}
#' to find averages and envelopes
#' @return summary functions listed in `funcs` variable
#' @export
relabel_summarize <- function(seed, pattern, funcs = c("K", "G", "F", "GXGH"),
                              dopant_formula, host_formula, k = 1,
                              maxKr = 10, nKr = 200, maxGr = 5, nGr = 1000,
                              maxGXGHr = 3, maxGXHGr = 8, nGXr = 1000, vside = 0.3,
                              K_cor = "trans", G_cor = "km", F_cor = "km",
                              GXGH_cor = "km", GXHG_cor = "km") {
  # set seed to make sure each random relabeling is different
  set.seed(seed)

  # calculate total host and dopant
  dopant_total <- sum(marks(pattern) == dopant_formula)
  host_total <- sum(marks(pattern) == host_formula)

  # relabel the input pattern, while maintaining the same proportion of host and dopant points
  relabeled <- rlabel(pattern,
    labels = as.factor(c(
      rep(dopant_formula, dopant_total),
      rep(host_formula, host_total)
    )), permute = TRUE
  )

  # select only the dopant type points
  dopant_relabeled <- subset(relabeled, marks == dopant_formula)

  # calculate summary functions
  if (is.pp3(pattern)) {
    if ("K" %in% funcs) {
      if (K_cor == "border") {
        rrl_K <- bK3est(dopant_relabeled, rmax = maxKr, nrval = nKr)
      } else {
        rrl_K <- K3est(dopant_relabeled, rmax = maxKr, nrval = nKr, correction = K_cor)
      }
    }
    if ("G" %in% funcs) {
      rrl_G <- G3est(dopant_relabeled, rmax = maxGr, nrval = nGr, correction = G_cor)
    }
    if ("F" %in% funcs) {
      rrl_F <- F3est(dopant_relabeled, rmax = maxGr, nrval = nGr, correction = F_cor, vside = vside)
    }
    if ("GXGH" %in% funcs) {
      rrl_GXGH <- G3cross(relabeled,
        i = dopant_formula, j = host_formula,
        rmax = maxGXGHr, nrval = nGXr, correction = GXGH_cor
      )
    }
    if ("GXHG" %in% funcs) {
      rrl_GXHG <- G3cross(relabeled,
        i = host_formula, j = dopant_formula,
        rmax = maxGXHGr, nrval = nGXr, correction = GXHG_cor
      )
    }
  } else if (is.ppp(pattern)) {
    if ("K" %in% funcs) {
      rrl_K <- Kest(dopant_relabeled, r = seq(0, maxKr, length.out = nKr), correction = K_cor)
    }
    if ("G" %in% funcs) {
      rrl_G <- Gest_nn(dopant_relabeled, r = seq(0, maxGr, length.out = nGr), correction = G_cor, k = k)
    }
    if ("F" %in% funcs) {
      rrl_F <- Fest(dopant_relabeled, r = seq(0, maxGr, length.out = nGr), correction = F_cor)
    }
    if ("GXGH" %in% funcs) {
      rrl_GXGH <- Gcross(relabeled,
        i = dopant_formula, j = host_formula,
        r = seq(0, maxGXGHr, length.out = nGXr), correction = GXGH_cor
      )
    }
    if ("GXHG" %in% funcs) {
      rrl_GXHG <- Gcross(relabeled,
        i = host_formula, j = dopant_formula,
        r = seq(0, maxGXHGr, length.out = nGXr), correction = GXHG_cor
      )
    }
  }


  all_funcs <- c("K", "G", "F", "GXGH", "GXHG")
  all_funcs %in% funcs
  relabs <- c("rrl_K", "rrl_G", "rrl_F", "rrl_GXGH", "rrl_GXHG")
  out <- lapply(relabs[all_funcs %in% funcs], get, envir = environment())
  names(out) <- relabs[all_funcs %in% funcs]
  return(out)
}

#' function to calculate all summary functions on point pattern
#' @param funcs vector of summary functions to calculate
#' @param pattern point pattern of type ppp or pp3
#' @param dopant_formula mark of points from which to calculate summary functions
#' @param host_formula marks of points to not use in summary functions (except cross type functions)
#' @param K_cor edge correction(s) to be used for K function
#' @param G_cor edge correction(s) to be used for G function
#' @param F_cor edge correction(s) to be used for F function
#' @param GXGH_cor edge correction(s) to be used for GXGH function
#' @param GXHG_cor edge correction(s) to be used for GXHG function
#' @param ... maxKr, nKr, maxGr, nGr, maxGXGHr, maxGXHGr, nGXr
#'
#' @description Cacluate that summary functions included in `funcs` on `pattern`, a point pattern
#' of class ppp or pp3.
#' @export
calc_summary_funcs <- function(pattern, funcs = c("K", "G", "F", "GXGH"),
                               dopant_formula, host_formula, k = 1,
                               maxKr = 10, nKr = 200, maxGr = 5, nGr = 1000,
                               maxGXGHr = 3, maxGXHGr = 8, nGXr = 1000, vside = 0.3,
                               K_cor = "trans", G_cor = "km", F_cor = "km",
                               GXGH_cor = "km", GXHG_cor = "km") {
  # select only the dopant type points
  pattern_dopant <- subset(pattern, marks == dopant_formula)

  if (is.pp3(pattern)) {
    # calculate summary functions
    if ("K" %in% funcs) {
      if (K_cor == "border") {
        K <- bK3est(pattern_dopant, rmax = maxKr, nrval = nKr)
      } else {
        K <- K3est(pattern_dopant, rmax = maxKr, nrval = nKr, correction = K_cor)
      }
    }
    if ("G" %in% funcs) {
      G <- G3est(pattern_dopant, rmax = maxGr, nrval = nGr, correction = G_cor)
    }
    if ("F" %in% funcs) {
      F <- F3est(pattern_dopant, rmax = maxGr, nrval = nGr, correction = F_cor, vside = vside)
    }
    if ("GXGH" %in% funcs) {
      GXGH <- G3cross(pattern,
        i = dopant_formula, j = host_formula,
        rmax = maxGXGHr, nrval = nGXr, correction = GXGH_cor
      )
    }
    if ("GXHG" %in% funcs) {
      GXHG <- G3cross(pattern,
        i = host_formula, j = dopant_formula,
        rmax = maxGXHGr, nrval = nGXr, correction = GXHG_cor
      )
    }
  } else if (is.ppp(pattern)) {
    if ("K" %in% funcs) {
      K <- Kest(pattern_dopant, r = seq(0, maxKr, length.out = nKr), correction = K_cor)
    }
    if ("G" %in% funcs) {
      G <- Gest_nn(pattern_dopant, r = seq(0, maxGr, length.out = nGr), correction = G_cor, k = k)
    }
    if ("F" %in% funcs) {
      F <- Fest(pattern_dopant, r = seq(0, maxGr, length.out = nGr), correction = F_cor)
    }
    if ("GXGH" %in% funcs) {
      GXGH <- Gcross(pattern,
        i = dopant_formula, j = host_formula,
        r = seq(0, maxGXGHr, length.out = nGXr), correction = GXGH_cor
      )
    }
    if ("GXHG" %in% funcs) {
      GXHG <- Gcross(pattern,
        i = host_formula, j = dopant_formula,
        r = seq(0, maxGXHGr, length.out = nGXr), correction = GXHG_cor
      )
    }
  }

  all_funcs <- c("K", "G", "F", "GXGH", "GXHG")
  relabs <- c("K", "G", "F", "GXGH", "GXHG")
  out <- lapply(relabs[all_funcs %in% funcs], get, envir = environment())
  names(out) <- relabs[all_funcs %in% funcs]
  return(out)
}



#' Plot summary functions
#' @param func Summary function to plot
#' @param observed_values observed summary function value
#' @param envelopes theoretical values found by function such as `calc_summary_funcs().` Should be list
#' @param pattern.colors colors and names for envelope lines
#' @param fill.colors colors and names for envelope fill.  Match names to `pattern.colors`
#' @param K_cor edge correction(s) to be used for K function
#' @param G_cor edge correction(s) to be used for G function
#' @param F_cor edge correction(s) to be used for F function
#' @param GXGH_cor edge correction(s) to be used for GXGH function
#' @param GXHG_cor edge correction(s) to be used for GXHG function
#'
#' @description Plot the observed value with envelopes for expected values for summary function
#' @export

plot_summary <- function(func = "K",
                         observed_values, envelopes, ...,
                         pattern.colors = c("95.0% AI" = "pink", "Median" = "black", "Observed" = "blue"),
                         fill.colors = pattern.colors, unit = "nm",
                         K_cor = "trans", G_cor = "km", F_cor = "km",
                         GXGH_cor = "km", GXHG_cor = "km") {
  if (func == "K") {
    long <- as.data.frame(envelopes[[1]]$rrl_K)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_K)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }

    baseline <- sqrt(envelopes[[1]]$rrl_K$mmean)

    observed <- as.data.frame(observed_values$K)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, K_cor],
      lo = observed[, K_cor], hi = observed[, K_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = sqrt(lo) - baseline,
      ymax = sqrt(hi) - baseline, color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)
    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(K)[g](r)))
  }

  ## plot for G
  else if (func == "G") {
    long <- as.data.frame(envelopes[[1]]$rrl_G)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    head(long)
    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_G)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }

    baseline <- (envelopes[[1]]$rrl_G$mmean)

    observed <- as.data.frame(observed_values$G)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, G_cor],
      lo = observed[, G_cor], hi = observed[, G_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = (lo) - baseline,
      ymax = (hi) - baseline, color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(G)[g](r)))
  }


  ### Plot for F
  else if (func == "F") {
    long <- as.data.frame(envelopes[[1]]$rrl_F)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_F)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }

    baseline <- (envelopes[[1]]$rrl_F$mmean)

    observed <- as.data.frame(observed_values$F)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, F_cor],
      lo = observed[, F_cor], hi = observed[, F_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = (lo) - baseline,
      ymax = (hi) - baseline, color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(F)[g](r)))
  }

  ### Plot for GXGH
  else if (func == "GXGH") {
    long <- as.data.frame(envelopes[[1]]$rrl_GXGH)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    head(long)
    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_GXGH)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }

    baseline <- (envelopes[[1]]$rrl_GXGH$mmean)

    observed <- as.data.frame(observed_values$GXGH)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, GXGH_cor],
      lo = observed[, GXGH_cor], hi = observed[, GXGH_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = (lo) - baseline,
      ymax = (hi) - baseline, color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(G)[gh](r)))
  }

  ### Plot for GXHG
  else if (func == "GXHG") {
    long <- as.data.frame(envelopes[[1]]$rrl_GXHG)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    head(long)
    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_GXHG)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }

    baseline <- (envelopes[[1]]$rrl_GXHG$mmean)

    observed <- as.data.frame(observed_values$GXHG)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, GXHG_cor],
      lo = observed[, GXHG_cor], hi = observed[, GXHG_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = (lo) - baseline,
      ymax = (hi) - baseline, color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(G)[hg](r)))
  }
}

#' Plot summary functions
#' @param func Summary function to plot
#' @param observed_values observed summary function value
#' @param envelopes theoretical values found by function such as `calc_summary_funcs().` Should be list
#' @param pattern.colors colors and names for envelope lines
#' @param fill.colors colors and names for envelope fill.  Match names to `pattern.colors`
#' @param K_cor edge correction(s) to be used for K function
#' @param G_cor edge correction(s) to be used for G function
#' @param F_cor edge correction(s) to be used for F function
#' @param GXGH_cor edge correction(s) to be used for GXGH function
#' @param GXHG_cor edge correction(s) to be used for GXHG function
#'
#' @description Plot the observed value with envelopes for expected values for summary function
#' @export
plot_summary_raw <- function(func = "K",
                             observed_values, envelopes, ...,
                             pattern.colors = c("95.0% AI" = "pink", "Median" = "black", "Observed" = "blue"),
                             fill.colors = pattern.colors, unit = "nm",
                             K_cor = "trans", G_cor = "km", F_cor = "km",
                             GXGH_cor = "km", GXHG_cor = "km") {
  if (func == "K") {
    long <- as.data.frame(envelopes[[1]]$rrl_K)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_K)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }


    observed <- as.data.frame(observed_values$K)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, K_cor],
      lo = observed[, K_cor], hi = observed[, K_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = lo,
      ymax = hi, color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)
    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(K[g](r)))
  }

  ## plot for G
  else if (func == "G") {
    long <- as.data.frame(envelopes[[1]]$rrl_G)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    head(long)
    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_G)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }


    observed <- as.data.frame(observed_values$G)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, G_cor],
      lo = observed[, G_cor], hi = observed[, G_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = (lo),
      ymax = (hi), color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(G[g](r)))
  }


  ### Plot for F
  else if (func == "F") {
    long <- as.data.frame(envelopes[[1]]$rrl_F)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_F)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }


    observed <- as.data.frame(observed_values$F)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, F_cor],
      lo = observed[, F_cor], hi = observed[, F_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = (lo),
      ymax = (hi), color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(F[g](r)))
  }

  ### Plot for GXGH
  else if (func == "GXGH") {
    long <- as.data.frame(envelopes[[1]]$rrl_GXGH)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    head(long)
    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_GXGH)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }

    baseline <- (envelopes[[1]]$rrl_GXGH$mmean)

    observed <- as.data.frame(observed_values$GXGH)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, GXGH_cor],
      lo = observed[, GXGH_cor], hi = observed[, GXGH_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = (lo),
      ymax = (hi), color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(G[gh](r)))
  }

  ### Plot for GXHG
  else if (func == "GXHG") {
    long <- as.data.frame(envelopes[[1]]$rrl_GXHG)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    head(long)
    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_GXHG)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }


    observed <- as.data.frame(observed_values$GXHG)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, GXHG_cor],
      lo = observed[, GXHG_cor], hi = observed[, GXHG_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = (lo),
      ymax = (hi), color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(G[hg](r)))
  }
}
