#' Generate data from the Lotka-Volterra model
#'
#' @param time numeric sequence of time values
#' @param rs numeric sequence of r values
#' @param timestep numeric given \delta t
#' @param X0 numeric given the initial conditions
#' @param noise_sd numeric giving the extent of noise
#' @param rho numeric giving the extent of correlation between noise
#' @param pbar boolean indicating whether a progress bar should be shown
#' @returns data from the Lotka-Volterra model
generate_data <- function(
  time, rs, timestep, X0 = c(6, 6, 0.50, 0.50), mu = 1,
  noise_sd = 1, rho = 0, pbar = TRUE, noisetype = 'additive') {
  
  # Define Parameters
  p <- 4
  C <- matrix(
    c(-.2, .04, -.2, -.2,
      .04, -.2, -.2, -.2,
      -.2, -.2, -.2, .04,
      -.2, -.2, .04, -.2),
    p, p, byrow = TRUE)
  
  # Storage
  nr_iter <- time / timestep
  X <- matrix(NA, nr_iter, 4)
  X[1, ] <- X0
  
  # Set progress bar
  if(pbar) {
    pb <- txtProgressBar(min = 2, max = nr_iter, initial = 0, char = '-', style = 3)
  }
  
  S <- matrix(rho, nrow = 4, ncol = 4)
  diag(S) <- noise_sd
  
  # Solve with Euler's method
  for (j in seq(2, nr_iter)) {
    Xprev <- X[j - 1, ]
    r <- c(1, 1, rs[j], rs[j])
    
    if (rho != 0) {
      # Simulate correlated noise
      epsilon <- mvtnorm::rmvnorm(1, c(0, 0, 0, 0), S) 
      
    } else {
      epsilon <- rnorm(4, 0, noise_sd)
    }
    
    
    for (i in 1:p) {
      if (noisetype == 'additive') {
        
        dX <- r[i] * Xprev[i] + Xprev[i] * sum(C[i, ] * Xprev) + mu + epsilon[i]
        
      } else {
        
        dX <- (r[i] + epsilon[i]) * Xprev[i] + Xprev[i] * sum(C[i, ] * Xprev) + mu
      }
      
      X[j, i] <- Xprev[i] + dX * timestep
      
    }
    
    if (pbar) setTxtProgressBar(pb, j)
  }
  
  time <- seq(0, time, length = nr_iter)
  out <- cbind(time, X)
  colnames(out) <- c('time', 'x1', 'x2', 'x3', 'x4')
  out
}


#' Numerically Integrate One Step Forward
#'
#' @returns a value for X after one time step
dXdt <- function(
  i, r, epsilon, C, mu, X, timestep = 1,
  noisetype = c('additive', 'multiplicative')
) {
  
  if (noisetype == 'additive') {
    dX <- r[i] * X[i] + X[i] * sum(C[i, ] * X) + mu + epsilon
  } else {
    dX <- (r[i] + epsilon) * X[i] + X[i] * sum(C[i, ] * X) + mu
  }
  
  dX * timestep
}


#' z-standardize the data
#' @param x numeric vector
#' @returns z-standardized vector
zstand <- function(x) {
  (x - mean(x)) / sd(x)
}


#' Computes Lag-1 Auto-correlation
#' @param x numeric vector
#' @returns autocorrelation
get_ar <- function(x) {
  cor(x[-1], x[-length(x)])
}


#' Computes the Skewness
#' @param x numeric vector
#' @returns skewness
skewness <- function(x) {
  m <- mean(x)
  n <- length(x)
  mean((x - m)^3) / mean((x - m)^2)^(3/2)
}


#' Computes the Kurtosis
#' @param x numeric vector
#' @returns kurtosis
kurtosis <- function(x) {
  m <- mean(x)
  n <- length(x)
  mean((x - m)^4) / mean((x - m)^2)^2
}


#' Computes Connectivity
#' @param x matrix
#' @returns mean absolute value of all cross-correlations
get_conn <- function(x) {
  mean(abs(cor(x)))
}


#' Computes the dominant eigenvalue
#' @param x matrix
#' @returns dominant eigenvalue
eigenvalue <- function(x) {
  eigen(cov(x))$values[1]
}


#' Computes the spatial variance
#' @param x
#' @returns spatial variance
spatial_variance <- function(x) {
  x <- as.matrix(x)
  mean((x - apply(x, 2, mean))^2)
}


#' Computes the spatial skewness
#' @param x matrix
#' @returns spatial skewness
spatial_skewness <- function(x) {
  x <- as.matrix(x)
  m <- apply(x, 2, mean)
  sigma <- sqrt(mean((x - apply(x, 2, mean))^2))
  #sigma <- sqrt(spatial_variance(x))
    
  mean((x - m)^3 / sigma^3)
}


#' Computes the spatial kurtosis
#' @param x matrix
#' @returns spatial skewness
spatial_kurtosis <- function(x) {
  x <- as.matrix(x)
  nm <- 4*nrow(x)
  m <- apply(x, 2, mean)
  mean((x - m)^4) / mean((x - m)^2)^2
}


#' Create the Lotka-Volterra Data
#' 
#' @param nr_transition_days numeric indicating the number of transition days
#' @param noise_sd numeric indicating the standard deviation of the Gaussian noise
#' @param nr_baselinedays numeric indicating the number of baseline days
#' @param noisetype character indicating the type of noise
#' @param rho numeric indicating the extent of the correlation between the errors
#' @returns Lotka-Volterra object
create_lotka_obj <- function(
  nr_transition_days, noise_sd, nr_baselinedays = 100, mu = 1,
  noisetype = 'additive', rho = 0, pbar = FALSE, r_b = 1.20
) {
  
  # One day allows 900 minutes of sampling
  nr_samples <- nr_baselinedays * 10 * 90
  nr_days_after <- 20
  transition <- nr_transition_days != 0
  
  if (transition) {
    rs <- c(
      rep(0.60, nr_samples - 1), # baseline
      seq(0.60, r_b, length.out = nr_transition_days * 900 + 2), # transition period
      rep(r_b, length.out = nr_days_after * 900 - 1)
    )
    
  } else {
    nr_transition_days <- 50 # longest period
    rs <- rep(0.60, nr_samples + nr_transition_days * 900 + nr_days_after * 900)
  }
  
  dat <- generate_data(
    time = length(rs) * 0.01, rho = rho, rs = rs, mu = mu, pbar = pbar,
    timestep = 0.01, noise_sd = noise_sd, noisetype = noisetype
  )
  
  list(
    'dat' = dat,
    'rs' = rs,
    'nr_transition_days' = ifelse(transition, nr_transition_days, NA),
    'nr_baselinedays' = nr_baselinedays,
    'noise' = noise_sd,
    'freq' = 1,
    'transition' = transition,
    'noisetype' = noisetype
  )
}


#' Subsample the Lotka-Volterra Data
#' 
#' @param lotka list Lotka-Volterra object
#' @param freq numeric indicating the frequenccy
#' @returns Subsampled Lotka-Volterra Object
subsample <- function(lotka, freq) {
  res <- lotka
  
  sel <- seq(1, nrow(res$dat), freq)
  res$dat <- res$dat[sel, ]
  res$rs <- res$rs[sel]
  res$freq <- res$freq * freq
  
  res
}


#' Subsample the Lotka-Volterra Data (Baseline)
#' @param lotka_dat matrix given the simulated data
#' @param new_baseline_days numeric given the extent of new baseline
#' @returns subsampled Lotka-Volterra Data (Baseline)
subsample_baseline <- function(lotka, new_baseline_days) {
  
  stopifnot(lotka$nr_baselinedays > new_baseline_days)
  
  baseline_samples <- lotka$nr_baselinedays * 900 / lotka$freq
  new_baseline_samples <- new_baseline_days * 900 / lotka$freq
  
  remove <- seq(baseline_samples - new_baseline_samples)
  lotka$dat <- lotka$dat[-remove, ]
  
  lotka$nr_baselinedays <- new_baseline_days
  lotka$rs <- lotka$rs[-remove]
  lotka
}


#' Rolling Window
#' @param y can be multidimensional
#' @param wsize
#' @param fn
#' @returns 
rolling_window <- function(y, wsize, fn) {
  y <- data.frame(y)
  res <- numeric(nrow(y) - wsize + 1)
  
  for (i in seq(length(res) - 1)) {
    res[i] <- fn(y[seq(i, wsize + i), ])
  }
  
  res
}


#' Detects whether the indicator coefficient
#' is x sigma larger than before (using a running mean)
#' 
#' @param indicator numeric early warning indicator
#' @param wsize numeric size of rolling window
#' @param sigma numeric decision threshold
detect_change <- function(indicator, wsize, sigma = 2) {
  changes <- sapply(seq(length(indicator) - wsize), function(i) {
    
    prev <- indicator[seq(wsize, wsize + i)]
    cur <- indicator[wsize + i]
    
    cur > (mean(prev) + sigma * sd(prev))
  })
  
  c(rep(NA, sum(is.na(indicator))), changes)
}



#' Detects whether the indicator coefficient
#' is x sigma larger than the baseline (using a running mean)
#' 
#' @param indicator numeric early warning indicator
#' @param baseline_ix numeric index of baseline
#' @param sigma numeric decision threshold
detect_change_from_baseline <- function(indicator, baseline_ix, sigma = 3) {
  
  mu <- mean(indicator[baseline_ix], na.rm = TRUE)
  std <- sd(indicator[baseline_ix], na.rm = TRUE)
  len_baseline <- sum(baseline_ix)
  
  # we do not compute changes at the baseline
  # so we add len_baseline numbers of NAs
  indicator_wb <- indicator[!baseline_ix]
  
  changes <- sapply(seq(length(indicator_wb)), function(i) {
    cur <- mean(indicator_wb[seq(1, i)])
    # cur <- indicator_wb[i]
    cur > (mu + sigma * std)
  })
  
  c(rep(NA, len_baseline), changes)
}


#' Creates ROC data from simulation result
#' 
#' @param res
#' @returns data frame including false positive and true positive rates
create_rocdat <- function(res, w = 1) {
  
  conv <- function(x) as.logical(as.numeric(as.character(x)))
  res <- mutate(res, transition_trial = !is.na(nr_transition_days))
  
  d1 <- res %>% 
    filter(transition_trial == 1) %>% 
    mutate(
      tp = conv(true_positive),
      tn = conv(true_negative),
      fp = conv(false_positive),
      fn = conv(false_negative)
    ) %>% 
    group_by(sigma) %>% 
    summarize(
      tpr = sum(tp, na.rm = TRUE) / (sum(tp, na.rm = TRUE) + sum(fn, na.rm = TRUE))
    ) %>% 
    mutate(
      tpr = ifelse(is.nan(tpr), 1, tpr)
    )
  
  d2 <- res %>% 
    filter(transition_trial == 0) %>% 
    mutate(
      tp = conv(true_positive),
      tn = conv(true_negative),
      fp = conv(false_positive),
      fn = conv(false_negative)
    ) %>% 
    group_by(sigma) %>% 
    summarize(
      fpr = sum(fp, na.rm = TRUE) / (sum(fp, na.rm = TRUE) + sum(tn, na.rm = TRUE))
    )
  
  dd <- d1
  dd$fpr <- d2$fpr
  
  dd
}



#' Simulate Data
#' 
#' @param times numeric given the number of simulation runs
#' @param nr_transition_days numeric given the number of transition days 
#' @param noise_sds numeric given the extent of the noise
#' @param nr_baselinedays numeric given the number of baseline days
#' @returns list of simulated data sets
simulate_data <- function(
  times, nr_transition_days, noise_sds, rho = 0,
  noisetype = 'additive', nr_baselinedays = 100, mu = 1, r_b = 1.20
) {
  # doParallel::registerDoParallel(detectCores() - 1)
  
  comb <- expand.grid(
    nr_transition_days = nr_transition_days,
    noise_sds = noise_sds, times = seq(times)
  )
  
  ncomb <- nrow(comb)
  
  simres <- foreach(i = seq(ncomb), .combine = 'rbind',
                    .export = c('create_lotka_obj', 'generate_data', 'subsample')) %dopar% {
    nr_transition_days <- comb[i, 1]
    noise_sd <- comb[i, 2]
    transition <- comb[i, 3]
    
    dat <- create_lotka_obj(
      nr_transition_days = nr_transition_days,
      nr_baselinedays = nr_baselinedays, rho = rho,
      noisetype = noisetype, noise_sd = noise_sd, mu = mu, r_b = r_b
    )
    
    # store only every 10 minute samples to save disk space
    subsample(dat, freq = 10)
  }
  
  simres
}


#' Finds the Transitioning Point
#' 
#' @param y numeric vector given the time-series
#' @returns the index of transition
find_transition <- function(y) {
  m <- changepoint::cpt.mean(y)
  
  # automatically returns the length of the data set
  # when no transition point is found!
  m@cpts[1]
}


#' Analyse the simulated data for ROC
#' 
#' @param simres list of simulation results
#' @param freqs numeric vector of frequencies to analyze
#' @param window_sizes numeric vector of rolling window days to analyze
#' @param nr_baselinedays numeric vector of baseline days to analyze
#' @param ews list of name of early warning indicators and function
#' @param sigmas numeric vector of decision thresholds to analyze the data with
#' @param detect character given type of detection (not used)
#' @returns analysis results in a data.frame
analyze_data <- function(
  simres, freqs, window_sizes,
  nr_baselinedays, ews, sigmas = 2, detect = 'baseline'
  ) {
  # doParallel::registerDoParallel(detectCores() - 1)
  
  nsims <- nrow(simres)
  res <- foreach(
		 i = seq(nsims), .combine = 'rbind', .packages = c('changepoint'),
		 .export = c('subsample', 'subsample_baseline', 'find_transition',
			     'skewness', 'kurtosis', 'eigenvalue', 'spatial_variance',
			     'spatial_kurtosis', 'spatial_skewness', 'get_conn', 'get_ar',
			     'rolling_window', 'detect_change_from_baseline', 'detect_change', 'zstand')
		) %dopar% {

    j <- 1
    res_loop <- matrix(
      NA, ncol = 12,
      nrow = (length(freqs) * length(window_sizes) *
              length(nr_baselinedays) * length(ews) * length(sigmas))
    )
  
    colnames(res_loop) <- c(
    	'nr_transition_days', 'noise', 'freq', 'window_size_days', 'nr_baselinedays',
    	'ews', 'true_positive', 'true_negative', 'false_positive', 'false_negative', 'sigma', 'transition_id'
    )
    
    sim <- simres[i, ]
    
    nr_transition_days <- sim$nr_transition_days
    simdat <- sim$dat
    noise <- sim$noise
    
    is_transition_trial <- sim$transition
    
    for (freq in freqs) {
      lotkaf <- subsample(sim, freq = freq)
      
      for (baselinedays in nr_baselinedays) {
        if (baselinedays < 100) {
          lotkab <- subsample_baseline(lotkaf, baselinedays)
          y <- lotkab$dat[, -1]
          rs <- lotkab$rs
        } else {
          y <- lotkaf$dat[, -1]
          rs <- lotkaf$rs
        }
        
        y1 <- y[, 1]
        
        # reduce the number of data points to where the transition
        # happens when there is indeed a transition
        if (is_transition_trial) {
          transition_idx <- find_transition(y1)
        } else {
          transition_idx <- nrow(y)
        }
        
        y <- y[seq(transition_idx), ]
        y1 <- y[, 1]
        
        for (wsize_day in window_sizes) {
          if (wsize_day < baselinedays) {
            # get window_sizes in number of samples, not days
            wsize <- wsize_day * 900 / (freq * sim$freq)
            
            for (l in seq(length(ews))) {
              fn <- ews[[l]]
              ews_name <- names(ews)[l]
              
              # connectivity requires multidimensional input!
              if (ews_name %in% c('Cross-correlation', 'CovEigen', 'Spatial-Variance', 'Spatial-Skewness', 'Spatial-Kurtosis')) {
                indicator <- rolling_window(y, wsize, fn)
                indicator <- c(rep(NA, nrow(y) - length(indicator)), indicator)
                
      	      } else if (ews_name == 'Combined-Indicator') {
      	        indicator_ar <- zstand(rolling_window(y1, wsize, get_ar))
      	        indicator_var <- zstand(rolling_window(y1, wsize, var))
      	        indicator_eigen <- zstand(rolling_window(y, wsize, eigenvalue))
      	        indicator_cross <- zstand(rolling_window(y, wsize, get_conn))
            		indicator <- indicator_ar + indicator_var + indicator_eigen + indicator_cross
                indicator <- c(rep(NA, length(y1) - length(indicator)), indicator)
                
      	      } else {
                indicator <- rolling_window(y1, wsize, fn)
                indicator <- c(rep(NA, length(y1) - length(indicator)), indicator)
              }
              
              # run the detection for various sigmas to compute ROC curves
              for (sigma in sigmas) {
                
                # if baseline, then test against baseline
                if (detect == 'baseline') {
                  # baseline_ix <- rs <= 0.60 # this only works for transition trials!
                  
                  # create an baseline index
                  baseline_ix <- rep(FALSE, length(y1))
                  baseline_ix[seq(baselinedays * 900 / (freq * sim$freq))] <- TRUE
                  
                  changes <- detect_change_from_baseline(indicator, baseline_ix, sigma = sigma)
                } else {
                  changes <- detect_change(indicator, wsize, sigma = sigma)
                }
                
                # number if detected, else NA
                first <- which(changes)[1]
                
                fn = fp = tp = tn = 0
                
                # true negatives make no sense here because we always transition
                # false positive make no sense here because a signal means we anticipated
                if (is_transition_trial) {
                  tn <- NA
                  fp <- NA
                  
                  # if we do not detect a change at all, false negative!
                  if (is.na(first)) {
                    fn <- 1
                    
                  } else {
                    
                    # transition happened earlier than first detection, false negative!
                    if (first > transition_idx) {
                      fn <- 1
                    } else if (first < transition_idx) {
                      tp <- 1
                    }
                  }
                } else {
                  # we never transition, so true positives make no sense, nor do false negatives
                  tp <- NA
                  fn <- NA
                  
                  # if we do not detect a change at all, true negative!
                  # otherwise it's a false positive
                  if (is.na(first)) {
                    tn <- 1
                  } else {
                    fp <- 1
                  }
                }
                
                res_loop[j, ] <- c(
                  nr_transition_days, noise, freq, wsize_day,
                  baselinedays, ews_name, tp, tn, fp, fn, sigma, transition_idx
                )
                j <- j + 1
              }
            }
          }
        }
      }
    }
    res_loop
  }
  
  res <- data.frame(res)
  res <- res[!is.na(res$nr_baselinedays), ] # remove NAs
  res
}


#' Analyse the simulated data for how far ews occur in advance (given true positive)
#' 
#' @param simres list of simulation results
#' @param freqs numeric vector of frequencies to analyze
#' @param window_sizes numeric vector of rolling window days to analyze
#' @param nr_baselinedays numeric vector of baseline days to analyze
#' @param ews list of name of early warning indicators and function
#' @param sigmas numeric vector of decision thresholds to analyze the data with
#' @param pred_days numeric vector specifying how far in advance a signal
#' is taken as a true positive (not used, i.e., set to Inf)
#' @param detect character given type of detection (not used)
#' @returns analysis results how far in advance ews occur in a data.frame
analyze_data_transitions <- function(
  simres, freqs, window_sizes, nr_baselinedays,
  ews, sigmas = 2, pred_days = Inf, detect = 'baseline'
) {
  # doParallel::registerDoParallel(detectCores() - 1)
  
  nsims <- nrow(simres)
  res <- foreach(
    i = seq(nsims), .combine = 'rbind', .packages = c('changepoint'),
    .export = c('subsample', 'subsample_baseline', 'find_transition',
                'skewness', 'kurtosis', 'eigenvalue', 'spatial_variance',
                'spatial_kurtosis', 'spatial_skewness', 'get_conn', 'get_ar',
                'rolling_window', 'detect_change_from_baseline', 'detect_change', 'zstand')
  ) %dopar% {
    
    j <- 1
    res_loop <- matrix(
      NA,
      ncol = 11,
      nrow = (length(freqs) * length(window_sizes) * length(pred_days) *
                length(nr_baselinedays) * length(ews) * length(sigmas))
    )
    
    colnames(res_loop) <- c(
      'nr_transition_days', 'noise', 'freq', 'window_size_days', 'pred_days',
      'nr_baselinedays', 'ews', 'false_positive', 'first_signal', 'sigma', 'transition_idx'
    )
    sim <- simres[i, ]
    
    nr_transition_days <- sim$nr_transition_days
    simdat <- sim$dat
    noise <- sim$noise
    
    is_transition_trial <- sim$transition
    
    for (freq in freqs) {
      lotkaf <- subsample(sim, freq = freq)
      
      for (baselinedays in nr_baselinedays) {
        if (baselinedays < 100) {
          lotkab <- subsample_baseline(lotkaf, baselinedays)
          y <- lotkab$dat[, -1]
          rs <- lotkab$rs
          time <- lotkab$dat[, 1]
        } else {
          y <- lotkaf$dat[, -1]
          rs <- lotkaf$rs
          time <- lotkaf$dat[, 1]
        }
        
        y1 <- y[, 1]
        
        # reduce the number of data points to where the transition
        # happens when there is indeed a transition
        if (is_transition_trial) {
          transition_idx <- find_transition(y1)
        } else {
          transition_idx <- nrow(y)
        }
        
        y <- y[seq(transition_idx), ]
        y1 <- y[, 1]
        
        
        for (wsize_day in window_sizes) {
          if (wsize_day < baselinedays) {
            # get window_sizes in number of samples, not days
            wsize <- wsize_day * 900 / (freq * sim$freq)
            
            for (l in seq(length(ews))) {
              fn <- ews[[l]]
              ews_name <- names(ews)[l]
              
              # connectivity requires multidimensional input!
              if (ews_name %in% c('Cross-correlation', 'CovEigen', 'Spatial-Variance', 'Spatial-Skewness', 'Spatial-Kurtosis')) {
                indicator <- rolling_window(y, wsize, fn)
                indicator <- c(rep(NA, nrow(y) - length(indicator)), indicator)
                
              } else if (ews_name == 'Combined-Indicator') {
                indicator_ar <- zstand(rolling_window(y1, wsize, get_ar))
                indicator_var <- zstand(rolling_window(y1, wsize, var))
                indicator_eigen <- zstand(rolling_window(y, wsize, eigenvalue))
                indicator_cross <- zstand(rolling_window(y, wsize, get_conn))
                indicator <- indicator_ar + indicator_var + indicator_eigen + indicator_cross
                indicator <- c(rep(NA, length(y1) - length(indicator)), indicator)
                
              } else {
                indicator <- rolling_window(y1, wsize, fn)
                indicator <- c(rep(NA, length(y1) - length(indicator)), indicator)
              }
              
              # run the detection for various sigmas to compute ROC curves
              for (sigma in sigmas) {
                
                # if baseline, then test against baseline
                if (detect == 'baseline') {
                  
                  # create a baseline index
                  baseline_ix <- rep(FALSE, length(y1))
                  baseline_ix[seq(baselinedays * 900 / (freq * sim$freq))] <- TRUE
                  
                  changes <- detect_change_from_baseline(indicator, baseline_ix, sigma = sigma)
                } else {
                  changes <- detect_change(indicator, wsize, sigma = sigma)
                }
                
                # number if detected, else NA
                signal <- which(changes)
                
                if (length(signal) == 0) {
                  fp <- NA
                  pred_day <- NA
                  first_signal <- NA
                  res_loop[j, ] <- c(
                    nr_transition_days, noise, freq, wsize_day, pred_day,
                    baselinedays, ews_name, fp, first_signal, sigma, transition_idx
                  )
                  j <- j + 1
                } else {
                  
                  # only get signals before the transition
                  signal <- signal[signal < transition_idx]
                  first_signal <- transition_idx - signal[1]
                  
                  for (pred_day in pred_days) {
                    pred_window <- pred_day * 900 / (sim$freq * freq)
                    
                    fp <- sum((signal + pred_window) < transition_idx)
                    res_loop[j, ] <- c(
                      nr_transition_days, noise, freq, wsize_day, pred_day,
                      baselinedays, ews_name, fp, first_signal, sigma, transition_idx
                    )
                    j <- j + 1
                  }
                }
              }
            }
          }
        }
      }
    }
    res_loop
  }
  
  res <- data.frame(res)
  res <- res[!is.na(res$nr_baselinedays), ] # remove NAs
  res
}


#' Create ROC Plot
#' 
#' @param dat data.frame given the AUC data
#' @param ...
#' @returns NULL
create_rocplot <- function(dat, cex.axis = 1, cex.legend = 1, xlab = '', ylab = '',...) {
  
  ews <- unique(dat$ews)
  noise <- unique(dat$noise)
  trans_days <- unique(dat$nr_transition_days)
  trans_days <- trans_days[!is.na(trans_days)]
  
  uniq_freq <- as.numeric(as.character(unique(dat$freq)))
  nfreq <- length(uniq_freq)
  
  if (length(ews) > 1 || length(noise) > 1 || length(trans_days) > 1) {
    stop('Function requires single specification.')
  }
  
  plot(
    0, 0, pch = 20, ylim = c(0, 1), xlim = c(0, 1),
    axes = FALSE, xlab = xlab, ylab = ylab, cex = 0, ...
  )
  
  axis(1, cex.axis = cex.axis)
  axis(2, las = 1, cex.axis = cex.axis)
  lines(c(0, 1), c(0, 1), lty = 1, col = 'gray76')
  
  l <- 1
  legend_names <- c()
  cols <- RColorBrewer::brewer.pal(n = nfreq, 'Set1')
      
  for (i in seq(nfreq)) {
    
    freq <- uniq_freq[i]
    
    rocdat <- create_rocdat(
      filter(dat, freq == !!freq)
    )
    
    day <- 900 / (freq * 10)
    legend_names <- c(legend_names, ifelse(day < 10, paste0('  ', day, 'x Day'), paste0(day, 'x Day')))
    points(rocdat$fpr, rocdat$tpr, col = cols[l], pch = 20, cex = 1.50)
    lines(rocdat$fpr, rocdat$tpr, col = cols[l], lty = 1, lwd = 2)
    
    # colour the popular 2\sigma rule in black
    if (any(rocdat$sigma == 2)) {
      u <- filter(rocdat, sigma == 2)
      points(u$fpr, u$tpr, col = 'black', pch = 20, cex = 2.2)
    }
    
    l <- l + 1
  }
      
  legend(
    'bottomright',
    legend = legend_names, lty = c(1, 1, 1), lwd = 2,
    col = cols[seq(l)], box.lty = 0, bty = 'n', cex = cex.legend
  )
}
