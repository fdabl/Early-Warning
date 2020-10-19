source('Code/helpers.R')

# Too big for Github, please send an email if you
# do not want to run all the simulations yourself
simres <- readRDS('Simulation-Results/simres-500.RDS')
roc_analysis <- FALSE

sigmas <- seq(.25, 6, .25)
freqs <- c(90, 180, 900) / 10
nr_baselinedays <- c(25, 50, 100)
rolling_window_sizes <- c(10, 25, 50)


if (roc_analysis) {
  ews <- list(
    'Autocorrelation' = get_ar,
    'Cross-correlation' = get_conn,
    'Variance' = var,
    'Skewness' = skewness,
    'Kurtosis' = kurtosis,
    'CovEigen' = eigenvalue,
    'Spatial-Variance' = spatial_variance,
    'Spatial-Skewness' = spatial_skewness,
    'Spatial-Kurtosis' = spatial_kurtosis,
    'Mean' = mean,
    'Combined-Indicator' = NULL
  )
  
  start <- Sys.time()
  doParallel::registerDoParallel(detectCores())
  dat <- analyze_data(
    simres, freqs, rolling_window_sizes,
    nr_baselinedays,ews, sigmas, detect = 'baseline'
  )
  end <- Sys.time()
  print(end - start)
  
  write.csv(dat, 'Simulation-Results/results-roc.csv', row.names = FALSE)
  
} else {
  ews <- list(
    'Autocorrelation' = get_ar,
    'Cross-correlation' = get_conn,
    'Variance' = var,
    'CovEigen' = eigenvalue,
    'Combined-Indicator' = NULL
  )
  
  dat <- analyze_data_transitions(
    simres, freqs, rolling_window_sizes,
    nr_baselinedays, ews, sigmas, detect = 'baseline'
  )
  
  write.csv(dat, 'Simulation-Results/results-advance.csv', row.names = FALSE)
}
