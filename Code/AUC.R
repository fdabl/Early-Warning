library('MESS')
library('dplyr')
library('readr')
library('doParallel')
source('Code/helpers.R')

# Too big for Github, please send an email if you
# do not want to run all the simulations yourself
dat <- read.csv('Simulation-results/results-roc.csv')

write_auc <- function(dat, iter) {
  config <- expand.grid(
    noise = c(4, 6, 8), freq = c(9, 18, 90), nr_transition_days = c(10, 25, 50),
    window_size_days = c(10, 25, 50), nr_baselinedays = c(25, 50, 100), ews = unique(dat$ews)
  )
  
  ews <- levels(config$ews)
  ewscur <- ews[iter]
  dews <- dplyr::filter(dat, ews == !!ewscur)
  config <- dplyr::filter(config, ews == !!ewscur)
  
  start <- Sys.time()
  res <- c()
  
  for (i in seq(nrow(config))) {
    row <- config[i, ]
    noise_sd <- row[1, 1]
    freq <- row[1, 2]
    transition_days <- row[1, 3]
    window_size_days <- row[1, 4]
    nr_baselinedays <- row[1, 5]
    
    if (window_size_days < nr_baselinedays) {
      d <- dplyr::filter(
        dews, noise == !!noise_sd, freq == !!freq, nr_transition_days %in% c(NA, transition_days), ews == !!ewscur,
        nr_baselinedays == !!nr_baselinedays, window_size_days == !!window_size_days 
      )
      res <- rbind(
        res,
        rbind(
          cbind(
            create_rocdat(d), noise = noise_sd,
            freq = freq, nr_transition_days = transition_days,
            window_size_days = window_size_days, nr_baselinedays = nr_baselinedays, ews = ewscur
            )
        )
      )
    }
  }
  
  end <- Sys.time()
  print(end - start)
  
  write.csv(res, paste0('Simulation-Results/AUC-', ewscur, '.csv'), row.names = FALSE)
}


# Compute AUC for each early warning indicator
registerDoParallel(cores = 7)
foreach(i = seq(1, 11)) %dopar% { write_auc(dat, i) }
