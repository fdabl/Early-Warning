library('doParallel')
source('Code/helpers.R')

test_lotka <- function(nr_transition_days = 20, new_baseline_days = 50) {
  
  lotka <- create_lotka_obj(
    nr_transition_days = nr_transition_days, nr_baselinedays = 100, noise_sd = 4
  )
  
  rs <- lotka$rs
  dat <- lotka$dat
  
  len_post <- 20 * 900
  len_change <- nr_transition_days * 900
  len_baseline <- 100 * 900
  len <- len_baseline + len_change + len_post
  
  # Test whether elements have the right length
  # 100 baseline days, 50 transition days, 20 days afterwards
  stopifnot(nrow(dat) == len)
  
  in_change <- !(rs %in% c(0.60, 1.20))
  
  # We should be len_change number of observations in transition
  stopifnot(len_change == sum(in_change))
  
  
  # Test baseline subsampling
  lotkas <- subsample_baseline(lotka, new_baseline_days = new_baseline_days)
  
  # Baseline should be shortened by new_baseline_days
  stopifnot(nrow(lotkas$dat) == len - (100 - new_baseline_days) * 900)
  stopifnot(lotkas$nr_baselinedays == 50)
  
  
  # Test frequency subsampling
  lotkaf <- subsample(lotka, freq = 90)
  
  stopifnot(nrow(lotkaf$dat) == nrow(lotka$dat) / 90)
  stopifnot(length(lotkaf$rs) == length(lotka$rs) / 90)
  stopifnot(lotkaf$freq == 90)
  
  
  # Test no transision
  lotkat <- create_lotka_obj(
    nr_transition_days = 0,
    nr_baselinedays = 100, noise_sd = 4
  )
  
  stopifnot(all(lotkat$rs == 0.60))
}


test_simulate_data <- function(times = 2) {
  noise_sds <- c(2, 4)
  nr_transition_days <- c(10, 20)
  simres <- simulate_data(
    times = times, nr_transition_days = nr_transition_days, noise_sds = noise_sds
  )
  
  stopifnot(nrow(simres) == 8)
  l1 <- simres[1, ]
  stopifnot(nrow(l1$dat) == ((100 + 10) * 900 + 20 * 900) / 10) # Subsample by 10
}

test_lotka()
test_simulate_data()
