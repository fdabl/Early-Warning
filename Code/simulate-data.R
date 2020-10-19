library('doParallel')
source('Code/helpers.R')

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

start <- Sys.time()
times <- 500
noise_sds <- c(4, 6, 8, 10)
nr_transition_days <- c(0, 10, 25, 50)
simres <- simulate_data(
  times = times, nr_transition_days = nr_transition_days,
  noise_sds = noise_sds, nr_baselinedays = 100, noisetype = 'additive'
)

saveRDS(simres, 'Simulation-Results/simres-500.RDS')
end <- Sys.time()
print(end - start)