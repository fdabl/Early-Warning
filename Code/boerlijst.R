library('doParallel')
# Simulate from the Boerlijst et al. (2013) model

simboer <- function(time, 
                    timestep,
                    muJ, muA, muP,
                    b = 1, cc = 1, 
                    noiseSD = 1, 
                    noise = TRUE, 
                    noisetype = 'A',
                    nt = 'additive') {
  
  nIter <- time / timestep
  X <- matrix(.5, nIter, 3)
  pb <- txtProgressBar(min = 2, max=nIter, initial=0, char="-", style = 3)
  
  for (i in 2:nIter) {
    
    J <- X[i-1, 1]
    A <- X[i-1, 2]
    P <- X[i-1, 3]
    
    # Independent noise added to the juveniles
    if (noisetype == 'A') {
      
      if (nt == 'additive') {
        X[i, 1] <- X[i-1, 1] + ((b*A - J / (1 + J^2) - muJ * J) + rnorm(1, 0, noiseSD)) * timestep
        X[i, 2] <- X[i-1, 2] + ((J / (1 + J^2) - A*P - muA * A)) * timestep
        X[i, 3] <- X[i-1, 3] + ((A*P*cc - muP[i] * P)) * timestep 
        
      } else {
        X[i, 1] <- X[i-1, 1] + (b*A - J / (1 + J^2) - (muJ + rnorm(1, 0, noiseSD)) * J) * timestep
        X[i, 2] <- X[i-1, 2] + (J / (1 + J^2) - A*P - muA * A) * timestep
        X[i, 3] <- X[i-1, 3] + (A*P*cc - muP[i] * P) * timestep 
      }
      
    # Independent noise added to the adults
    } else if (noisetype == 'B') {
      
      if (nt == 'additive') {
        X[i, 1] <- X[i-1, 1] + (b*A - J / (1 + J^2) - muJ * J) * timestep
        X[i, 2] <- X[i-1, 2] + ((J / (1 + J^2) - A*P - muA * A) + rnorm(1, 0, noiseSD)) * timestep
        X[i, 3] <- X[i-1, 3] + (A*P*cc - muP[i] * P) * timestep 
        
      } else {
        X[i, 1] <- X[i-1, 1] + (b*A - J / (1 + J^2) - muJ * J) * timestep
        X[i, 2] <- X[i-1, 2] + (J / (1 + J^2) - A*P - (muA + rnorm(1, 0, noiseSD)) * A) * timestep
        X[i, 3] <- X[i-1, 3] + (A*P*cc - muP[i] * P) * timestep 
      }
      
    # Independent noise added to all three population
    } else if (noisetype == 'C') {
      if (nt == 'additive') {
        X[i, 1] <- X[i-1, 1] + ((b*A - J / (1 + J^2) - muJ * J) + rnorm(1, 0, noiseSD)) * timestep 
        X[i, 2] <- X[i-1, 2] + ((J / (1 + J^2) - A*P - muA * A) + rnorm(1, 0, noiseSD)) * timestep
        X[i, 3] <- X[i-1, 3] + ((A*P*cc - muP[i] * P) + rnorm(1, 0, noiseSD)) * timestep 
        
      } else {
        X[i, 1] <- X[i-1, 1] + (b*A - J / (1 + J^2) - (muJ + rnorm(1, 0, noiseSD)) * J) * timestep
        X[i, 2] <- X[i-1, 2] + (J / (1 + J^2) - A*P - (muA + rnorm(1, 0, noiseSD)) * A) * timestep
        X[i, 3] <- X[i-1, 3] + (A*P*cc - (muP[i] + rnorm(1, 0, noiseSD)) * P) * timestep
      }
      
    # Same noise added to the death rate of all three population
    } else if (noisetype == 'D') {
      noise <- rnorm(1, 0, noiseSD)
      
      if (nt == 'additive') {
        X[i, 1] <- X[i-1, 1] + ((b*A - J / (1 + J^2) - muJ * J) + noise) * timestep
        X[i, 2] <- X[i-1, 2] + ((J / (1 + J^2) - A*P - muA * A) + noise) * timestep
        X[i, 3] <- X[i-1, 3] + ((A*P*cc - muP[i] * P) + noise) * timestep
        
      } else {
        X[i, 1] <- X[i-1, 1] + (b*A - J / (1 + J^2) - (muJ + noise) * J) * timestep
        X[i, 2] <- X[i-1, 2] + (J / (1 + J^2) - A*P - (muA + noise) * A) * timestep
        X[i, 3] <- X[i-1, 3] + (A*P*cc - (muP[i] + noise) * P) * timestep
      }
      
    } else {
      stop('Hey!')
    }
    
    setTxtProgressBar(pb, i)
  }
  
  v_time <- seq(0, time, length = nIter)
  out <- cbind(v_time, X)
  colnames(out) <- c('time', 'J', 'A', 'P')
  out
}


simulate_boerlijst <- function(noiseSD = 0.005, noisetype = 'C', nt = 'additive', cores = 8) {
  
  muA <- 0.1
  muJ <- 0.05
  muP <- seq(0.45, 0.553, length.out = 100)
  ars <- matrix(NA, nrow = length(muP), ncol = 3)
  
  registerDoParallel(cores = cores)
  
  dat <- foreach(i = seq(length(muP))) %dopar% {
    
    mup <- muP[i]
    dat <- simboer(
      time = 60000, timestep = 0.1, muJ, muA, rep(mup, 600000),
      noiseSD = noiseSD, noisetype = noisetype, nt = nt
    )
    
    # remove first 10000 as burnin
    sel <- seq(10000, nrow(dat), 20)
    dat <- dat[sel, ]
    dat <- cbind(dat, muP = mup)
    dat
  }
  
  dat
}


dat <- simulate_boerlijst(nt = 'multiplicative')
res <- matrix(NA, nrow = length(dat), ncol = 10)

for (i in seq(length(dat))) {
  dd <- dat[[i]]
  d <- dd[, -c(1, 5)]
  
  if (any(is.na(d)) | any(is.nan(d)) | any(d[, 3] < .1)) {
    print('heya!')
    res[i, ] <- c(rep(NaN, 9), dd[1, 5])
  } else {
    res[i, ] <- c(
      get_ar(d[, 1]), get_ar(d[, 2]), get_ar(d[, 3]),
      sd(d[, 1]), sd(d[, 2]), sd(d[, 3]),
      get_conn(d), eigenvalue(d), spatial_variance(d), dd[1, 5]
      )
  }
}

res <- data.frame(res)[-seq(500-4, 500), ] # last four are NaNs (already transitioned!)
names(res) <- c(
  'AR_Juvenile', 'AR_Adult', 'AR_Prey',
  'SD_Juvenile', 'SD_Adult', 'SD_Prey',
  'Cross_correlation', 'Eigenvalue', 'Spatial_variance', 'muP'
)

write.csv(res, 'Simulation-Results/boerlijst.csv', row.names = FALSE)