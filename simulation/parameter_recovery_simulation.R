### Hidden Markov model parameter recovery ###

# Author: Malte Lüken
# Date: 06.05.2020


# Global parameters -------------------------------------------------------

library(depmixS4)
library(parallel)
library(tidyverse)
library(psych)
source("simulation/model_simulation.R")
source("algorithm/model_helper_functions.R")


# Number of data sets to be generated

D <- 5


# Default response parameter values

par.fix.true <- list(fix = list(vel = c(3, 0.35), acc = c(3, 0.15)), 
                     sac = list(vel = c(3, 10), acc = c(3, 3), angle = c(0, 1)),
                     pso = list(vel = c(3, 1), acc = c(3, 3), angle = c(pi, 1)),
                     sp = list(vel = c(3, 1), acc = c(3, 0.15), angle = c(0, 1)))


# Part 1 ------------------------------------------------------------------

# Sample size

N <- 2500


# Parameter intervals

par.int.true <- cbind(seq(0.01, 0.99, length.out = D),
                      seq(1, 5, length.out = D),
                      seq(0.1, 0.6, length.out = D),
                      seq(1, 5, length.out = D),
                      seq(0.05, 0.25, length.out = D),
                      seq(1, 5, length.out = D),
                      seq(5, 15, length.out = D),
                      seq(1, 5, length.out = D),
                      seq(1, 5, length.out = D),
                      1/seq(0.1, 10, length.out = D),
                      seq(1, 5, length.out = D),
                      seq(0.5, 1.5, length.out = D),
                      seq(1, 5, length.out = D),
                      seq(1, 5, length.out = D),
                      1/seq(0.1, 10, length.out = D),
                      seq(1, 5, length.out = D),
                      seq(0.5, 1.5, length.out = D),
                      seq(1, 5, length.out = D),
                      seq(0.05, 0.25, length.out = D),
                      1/seq(0.1, 10, length.out = D))


# Parameter indicators

par.states <- c(0, rep(1, 4), rep(2:4, each = 5))
par.names <- c(0, 1, 2, 1, 2, rep(c(1, 2, 1, 2, 2), 3))
par.resps <- c(0, rep(1:2, each = 2), rep(c(1, 1, 2, 2, 3),  3))


# Parameter recovery simulation

estimates.1 <- list()

for (ns in 2:4) {
  for (p in 1:sum(par.states <= ns)) {
    
    # Prepare parallel processing
    
    ncores <- detectCores()
    clust <- makeCluster(ncores)
    clusterExport(clust, list("HMM_simulate", "cohen.kappa", "N", "par.fix.true", "par.states", "par.names", "par.resps", "p", "ns"))
    
    
    # Iterate over intervals
    
    estimates.1[[as.character(ns)]][[p]] <- parLapply(clust, par.int.true[,p], function(x, samples = N, k = ns, 
                                                                                        par.state = par.states[p], 
                                                                                        par.name = par.names[p], 
                                                                                        par.resp = par.resps[p], 
                                                                                        trueResp = par.fix.true) {
      
      # Set true parameter values
      
      if(par.name > 0 && par.resp > 0) {
        
        trueResp[[par.state]][[par.resp]][par.name] <- x
        
        trueTr <- matrix(c(0.9, 0.1, 
                           0.1, 0.9), nrow = k, ncol = k)
      } else {
        
        trueTr <- matrix(NA, nrow = k, ncol = k)
        
        diag(trueTr) <- x
        
        trueTr[is.na(trueTr)] <- (1-x)/(k-1)
      }
      
      trueIn <- rep(1/k, k)
      
      
      # Simulate data with model
      
      model <- HMM_simulate(n = samples, nstates = k, trueresp = trueResp, truetr = trueTr, truein = trueIn)
      
      model.sim <- simulate(model)
      
      class(model.sim) <- "depmix"
      
      
      # Generate and set starting values
      
      respStart <- lapply(trueResp, function(x) {
        
        lapply(x, function(y) {
          
          vec <- log(c(rgamma(1, shape = 3, scale = y[1]/2), rgamma(1, shape = 3, scale = y[2]/2)))
          
          return(vec) 
        })
      })
      
      respStart[["fix"]][["angle"]] <- c(0, 2*pi)
      respStart[["sac"]][["angle"]][1] <- 0
      respStart[["pso"]][["angle"]][1] <- pi
      respStart[["sp"]][["angle"]][1] <- 0
      
      inStart <- trueIn
      trStart <- matrix(1/k, nrow = k, ncol = k)
      trStart <- apply(trStart, 1, model.sim@transition[[1]]@family$linkfun, base = model.sim@transition[[1]]@family$base)
      
      start <- c(inStart, trStart, unlist(respStart)[1:(6*k)])
      
      names(start) <- NULL
      
      model.start <- setpars(model.sim, start)
      
      
      # Fit model
      
      model.fit <- try(fit(model.start, emc = em.control(maxit = 5000, random.start = F)))
      
      
      # Calculate accuracy
      
      acc <- try(cohen.kappa(data.frame(x = model.sim@states, y = model.fit@posterior$state))$kappa)
      
      output <- list(pars.true = getpars(model), pars.start = getpars(model.start), pars.est = try(getpars(model.fit)), accuracy = acc)
      
      return(output)
    })
    
    stopCluster(clust)
    
  }
}

beepr::beep(sound = 1)

save("estimates.1", file = "simulation/part1.Rdata")


# Part 2 ------------------------------------------------------------------

# Sample size

N <- c(500, 2500, 10000)


# Noise intervals

noise.sigma.int <- seq(1, 5, length.out = D)
noise.kappa.int <- 1/seq(0.1, 10, length.out = D)


# Parameter recovery simulation

estimates.2 <- list()

for (ns in 2:4) {
  for (ss in N) {
    
    # Prepare parallel processing
    
    ncores <- detectCores()
    clust <- makeCluster(ncores)
    clusterExport(clust, list("HMM_simulate", "cohen.kappa", "D", "par.fix.true", "noise.sigma.int", "noise.kappa.int", "ss", "ns"))
    
    
    # Iterate over intervals
    
    estimates.2[[as.character(ns)]][[as.character(ss)]] <- parLapply(clust, 1:D, function(x, samples = ss, k = ns, 
                                                                                          noise.sigma = noise.sigma.int,
                                                                                          noise.kappa = noise.kappa.int,
                                                                                          trueResp = par.fix.true) {
      
      # Set true parameter values
      
      trueTr <- matrix(c(0.9, 0.1, 
                         0.1, 0.9), nrow = k, ncol = k)
      
      trueIn <- rep(1/k, k)
      
      
      # Simulate data with model
      
      model <- HMM_simulate(n = samples, nstates = k, trueresp = trueResp, truetr = trueTr, truein = trueIn)
      
      model.sim <- simulate(model)
      
      class(model.sim) <- "depmix"
      
      
      # Add noise to data

      noise.angle <- rnorm(ntimes(model.sim), 0, noise.kappa[x])
      
      gamma_noise <- function(x, noise) {
        
        return(rgamma(1, shape = 3, scale = x/2*noise))
      }
      
      for (s in 1:nstates(model.sim)) {
          
        model.sim@response[[s]][[1]]@y <- as.matrix(apply(model.sim@response[[s]][[1]]@y, 1, 
                                                    gamma_noise, noise = noise.sigma[x]), ncol = 1)
        model.sim@response[[s]][[2]]@y <- as.matrix(apply(model.sim@response[[s]][[2]]@y, 1, 
                                                    gamma_noise, noise = noise.sigma[x]), ncol = 1)
        model.sim@response[[s]][[3]]@y <- model.sim@response[[s]][[3]]@y + noise.angle
        
      }
      
      
      # Generate and set starting values
      
      respStart <- lapply(trueResp, function(x) {
        
        lapply(x, function(y) {
          
          vec <- log(c(rgamma(1, shape = 3, scale = y[1]/2), rgamma(1, shape = 3, scale = y[2]/2)))
          
          return(vec) 
        })
      })
      
      respStart[["fix"]][["angle"]] <- c(0, 2*pi)
      respStart[["sac"]][["angle"]][1] <- 0
      respStart[["pso"]][["angle"]][1] <- pi
      respStart[["sp"]][["angle"]][1] <- 0
      
      inStart <- trueIn
      trStart <- matrix(1/k, nrow = k, ncol = k)
      trStart <- apply(trStart, 1, model.sim@transition[[1]]@family$linkfun, base = model.sim@transition[[1]]@family$base)
      
      start <- c(inStart, trStart, unlist(respStart)[1:(6*k)])
      
      names(start) <- NULL
      
      model.start <- setpars(model.sim, start)
      
      
      # Fit model
      
      model.fit <- try(fit(model.start, emc = em.control(maxit = 5000, random.start = F)))
      
      
      # Calculate accuracy
      
      acc <- try(cohen.kappa(data.frame(x = model.sim@states, y = model.fit@posterior$state))$kappa)
      
      output <- list(pars.true = getpars(model), pars.start = getpars(model.start), pars.est = try(getpars(model.fit)), accuracy = acc)
      
      return(output)
    })
    
    stopCluster(clust)
    
  }
}

beepr::beep(sound = 1)

save("estimates.2", file = "simulation/part2.Rdata")


# Part 3 ------------------------------------------------------------------

# Sample size

N <- 2500


# Noise variability

b.int <- 1:3


# Parameter recovery simulation

estimates.3 <- list()

for (ns in 2:4) {
  for (b in b.int) {
    
    # Prepare parallel processing
    
    ncores <- detectCores()
    clust <- makeCluster(ncores)
    clusterExport(clust, list("HMM_simulate", "cohen.kappa", "N", "D", "par.fix.true", "ns", "b"))
    
    
    # Iterate over intervals
    
    estimates.3[[as.character(ns)]][[as.character(b)]] <- parLapply(clust, 1:D, function(x, samples = N, k = ns, B = b,
                                                                                               trueResp = par.fix.true) {
      
      # Set true parameter values
      
      trueTr <- matrix(c(0.9, 0.1, 
                         0.1, 0.9), nrow = k, ncol = k)
      
      trueIn <- rep(1/k, k)
      
      
      # Simulate data with model
      
      model <- HMM_simulate(n = samples, nstates = k, trueresp = trueResp, truetr = trueTr, truein = trueIn)
      
      model.sim <- simulate(model)
      
      class(model.sim) <- "depmix"
      
      
      # Generate and set starting values
      
      respStart <- lapply(trueResp, function(x) {
        
        lapply(x, function(y) {
          
          vec <- log(c(rgamma(1, shape = 3, scale = (y[1]/2)*B), rgamma(1, shape = 3, scale = (y[2]/2)*B)))
          
          return(vec) 
        })
      })
      
      respStart[["fix"]][["angle"]] <- c(0, 2*pi)
      respStart[["sac"]][["angle"]][1] <- 0
      respStart[["pso"]][["angle"]][1] <- pi
      respStart[["sp"]][["angle"]][1] <- 0
      
      inStart <- trueIn
      trStart <- matrix(1/k, nrow = k, ncol = k)
      trStart <- apply(trStart, 1, model.sim@transition[[1]]@family$linkfun, base = model.sim@transition[[1]]@family$base)
      
      start <- c(inStart, trStart, unlist(respStart)[1:(6*k)])
      
      names(start) <- NULL
      
      model.start <- setpars(model.sim, start)
      
      
      # Fit model
      
      model.fit <- try(fit(model.start, emc = em.control(maxit = 5000, random.start = F)))
      
      
      # Calculate accuracy
      
      acc <- try(cohen.kappa(data.frame(x = model.sim@states, y = model.fit@posterior$state))$kappa)
      
      output <- list(pars.true = getpars(model), pars.start = getpars(model.start), pars.est = try(getpars(model.fit)), accuracy = acc)
      
      return(output)
    })
    
    stopCluster(clust)
    
  }
}

beepr::beep(sound = 1)

save("estimates.3", file = "simulation/part3.Rdata")


# Part 4 ------------------------------------------------------------------

# Sample size

N <- 2500


# Number of missing data intervals

nmiss <- c(1, 3, 5)


# Length of missing data interval

miss.int <- floor(seq(1, 200, length.out = D))


# Parameter recovery simulation

estimates.4 <- list()

for (ns in 2:4) {
  for (ms in nmiss) {
    
    # Prepare parallel processing
    
    ncores <- detectCores()
    clust <- makeCluster(ncores)
    clusterExport(clust, list("HMM_simulate", "cohen.kappa", "N", "par.fix.true", "miss.int", "ns", "ms"))
    
    
    # Iterate over intervals
    
    estimates.4[[as.character(ns)]][[as.character(ms)]] <- parLapply(clust, miss.int, function(x, samples = N, k = ns, 
                                                                                          nint = ms,
                                                                                          trueResp = par.fix.true) {
      
      # Set true parameter values
      
      trueTr <- matrix(c(0.9, 0.1, 
                         0.1, 0.9), nrow = k, ncol = k)
      
      trueIn <- rep(1/k, k)
      
      
      # Simulate data with model
      
      model <- HMM_simulate(n = samples, nstates = k, trueresp = trueResp, truetr = trueTr, truein = trueIn)
      
      model.sim <- simulate(model)
      
      class(model.sim) <- "depmix"
      
      
      # Set intervals of data to missing
      
      na <- numeric(ntimes(model.sim))
      pos <- 1:ntimes(model.sim)
      
      while(mean(na == 0) > 1-((nint*x)/ntimes(model.sim))) {
        
        na <- numeric(ntimes(model.sim))
        
        for (i in 1:nint) {
          
          st <- sample(pos[na == 0 & pos <= max(pos)-x], 1)
          
          na[st:(st+x)] <- i
          
        }
      }
        
        
      for (s in 1:nstates(model.sim)) {
        
        model.sim@response[[s]][[1]]@y[na > 0,] <- NA
        model.sim@response[[s]][[2]]@y[na > 0,] <- NA
        model.sim@response[[s]][[3]]@y[na > 0,] <- NA
        
      }
      
      
      # Generate and set starting values
      
      respStart <- lapply(trueResp, function(x) {
        
        lapply(x, function(y) {
          
          vec <- log(c(rgamma(1, shape = 3, scale = y[1]/2), rgamma(1, shape = 3, scale = y[2]/2)))
          
          return(vec) 
        })
      })
      
      respStart[["fix"]][["angle"]] <- c(0, 2*pi)
      respStart[["sac"]][["angle"]][1] <- 0
      respStart[["pso"]][["angle"]][1] <- pi
      respStart[["sp"]][["angle"]][1] <- 0
      
      inStart <- trueIn
      trStart <- matrix(1/k, nrow = k, ncol = k)
      trStart <- apply(trStart, 1, model.sim@transition[[1]]@family$linkfun, base = model.sim@transition[[1]]@family$base)
      
      start <- c(inStart, trStart, unlist(respStart)[1:(6*k)])
      
      names(start) <- NULL
      
      model.start <- setpars(model.sim, start)
      
      
      # Fit model
      
      model.fit <- try(fit(model.start, emc = em.control(maxit = 5000, random.start = F)))
      
      
      # Calculate accuracy
      
      acc <- try(cohen.kappa(data.frame(x = model.sim@states, y = model.fit@posterior$state))$kappa)
      
      output <- list(pars.true = getpars(model), pars.start = getpars(model.start), pars.est = try(getpars(model.fit)), accuracy = acc)
      
      return(output)
    })
    
    stopCluster(clust)
    
  }
}

beepr::beep(sound = 1)

save("estimates.4", file = "simulation/part4.Rdata")
