### Hidden Markov model parameter recovery exploration ###

# Author: Malte Lüken
# Date: 31.05.2020


# Global parameters -------------------------------------------------------

library(depmixS4)
library(parallel)
library(tidyverse)
library(psych)
source("simulation/model_simulation.R")
source("algorithm/model_helper_functions.R")


# Record starting time

time.begin <- Sys.time()


# Number of data sets to be generated

D <- 100


# Default response parameter values

par.fix.true <- list(fix = list(vel = c(3, 0.35), acc = c(3, 0.15)), 
                     sac = list(vel = c(3, 10), acc = c(3, 3), angle = c(0, 1)),
                     pso = list(vel = c(3, 1), acc = c(3, 3), angle = c(pi, 1)),
                     sp = list(vel = c(3, 1), acc = c(3, 0.15), angle = c(0, 1)))


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
      
      # Generate data seed
      
      seed.data <- k*1e6 + b*1e4 + x*1e2 + 2e7
      
      
      # Set true parameter values
      
      trueTr <- matrix(c(0.9, 0.1, 
                         0.1, 0.9), nrow = k, ncol = k)
      
      trueIn <- rep(1/k, k)
      
      
      # Simulate data with model
      
      model <- HMM_simulate(n = samples, nstates = k, trueresp = trueResp, truetr = trueTr, truein = trueIn)
      
      model.sim <- simulate(model, seed = seed.data)
      
      class(model.sim) <- "depmix"
      
      
      # Generate and set starting values
      
      respStart <- lapply(1:length(trueResp), function(x, event = trueResp) {
        
        lapply(1:length(event[[x]]), function(y, resp = event[[x]]) {
          
          set.seed(seed.data + x*10 + y)
          
          vec <- log(c(rgamma(1, shape = 3, scale = resp[[y]][1]/2), rgamma(1, shape = 3, scale = resp[[y]][2]/2)))
          
          return(vec) 
        })
      })
      
      respStart[[1]][[3]] <- c(0, 2*pi)
      respStart[[2]][[3]][1] <- 0
      respStart[[3]][[3]][1] <- pi
      respStart[[4]][[3]][1] <- 0
      
      inStart <- trueIn
      trStart <- matrix(1/k, nrow = k, ncol = k)
      trStart <- apply(trStart, 1, model.sim@transition[[1]]@family$linkfun, base = model.sim@transition[[1]]@family$base)
      
      start <- c(inStart, trStart, unlist(respStart)[1:(6*k)])
      
      names(start) <- NULL
      
      model.start <- setpars(model.sim, start)
      
      
      # Fit model
      
      model.fit <- try(fit(model.start, emc = em.control(maxit = 5000, random.start = F)))
      
      
      # Calculate accuracy
      
      states <- try(data.frame(x = model.sim@states, y = model.fit@posterior$state))
      
      output <- list(pars.true = getpars(model), pars.start = getpars(model.start), pars.est = try(getpars(model.fit)), states = states)
      
      return(output)
    })
    
    stopCluster(clust)
    
  }
}

beepr::beep(sound = 1)

save("estimates.3", file = "simulation/part3_expl.Rdata")