### Hidden Markov model parameter recovery ###

# Author: Malte Lüken
# Date: 06.05.2020

library(depmixS4)
library(parallel)
source("simulation/model_simulation.R")
source("algorithm/model_helper_functions.R")


# Part 1 ------------------------------------------------------------------

D <- 5
N <- 2500
K <- 2:4

aij <- seq(0.01, 0.99, length.out = D)

par.int.true <- cbind(seq(1, 5, length.out = D),
                      seq(0.1, 0.6, length.out = D),
                      seq(1, 5, length.out = D),
                      seq(0.05, 0.25, length.out = D),
                      seq(1, 5, length.out = D),
                      seq(5, 25, length.out = D),
                      seq(1, 5, length.out = D),
                      seq(1, 5, length.out = D),
                      1/seq(0.1, 10, length.out = D),
                      seq(1, 5, length.out = D),
                      seq(1, 5, length.out = D),
                      seq(1, 5, length.out = D),
                      seq(1, 3, length.out = D),
                      1/seq(0.1, 10, length.out = D),
                      seq(1, 5, length.out = D),
                      seq(1, 5, length.out = D),
                      seq(1, 3, length.out = D),
                      seq(1, 3, length.out = D),
                      1/seq(0.1, 10, length.out = D))


par.fix.true <- list(fix = list(vel = c(3, 0.35), acc = c(3, 0.15)), 
                     sac = list(vel = c(3, 15), acc = c(3, 3), angle = c(0, 1)),
                     pso = list(vel = c(3, 3), acc = c(2, 2), angle = c(pi, 1)),
                     sp = list(vel = c(3, 1.5), acc = c(3, 0.15), angle = c(0, 1)))

par.states <- c(rep(1, 4), rep(2:4, each = 5))
par.names <- c("alpha", "beta", "alpha", "beta", rep(c("alpha", "beta", "alpha", "beta", "kappa"), 3))
par.resp <- c(rep(1:2, each = 2), rep(c(1, 1, 2, 2, 3),  3))

ncores <- detectCores()
clust <- makeCluster(ncores)
clusterExport(clust, list("HMM_simulate", "N", "par.fix.true"))
                               
estimates <- parLapply(clust, par.int.true[,1], function(x, N = 500, K = 2, k = 2, par = "alpha", met = 1, true = par.fix.true) {
  
  par.fix.true[[k]][[met]][1] <- 3
  
  instart <- rep(1/K, K)
  
  trstart <- matrix(1/K, nrow = K, ncol = K)
  
  model <- HMM_simulate(n = N, nstates = k, respstart = true, trstart = trstart, instart = instart)
  
  sim <- simulate(model)
  
  class(sim) <- "depmix"
  
  out <- fit(sim, emc = em.control(random.start = F))
  
  return(out)
})


stopCluster(clust)

exp(summary(estimates[[2]]))
