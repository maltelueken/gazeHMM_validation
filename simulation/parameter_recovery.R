### Hidden Markov model parameter recovery ###

# Author: Malte Lüken
# Date: 06.05.2020

library(depmixS4)
library(parallel)
source("simulation/model_simulation.R")


# Part 1 ------------------------------------------------------------------

D <- 100
N <- 2500
K <- 2:4

aij <- seq(0.01, 0.99, length.out = D)

alpha.vel <- matrix(NA, D, max(K))
beta.vel <- matrix(NA, D, max(K))
alpha.acc <- matrix(NA, D, max(K))
beta.acc <- matrix(NA, D, max(K))

alpha.vel[,1] <- seq(1, 5, length.out = D)
beta.vel[,1] <- seq(0.1, 0.6, length.out = D)
alpha.acc[,1] <- seq(1, 5, length.out = D)
beta.acc[,1] <- seq(0.05, 0.25, length.out = D)

alpha.vel[,2] <- seq(1, 5, length.out = D)
beta.vel[,2] <- seq(5, 25, length.out = D)
alpha.acc[,2] <- seq(1, 5, length.out = D)
beta.acc[,2] <- seq(1, 5, length.out = D)

alpha.vel[,3] <- seq(1, 5, length.out = D)
beta.vel[,3] <- seq(1, 5, length.out = D)
alpha.acc[,3] <- seq(1, 5, length.out = D)
beta.acc[,3] <- seq(1, 3, length.out = D)

alpha.vel[,4] <- seq(1, 5, length.out = D)
beta.vel[,4] <- seq(1, 5, length.out = D)
alpha.acc[,4] <- seq(1, 3, length.out = D)
beta.acc[,4] <- seq(1, 3, length.out = D)

kappa <- 1/seq(0.1, 10, length.out = D)

alpha.vel.fix <- c(3, 3, 3, 3)
beta.vel.fix <- c(0.35, 15, 3, 1.5)
alpha.acc.fix <- c(3, 3, 3, 3)
beta.acc.fix <- c(0.15, 3, 2, 0.15)

mu.fix <- c(0, pi, 0)
kappa.fix <- 1

par.true <- cbind(seq(1, 5, length.out = D),
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

respstart <- list(fix = list(vel = c(alpha.vel.fix[1], beta.vel.fix[1]), 
                                          acc = c(alpha.acc.fix[1], beta.acc.fix[1])), 
                               sac = list(vel = c(alpha.vel.fix[2], beta.vel.fix[2]), 
                                          acc = c(alpha.acc.fix[2], beta.acc.fix[2]),
                                          angle = c(mu.fix[1], kappa.fix)))

trstart <- matrix(1/2, nrow = 2, ncol = 2)
instart <- rep(1/2, 2)

ncores <- detectCores()
clust <- makeCluster(ncores)
clusterExport(clust, list("HMM_simulate", "N", "respstart", "trstart", "instart"))
                               
estimates <- parLapply(clust, alpha.vel[,1], function(x, N = 500, k = 2, 
                                                      start1 = respstart, start2 = trstart, start3 = instart) {
  
  respstart[[1]][[1]][1] <- x
  
  model <- HMM_simulate(n = N, nstates = k, respstart = start1, trstart = start2, instart = start3)
  
  sim <- simulate(model)
  
  class(sim) <- "depmix"
  
  out <- fit(sim, emc = em.control(random.start = F))
  
  return(out)
})


stopCluster(clust)

exp(summary(estimates[[2]]))
