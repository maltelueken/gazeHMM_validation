### Hidden Markov model parameter recovery ###

# Author: Malte Lüken
# Date: 06.05.2020

library(depmixS4)
library(parallel)
source("simulation/model_simulation.R")
source("algorithm/model_helper_functions.R")


# Part 1 ------------------------------------------------------------------

# Global parameters

D <- 5
N <- 2500


# Transition probability interval

aij <- seq(0.01, 0.99, length.out = D)


# Response parameter intervals

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


# Fixed response parameter values

par.fix.true <- list(fix = list(vel = c(3, 0.35), acc = c(3, 0.15)), 
                     sac = list(vel = c(3, 15), acc = c(3, 3), angle = c(0, 1)),
                     pso = list(vel = c(3, 3), acc = c(2, 2), angle = c(pi, 1)),
                     sp = list(vel = c(3, 1.5), acc = c(3, 0.15), angle = c(0, 1)))


# Parameter indicators

par.states <- c(rep(1, 4), rep(2:4, each = 5))
par.names <- c(1, 2, 1, 2, rep(c(1, 2, 1, 2, 2), 3))
par.resps <- c(rep(1:2, each = 2), rep(c(1, 1, 2, 2, 3),  3))


# Prepare parallel processing

estimates <- list()

for (p in 1:sum(par.states <= 2)) {

  ncores <- detectCores()
  clust <- makeCluster(ncores)
  clusterExport(clust, list("HMM_simulate", "N", "par.fix.true", "par.states", "par.names", "par.resps", "p"))
  
  estimates[[p]] <- parLapply(clust, par.int.true[,p], function(x, samples = N, k = 2, par.state = par.states[p], 
                                                           par.name = par.names[p], par.resp = par.resps[p], 
                                                           trueResp = par.fix.true) {
    
    # Set true parameter values
    
    trueResp[[par.state]][[par.resp]][par.name] <- x
    
    trueIn <- rep(1/k, k)
    
    trueTr <- matrix(c(0.9, 0.1, 
                       0.1, 0.9), nrow = k, ncol = k)
    
    
    # Simulate data with model
    
    model <- HMM_simulate(n = N, nstates = k, trueresp = trueResp, truetr = trueTr, truein = trueIn)
    
    model.sim <- simulate(model)
    
    class(model.sim) <- "depmix"
    
    
    # Generate and set starting values
    
    respStart <- lapply(trueResp, function(x, beta = 1) {
      
      lapply(x, function(y) {
        
        vec <- log(c(rgamma(1, shape = 3, scale = y[1]), rgamma(1, shape = 3, scale = y[2])))
        # vec <- log(y)
        
        return(vec) 
      })
    })
    
    respStart[["fix"]][["angle"]] <- c(0, 2*pi)
    respStart[["sac"]][["angle"]][1] <- 0
    respStart[["pso"]][["angle"]][1] <- pi
    respStart[["sp"]][["angle"]][1] <- 0
    
    inStart <- trueIn
    trStart <- logit(matrix(1/k, nrow = k, ncol = k))
    
    start <- c(inStart, trStart, unlist(respStart)[1:12])
    
    model.start <- setpars(model.sim, start)
    
    model.fit <- try(fit(model.start, emc = em.control(maxit = 5000, random.start = F)))
    
    acc <- try(mean(model.sim@states == model.fit@posterior$state))
    
    output <- list(pars.true = getpars(model), pars.start = getpars(model.start), 
                   start = unlist(respStart)[1:12], pars.est = model.fit, accuracy = acc)
    
    return(output)
  })
  
  stopCluster(clust)
  
}

exp(summary(estimates[[2]]))

di <- lapply(list(u = list("1" = 1, "2" = 2), p = list(3, 4)), function(x) {lapply(x, function(y) {print(names(x))})})

ma <- logit(matrix(0.5, 2, 2))

exp(ma)/(1+exp(ma))
