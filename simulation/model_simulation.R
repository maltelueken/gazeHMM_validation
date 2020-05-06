### Simulation from hidden Markov model ###

# Author: Malte Lüken
# Date: 06.05.2020

HMM_simulate <- function(n, nstates, respstart, trstart, instart) {
  
  require(depmixS4)
  source("algorithm/model_helper_functions.R")
  
  
  data <- data.frame(vel = numeric(n),
                     acc = numeric(n),
                     angle = seq(0, 2*pi, length.out = n))
  
  
  # Create response model
  
  resp <- list(list(altGamma(data$vel, pstart = respstart[["fix"]][["vel"]]),
                    altGamma(data$acc, pstart = respstart[["fix"]][["acc"]]),
                    unif(data$angle)))
  
  for (s in 2:nstates) {
    
    resp[[s]] <- list(altGamma(data$vel, pstart = respstart[[s]][["vel"]]),
                      altGamma(data$acc, pstart = respstart[[s]][["acc"]]),
                      vMF(data$angle, pstart = respstart[[s]][["angle"]]))
    
  }
  
  
  # Create state transition model
  
  trans <- list()
  
  for (s in 1:nstates) {
    
    trans[[s]] <- transInit(~ 1, nstates = nstates, data = data, pstart = c(trstart[s,]))
    
  }
  
  
  # Create initial state model
  
  init <- transInit(~ 1, ns = nstates, ps = instart, family = multinomial("identity"), data = data.frame(1))
  
  
  # Combine models
  
  model <- makeDepmix(response = resp, transition = trans, prior = init, homogeneous = F)
  
  return(model)
}
