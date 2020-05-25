### Classification by hidden Markov model ###

# Author: Malte Lüken
# Date: 25.05.2020


HMM_classify <- function(data, nstates, respstart, trstart, instart, 
                         sf = c(10, 10),
                         fit.control = em.control(maxit = 5000, random.start = F)) {
  
  require(depmixS4)
  source("algorithm/model_helper_functions.R")
  
  
  # Downsample velocity and acceleration data
  
  if(all(sf > 0)) {
    
    data$vel <- data$vel/sf[1]
    data$acc <- data$acc/sf[2]
    
  }
  
  
  # Create response model
  
  resp <- list(list(altGamma(data$vel, pstart = respstart[[1]][[1]]),
                    altGamma(data$acc, pstart = respstart[[1]][[2]]),
                    unif(data$angle)))
  
  for (s in 2:nstates) {
    
    resp[[s]] <- list(altGamma(data$vel, pstart = respstart[[s]][[1]]),
                      altGamma(data$acc, pstart = respstart[[s]][[2]]),
                      vMF(data$angle, pstart = respstart[[s]][[3]]))
    
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
  
  
  # Fit model
  
  time <- Sys.time()
  
  model.fit <- fit(model, emcontrol = fit.control)
  
  cat("Computation time: ", Sys.time() - time, "\n")
  
  
  return(model.fit)
}
