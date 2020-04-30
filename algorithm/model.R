### Classification by hidden Markov model ###

# Author: Malte Lüken
# Date: 25.05.2020


HMM_classify <- function(data, nstates, respstart, trstart, instart, 
                         sf = c(10, 10),
                         fit.control = em.control(maxit = 5000, random.start = F)) {
  
  require(depmixS4)
  source("algorithm/model_helper_functions.R")
  
  
  # Check if starting values for response model supplied, set to default if not
  
  if(missing(respstart)) {
    if(nstates == 1) {
      
      respstart <- list(fix = list(vel = c(1, 1), acc = c(1, 1)))
      
    } else if(nstates == 2) {
      
      respstart <- list(fix = list(vel = c(1, 1), acc = c(1, 1)), 
                        sac = list(vel = c(5, 5), acc = c(5, 5), angle = c(0, 10)))
      
    } else if(nstates == 3) {
      
      respstart <- list(fix = list(vel = c(1, 1), acc = c(1, 1)), 
                        sac = list(vel = c(5, 5), acc = c(5, 5), angle = c(0, 10)),
                        pso = list(vel = c(5, 5), acc = c(5, 5), angle = c(pi, 10)))
      
    } else if(nstates == 4) {
      
      respstart <- list(fix = list(vel = c(1, 1), acc = c(1, 1)), 
                        sac = list(vel = c(5, 5), acc = c(5, 5), angle = c(0, 10)),
                        pso = list(vel = c(5, 5), acc = c(5, 5), angle = c(pi, 10)),
                        sp = list(vel = c(2, 2), acc = c(2, 2), angle = c(0, 2)))
      
    } else if(nstates == 5) {
      
      respstart <- list(fix = list(vel = c(1, 1), acc = c(1, 1)), 
                        sac = list(vel = c(5, 5), acc = c(5, 5), angle = c(0, 10)),
                        pso = list(vel = c(5, 5), acc = c(5, 5), angle = c(pi, 10)),
                        sp = list(vel = c(2, 2), acc = c(2, 2), angle = c(0, 2)),
                        mic = list(vel = c(2, 2), acc = c(5, 5), angle = c(0, 10)))
      
    } 
  }
  
  
  # Check if starting values for transition model are supplied, set to default if not
  
  if(missing(trstart)) trstart <- matrix(1/nstates, nrow = nstates, ncol = nstates)
  
  
  # Check if starting values for initial state model are supplied, set to default if not
  
  if(missing(instart)) instart <- rep(1/nstates, nstates)
  
  
  # Downsample velocity and acceleration data
  
  if(all(sf > 0)) {
    
    data$vel <- data$vel/sf[1]
    data$acc <- data$acc/sf[2]
    
  }
  
  
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
  
  
  # Fit model
  
  time <- Sys.time()
  
  model.fit <- fit(model, emcontrol = fit.control)
  
  cat("Computation time: ", Sys.time() - time, "\n")
  
  
  return(model.fit)
}
