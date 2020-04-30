### gazeHMM - full algorithm ###

# Author: Malte Lüken
# Date: 24.04.2020


classify_gaze_data <- function(x, y, t, unit = "px", res, dim, dist, fr, blink = NULL,
                               nstates, respstart, trstart, instart, sf = c(10, 10),
                               fit.control = em.control(maxit = 5000, random.start = F)) {
  
  source("algorithm/preprocessing.R")
  source("algorithm/model.R")
  source("algorithm/postprocessing.R")
  
  
  # Preprocessing
  
  prep <- preprocess(x, y, t, unit, res, dim, dist, fr, blink)
  
  
  # Model classification
  
  model.fit <- HMM_classify(prep, nstates, respstart, trstart, instart, sf, fit.control)
  
  
  # Postprocessing
  
  post <- postprocess(prep, model.fit@posterior$state)
  
  
  # Output: Data frame with samples, data frames with events, fitted model
  
  output <- list(samples = post$samples, events = post$events, model = model.fit)
  
  return(output)
}
