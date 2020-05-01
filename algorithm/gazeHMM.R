### gazeHMM - full algorithm ###

# Author: Malte Lüken
# Date: 24.04.2020


classify_gaze_data <- function(x, y, t, unit = "px", res, dim, dist, fr, blink = NULL,
                               sg.order = 3, sg.length = 5, 
                               nstates, respstart, trstart, instart, sf = c(10, 10),
                               fit.control = em.control(maxit = 5000, random.start = F),
                               min.sac = 0.01) {
  
  source("algorithm/preprocessing.R")
  source("algorithm/model.R")
  source("algorithm/postprocessing.R")
  
  
  # Preprocessing
  
  prep <- preprocess(x, y, t, unit, res, dim, dist, fr, blink, sg.order, sg.length)
  
  
  # Model classification
  
  model.fit <- HMM_classify(prep, nstates, respstart, trstart, instart, sf, fit.control)
  
  
  # Postprocessing
  
  post <- postprocess(prep, model.fit@posterior$state, min.sac)
  
  
  # Output: Data frame with samples, data frames with events, fitted model
  
  output <- list(samples = post$samples, events = post$events, model = model.fit)
  
  return(output)
}
