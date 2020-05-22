### gazeHMM - full algorithm ###

# Author: Malte Lüken
# Date: 24.04.2020


gazeHMM <- function(x, y, t, unit = "px", res, dim, dist, fr, blink = NULL,
                               sg.order = 3, sg.length = 5, 
                               nstates, respstart, trstart, instart, sf = c(10, 10),
                               fit.control = em.control(maxit = 5000, random.start = F),
                               min.sac = 0.01) {
  
  # Store settings
  
  settings <- list(x, y, t, unit, res, dim, dist, fr, blink, sg.order, sg.length, nstates, sf, fit.control, min.sac)
  
  
  # Source algorithm parts
  
  source("algorithm/preprocessing.R")
  source("algorithm/model.R")
  source("algorithm/postprocessing.R")
  
  
  # Validate arguments
  
  if(!is.numeric(x)) stop("'x' must be numeric")
  if(!is.numeric(y)) stop("'y' must be numeric")
  if(!is.numeric(t)) stop("'t' must be numeric")
  
  if(!all.equal(length(x), length(y), length(t))) stop("'x', 'y', and 't' must have the same length")
  
  if(t < 0) stop("'t' must be zero or positive")
  
  if(!(unit %in% c("px", "va"))) stop("'unit' must be either 'px' or 'va'")
  
  if(!is.numeric(res)) stop("'res' must be numeric")
  if(length(res) != 2 || !is.vector(res)) stop("'res' must be a vector of length 2")
  
  if(!is.numeric(dim)) stop("'dim' must be numeric")
  if(length(dim) != 2 || !is.vector(dim)) stop("'dim' must be a vector of length 2")
  
  if(!is.numeric(dist)) stop("'dist' must be numeric")
  if(length(dist) != 1) stop("'dist' must be a single value")
  if(dist <= 0) stop("'dist' must be greater than zero")
  
  if(!is.numeric(fr)) stop("'fr' must be numeric")
  if(length(fr) != 1) stop("'fr' must be a single value")
  if(dist <= 0) stop("'fr' must be positive and greater than 0")
  
  if(!(is.numeric(blink) && is.vector(blink) && length(blink) == 2) && 
     !(is.logical(blink) && is.vector(blink) && length(blink) == length(t))) stop(
       "'blink' must be either a numeric vector of length 2 or a logical vector with the same length as 't'"
       )
  
  if(!is.integer(nstates)) stop("'nstates must be integer'")
  if(length(nstates) != 1) stop("'nstates' must be a single value")
  if(nstates < 2 || nstates > 5) stop("'nstates' must have a value between 2 and 5")
  
  if(!missing(respstart)) {
    if(!is.list(respstart)) stop("'respstart' must be a list")
    if(length(respstart) != nstates) stop("'respstart' must have length equal to 'nstates'")
    for(i in 1:length(respstart)) {
      if(length(respstart[[i]]) != 3) stop(
        paste("'respstart[['", i, "]] must contain a sublist for each response variable", sep = "")
      )
      for (j in 1:length(respstart[[i]])) {
        if(length(respstart[[i]][[i]] != 2)) stop(
          paste("'respstart[['", i, "]][[", j, "]] must have two parameter starting values", sep = "")
        )
        if(!is.numeric(respstart[[i]][[i]])) stop(
          paste("'respstart[['", i, "]][[", j, "]] must be numeric", sep = "")
        )
      }
    }
  }
  
  if(!missing(trstart)) {  
    if(!is.matrix(trstart)) stop("'trstart' must be a matrix")
    if(!is.numeric(trstart)) stop("'trstart' must be numeric")
    if(dim(trstart) != c(nstates, nstates)) stop("'trstart' must have dimensions equal to 'nstates'")
  }
  
  if(!missing(instart)) {
    if(!is.vector(instart)) stop("'instart' must be a vector")
    if(!is.numeric(instart)) stop("'instart' must be numeric")
    if(length(instart) != nstates) stop("'instart' must have length equal to 'nstates'")
  }
    
  if(!is.numeric(sf)) stop("'sf' must be numeric")
  if(!is.vector(sf)) stop("'sf' must be a vector")
  if(length(sf) != 2) stop("'sf' must have length 2")
  if(any(sf <= 0)) stop("'sf' must be greater than zero")
  
  if(!is.numeric(min.sac)) stop("'min.sac' must be numeric")
  if(length(min.sac) != 1) stop("'min.sac' must be a single value")
  if(min.sac <= 0) stop("'min.sac' must be greater than zero")
  
  
  # Preprocessing
  
  prep <- preprocess(x, y, t, unit, res, dim, dist, fr, blink, sg.order, sg.length)
  
  
  # Model classification
  
  model.fit <- HMM_classify(prep, nstates, respstart, trstart, instart, sf, fit.control)
  
  
  # Postprocessing
  
  post <- postprocess(prep, model.fit@posterior$state, min.sac)
  
  
  # Output: Data frame with samples, data frames with events, fitted model
  
  output <- list(samples = post$samples, events = post$events, model = model.fit, settings = settings)
  
  class(output) <- "gazeHMM"
  
  return(output)
}
