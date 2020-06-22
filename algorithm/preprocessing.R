### Preprocessing ###

# Author: Malte Lüken
# Date: 07.04.2020


preprocess <- function(x, y, t, unit = "px", res, dim, dist, fr, blink = NULL, b.win = 0.05,
                       sg.order, sg.length) {
  
  require(signal)
  source("algorithm/preprocessing_helper_functions.R")
  
  
  # Initialize sample metrics
  
  vel <- rep(NA, length(x))
  acc <- rep(NA, length(x))
  angle <- rep(NA, length(x))
  label <- rep(NA, length(x)) # for blinks
  
  
  # Check for NA and Inf and (0,0)
  
  valid <- ifelse(is.na(x) | is.na(y) | 
                    is.infinite(x) | is.infinite(y) |
                    (x == 0 & y == 0), F, T)
  
  
  # Initial label blinks
  
  inib <- rep(NA, length(x))
  
  if(is.numeric(blink) && length(blink) == 2) {
    
    inib[x == blink[1] & y == blink[2]] <- 0
    
  } else if(is.logical(blink) && length(blink) == length(t)) {
    
    inib[blink] <- 0
    
  }
  
  
  # Convert pixels to degrees of visual angle
  
  res.va <- c(px2va(res[1], res = res[1], dim = dim[1], dist = dist),
              px2va(res[2], res = res[2], dim = dim[2], dist = dist))
  
  if(unit == "px") {
    
    x.va <- px2va(x, res = res[1], dim = dim[1], dist = dist)
    y.va <- px2va(y, res = res[2], dim = dim[2], dist = dist)
    
  } else {
    
    x.va <- x
    y.va <- y
    
  }
  
  
  # Check if samples are on screen
  
  valid[valid] <- ifelse(x.va[valid] < -res.va[1] | x.va[valid] > res.va[1] | 
                    y.va[valid] < -res.va[2] | y.va[valid] > res.va[2], F, valid[valid])
  
  
  # Apply Butterworth filter

  # x.va[valid] <- as.numeric(signal::filter(butter(3, W = 0.5), x = x.va[valid]))
  # y.va[valid] <- as.numeric(signal::filter(butter(3, W = 0.5), x = y.va[valid]))
  
  
  # Calculate velocity and acceleration
  
  x.vel <- sgolayfilt(x.va[valid], m = 1, p = sg.order, n = sg.length)
  y.vel <- sgolayfilt(y.va[valid], m = 1, p = sg.order, n = sg.length)
  vel[valid] <- sqrt(x.vel^2 + y.vel^2) * fr
  
  x.acc <- sgolayfilt(x.va[valid], m = 2, p = sg.order, n = sg.length)
  y.acc <- sgolayfilt(y.va[valid], m = 2, p = sg.order, n = sg.length)
  acc[valid] <- sqrt(x.acc^2 + y.acc^2) * fr
  
  valid[valid] <- ifelse(is.na(vel[valid]) | is.na(acc[valid]), F, valid[valid])
  
  
  # Nudge zero velocities and accelerations
  
  nudge <- 0.01
  
  vel <- ifelse(vel == 0, nudge, vel)
  acc <- ifelse(acc == 0, nudge, acc)
  
  
  # Set +/- 50 ms around blinks to invalid
  
  for(i in 1:length(inib)) {
    if(!is.na(inib[i])) {
      
      valid[t > (t[i]-b.win) & t <= (t[i]+b.win)] <- F
      label[t > (t[i]-b.win) & t <= (t[i]+b.win)] <- 0
      
    }  
  }
  
  
  # Remove outliers
  
  max.vel <- 1e3
  max.acc <- 1e2
  
  valid[valid] <- ifelse(vel[valid] > max.vel | acc[valid] > max.acc, F, valid[valid])
  
  
  # Calculate angle
  
  angle[valid] <- calc_theta(x.va[valid], y.va[valid])
  
  valid[valid] <- ifelse(is.na(angle[valid]), F, valid[valid])
  
  
  # Print number of excluded frames
  
  warning(paste(sum(!valid), "samples were labeled as 'invalid' during preprocessing!"))
  
  
  # Set invalid samples to NA (if not already NA)
  
  vel[!valid] <- NA
  acc[!valid] <- NA
  angle[!valid] <- NA
  
  
  # Combine metrics in output data frame
  
  output <- data.frame(frame = 1:length(x), x = x.va, y = y.va, t = t, vel, acc, angle, label) # [valid,]
  
  return(output)
}
