### Postprocessing algorithm ###

# Author: Malte Lüken
# Date: 20.04.2020


postprocess <- function(data, state, min.sac) { # data frame with dependent variables, vector with sample states, minimum saccade duration
  
  source("algorithm/postprocessing_helper_functions.R")
  
  
  # Assign labels to events
  
  FIX <- 1
  SAC <- 2
  PSO <- 3
  SP <- 4
  
  
  # Mark invalid samples
  
  valid <- !is.na(data$vel) & !is.na(data$acc) & !is.na(data$angle)
  
  
  # Relabel samples
  
  df <- reclassify(state, data$t, min.sac)
  
  label <- df$label
  number <- df$number
  
  warning(paste(sum(label != state), " samples were relabeled during postprocessing!"))
  
  
  # Initialize output event metrics 
  
  dur <- numeric(max(number))
  event <- integer(max(number))
  
  x <- numeric(max(number))
  y <- numeric(max(number))
  
  start.x <- numeric(max(number))
  end.x <- numeric(max(number))
  start.y <- numeric(max(number))
  end.y <- numeric(max(number))
  
  amp <- numeric(max(number))
  
  max.vel <- numeric(max(number))
  avg.vel <- numeric(max(number))
  max.acc <- numeric(max(number))
  avg.acc <- numeric(max(number))
  
  dir <- numeric(max(number))
  
  
  # Calculate output event metrics
  
  for (n in unique(number)) {
    
    # Duration
    
    if(n > 1) {
      dur[n] <- max(data$t[number == n], na.rm = T) - max(data$t[number == (n-1)], na.rm = T)
    } else {
      dur[n] <- max(data$t[number == n], na.rm = T) - min(data$t[number == n], na.rm = T)
    }
    event[n] <- max(label[number == n])
    
    if(event[n] == FIX) {
      
      # Fixation position

      x[n] <- ifelse(sum(valid[number == n]) == 0, NA, mean(data$x[number == n & valid], trim = 0.2, na.rm = T))
      y[n] <- ifelse(sum(valid[number == n]) == 0, NA, mean(data$y[number == n & valid], trim = 0.2, na.rm = T))
      
    } else {
      
      # Start and end position
      
      start.x[n] <- ifelse(sum(valid[number == n]) == 0, NA, data$x[number == n & valid][1])
      end.x[n] <- ifelse(sum(valid[number == n]) == 0, NA, data$x[number == n & valid][length(data$x[number == n & valid])])
      start.y[n] <- ifelse(sum(valid[number == n]) == 0, NA, data$y[number == n & valid][1])
      end.y[n] <- ifelse(sum(valid[number == n]) == 0, NA, data$y[number == n & valid][length(data$y[number == n & valid])])
      
      dx <- end.x[n] - start.x[n]
      dy <- end.y[n] - start.y[n]
      
      
      # Amplitude
      
      amp[n] <- ifelse(sum(valid[number == n]) == 0, NA, sqrt(dx^2 + dy^2))
      
      
      # Velocity and acceleration
      
      max.vel[n] <- ifelse(sum(valid[number == n]) == 0, NA, max(data$vel[number == n], na.rm = T))
      avg.vel[n] <- ifelse(sum(valid[number == n]) == 0, NA, mean(data$vel[number == n], na.rm = T))
      max.acc[n] <- ifelse(sum(valid[number == n]) == 0, NA, max(data$acc[number == n], na.rm = T))
      avg.acc[n] <- ifelse(sum(valid[number == n]) == 0, NA, mean(data$acc[number == n], na.rm = T))
      
      
      # Direction
      
      dir[n] <- ifelse(sum(valid[number == n]) == 0, NA, atan2(dy, dx))
      
    }
  }
  
  data$state <- state
  data$label[is.na(data$label)] <- label[is.na(data$label)] # only labels for valid samples
  data$event <- number
  
  
  # Prepare output
  
  fixations <- data.frame(cbind(x, y, dur)[event == FIX,])
  saccades <- data.frame(cbind(dur, amp, max.vel, max.acc, avg.vel, avg.acc, dir, start.x, start.y, end.x, end.y)[event == SAC,])
  psos <- data.frame(cbind(dur, amp, max.vel, max.acc, avg.vel, avg.acc, start.x, start.y, end.x, end.y)[event == PSO,])
  sps <- data.frame(cbind(dur, amp, max.vel, max.acc, avg.vel, avg.acc, dir, start.x, start.y, end.x, end.y)[event == SP,])
  
  output <- list(samples = as.data.frame(data), 
                 events = list(fixations = fixations, saccades = saccades, PSOs = psos, smooth.pursuits = sps))
  
  return(output)
}
