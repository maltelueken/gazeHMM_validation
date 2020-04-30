### Postprocessing helper functions ###

# Author: Malte Lüken
# Date: 23.04.2020


# Function to retrieve last (different) state

last_state <- function(state, s) {
  if(state[s-1] != state[s] || (s-1) == 1) {
    
    return(state[s-1])
  } else {
    
    return(last_state(state, s-1))
  }
}


# Function to relabel samples

reclassify <- function(state, t) {
  
  # Assign labels to events
  
  FIX <- 1
  SAC <- 2
  PSO <- 3
  SP <- 4
  
  
  # Relabel single sample fixations and smooth pursuit
  
  number <- integer(length(state))
  
  number[1] <- 1
  
  counter <- 1
  
  N <- length(state)
  
  for (s in 2:(N-1)) {
    
    if(state[s] %in% c(FIX, SP) && state[s-1] != state[s] && state[s+1] != state[s]) {
      
      state[s] <- state[s-1]
      
    } 
    
    
    # Count consecutive samples with same label (events)
    
    if(state[s] != state[s-1]) counter <- counter + 1
    
    number[s] <- counter
  }
  
  number[length(number)] <- max(number)
  
  
  # Relabel saccades with duration below threshold
  
  dur <- numeric(max(number))
  event <- integer(max(number))
  
  for (n in unique(number)) {
    
    if(n > 1) {
      dur[n] <- max(t[number == n], na.rm = T) - max(t[number == (n-1)], na.rm = T)
    } else {
      dur[n] <- max(t[number == n], na.rm = T) - min(t[number == n], na.rm = T)
    }
    event[n] <- max(state[number == n])
    
    if(event[n] == SAC && dur[n] < 0.01 && n > 1) {
      
      state[number == n] <- last_state(event, n)
      
    }
  }
  
  
  # Relabel PSOs that do not occur immediately after saccade
  
  for (s in 2:(N-1)) {
    if(state[s] == PSO && state[s+1] == SAC) {
      
      state[s] <- state[s+1]
      
    } else if (state[s] == PSO && state[s-1] %in% c(FIX, SP) && state[s+1] != SAC) {
      
      state[s] <- state[s-1]
      
    }
  }
  
  
  # Check whether all relabeling conditions are met
  
  noSingleSamples <- T
  allPSOsafterSaccades <- T
  minDuration <- T
  
  number <- integer(length(state))
  number[1] <- 1
  counter <- 1
  
  for (s in 2:(N-1)) {
    if(state[s] %in% c(FIX, SP ) && state[s-1] != state[s] && state[s+1] != state[s]) {
      
      noSingleSamples <- F
      
    } else if(state[s] == PSO && state[s+1] == SAC || state[s] == PSO && state[s-1] %in% c(FIX, SP) && state[s+1] != SAC) {
      
      allPSOsafterSaccades <- F
      
    }
    
    if(state[s] != state[s-1]) counter <- counter + 1
    
    number[s] <- counter
  }
  
  number[length(number)] <- max(number)
  
  dur <- numeric(max(number))
  event <- integer(max(number))
  
  for (n in unique(number)) {
    
    if(n > 1) {
      dur[n] <- max(t[number == n], na.rm = T) - max(t[number == (n-1)], na.rm = T)
    } else {
      dur[n] <- max(t[number == n], na.rm = T) - min(t[number == n], na.rm = T)
    }
    event[n] <- max(state[number == n])
    
    if(event[n] == SAC && dur[n] < 0.01 && n > 1) {
      
      minDuration <- F
      
    }
  }
  
  
  # Return output if conditions are met, otherwise repeat relabeling
  
  output <- data.frame(label = state, number = number)
  
  if(all(noSingleSamples, allPSOsafterSaccades, minDuration)) {
    
    return(output)
    
  } else {
    
    return(reclassify(state, t))
    
  }
}