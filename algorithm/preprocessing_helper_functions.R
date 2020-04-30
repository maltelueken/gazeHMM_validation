### Preprocessing helper functions ###

# Author: Malte Lüken
# Date: 01.04.2020


# Function to convert gaze coordinates from pixels into degrees of visual angle

px2va <- function(x, dim, res, dist) { # x- or y-pos in px, dim in mm, res in px, dist in mm
  
  x_centered <- x - res/2 # transform to centered coordinate system (0,0)
  
  radian <- atan(x_centered/(dist*(res/dim))) # visual angle in radians
  
  degree <- radian*(180/pi) # visual angle in degrees
  
  return(degree) 
}
  

# Function to compute the relative angle between subsequent samples 

calc_theta <- function(x, y) { # x-pos, y-pos
  
  require(dplyr)
  
  diff_x <- x - lag(x) # x_t - x_t-1
  diff_y <- y - lag(y) # y_t - y_t-1
  
  angle <- atan2(diff_y, diff_x) # absolute angle of vector xy_t - xy_t-1
  theta <- lead(angle) - angle # relative angle of two vectors at sample t
  theta.mirrored <- ifelse(theta < 0, theta + 2*pi, theta) # mirror negative relative angles to positive range
  
  return(theta.mirrored)
}