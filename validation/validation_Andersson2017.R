### Validation Andersson et al. (2017) ###

# Author: Malte Lüken
# Date:17.05.2020

library(tidyverse)
library(depmixS4)
library(parallel)

source("algorithm/gazeHMM.R")
source("algorithm/model_helper_functions.R")
source("algorithm/plotting_functions.R")
source("algorithm/preprocessing.R")


# Function to fit 1-state HMM

onestate_HMM <- function(data, respstart, sf = c(10, 10)) {
  
  source("algorithm/model_helper_functions.R")
  
  if(missing(respstart)) {
    
    respstart <- list(fix = list(vel = c(1, 1), acc = c(1, 1)))
    
  }
  
  if(all(sf > 0)) {
    
    data$vel <- data$vel/sf[1]
    data$acc <- data$acc/sf[2]
    
  }
  
  velresp <- fit(altGamma(data$vel, pstart = respstart[["fix"]][["vel"]]))
  accresp <- fit(altGamma(data$acc, pstart = respstart[["fix"]][["acc"]]))
  angresp <- fit(unif(data$angle))
  
  ll <- sum(dens(velresp, log = T) + dens(accresp, log = T) + dens(angresp, log = T))
  
  n <- nrow(data)
  
  return(list(ll, n, response = list(velresp, accresp, angresp)))
}


# Load data

dir <- "~/Uni/Psychologie Master/Internship/Data/Test data/Andersson et al. (2017)/annotated_data/data used in the article"

A2017 <- list()

for(i in c("dots", "img", "video")) {
  
  filenames.MN <- list.files(path = paste(dir, "/", i, "/", sep = ""), pattern = "_MN")
  filenames.RA <- list.files(path = paste(dir, "/", i, "/", sep = ""), pattern = "_RA")
  
  data.MN <- lapply(paste(dir, "/", i, "/", filenames.MN, sep = ""), function(x) readMat(x)$ETdata[[1]])
  data.RA <- lapply(paste(dir, "/", i, "/", filenames.RA, sep = ""), function(x) readMat(x)$ETdata[[1]])
  
  A2017[[i]] <- lapply(1:length(data.MN), function(i, x, y) {

    print(length(x[[i]]) - length(y[[i]]))

    comb <- as.data.frame(x[[i]]) %>% left_join(as.data.frame(y[[i]]), by = names(.)[1:5])

    names(comb) <- c("t", "h_pupil", "v_pupil", "x", "y", "label_MN", "label_RA")

    return(comb)
  }, x = data.MN, y = data.RA)
}


# Classify data

res <- c(1024, 768)
dim <- c(380, 300)
dist <- 670
fr <- 500

ncores <- detectCores()
clust <- makeCluster(ncores)
clusterExport(clust, list("classify_gaze_data", "onestate_HMM", "preprocess", "res", "dim", "dist", "fr", "A2017"))

A2017.fit <- lapply(A2017, function(stim) parLapply(clust, stim, function(df) {
  
  if(all(is.nan(df$t) || df$t == 0)) {
    
    df$t <- seq(0, length(df$t)-1)/fr
      
  } else {
    
    df$t <- (df$t - df$t[1])/1e6
    
    df <- df[df$t >= 0,]
    
  }
  
  fit <- list()
  
  fit[["1"]] <- try(onestate_HMM(preprocess(x = df$x, y = df$y, t = df$t, unit = "px",
                                            res = res, dim = dim, dist = dist, fr = fr, blink = c(0, 0),
                                            sg.order = 3, sg.length = 5)))
  
  for (k in 2:5) {
    
    fit[[as.character(k)]] <- try(classify_gaze_data(x = df$x, y = df$y, t = df$t, unit = "px",
                                                     res = res, dim = dim, dist = dist, fr = fr, blink = c(0, 0),
                                                     nstates = k))
    
  }
  
  return(fit)
}))

stopCluster(clust)

save("A2017.fit", file = "validation/Andersson2017_fitted.Rdata")
