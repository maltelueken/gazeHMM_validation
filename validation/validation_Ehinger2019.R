### Validation Ehinger et al. (2019) ###

# Author: Malte Lüken
# Date: 18.05.2020

library(tidyverse)
library(depmixS4)
library(parallel)
library(FDBeye)
library(eyelinker)

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


# Load Eyelink data

dir <- "~/Uni/Psychologie Master/Internship/Data/Test data/Ehinger et al. (2019)/"


# # Convert .EDF files into .asc files
# 
# setwd("~/Uni/Psychologie Master/Internship/Data/Test data/Ehinger et al. (2019)/")
# 
# Sys.getenv("C:/Programme (x86)/SR Research/EyeLink/EDF_Access_API/Example", names = T)
# 
# edf.filenames <- list.files(path = getwd(), pattern = ".EDF")
# 
# converted <- edf2asc(edf.filenames)

  
filenames.el <- list.files(path = dir, pattern = ".asc")
  
E2019 <- lapply(paste(dir, "/", filenames.el, sep = ""), function(x) { 
  
  data <- read.asc(x) 
  
  task4 <- list()
  task5 <- list()
  
  for (b in 1:6) {
    
    task4.start <- data$msg$time[str_detect(data$msg$text, "MICROSACC start") & data$msg$block == b][2]
    task4.stop <- data$msg$time[str_detect(data$msg$text, "MICROSACC stop") & data$msg$block == b][2]
    
    task4[[b]] <- data$raw[data$raw$time > task4.start & data$raw$time < task4.stop & data$raw$block == b,]
    
    task5.start <- data$msg$time[str_detect(data$msg$text, "BLINK start") & data$msg$block == b][2]
    task5.stop <- data$msg$time[str_detect(data$msg$text, "BLINK stop") & data$msg$block == b][2]
    
    task5[[b]] <- data$raw[data$raw$time > task5.start & data$raw$time < task5.stop & data$raw$block == b,]

  } 
  
  return(list(task4, task5))
})
  

# Classify data

res <- c(1920, 1080)
dim <- c(531, 299)
dist <- 600
fr <- 500

ncores <- detectCores()
clust <- makeCluster(ncores)
clusterExport(clust, list("classify_gaze_data", "onestate_HMM", "preprocess", "res", "dim", "dist", "fr", "E2019"))

E2019.fit <- parLapply(clust, E2019, function(subj) lapply(subj, function(task) lapply(task, function(df) {
  
  df$time <- (df$time - df$time[1])/1e3

  
  fit <- list()
  
  fit[["1"]] <- try(onestate_HMM(preprocess(x = df$xp, y = df$yp, t = df$time, unit = "px",
                                            res = res, dim = dim, dist = dist, fr = fr, blink = is.na(df$ps),
                                            sg.order = 3, sg.length = 5)))
  
  for (k in 2:5) {

    fit[[as.character(k)]] <- try(classify_gaze_data(x = df$xp, y = df$yp, t = df$time, unit = "px",
                                                 res = res, dim = dim, dist = dist, fr = fr, blink = is.na(df$ps),
                                                 nstates = k))

  }
  
  return(fit)
})))

stopCluster(clust)

save("E2019.fit", file = "validation/Ehinger2019_fitted.Rdata")