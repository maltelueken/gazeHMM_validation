### Validation Andersson et al. (2017) ###

# Author: Malte Lüken
# Date:17.05.2020

library(tidyverse)
library(depmixS4)
library(parallel)
library(R.matlab)

source("algorithm/gazeHMM.R")
source("algorithm/model_helper_functions.R")
source("algorithm/plotting_functions.R")
source("algorithm/preprocessing.R")


# Function to fit 1-state HMM

onestate_HMM <- function(data, respstart, sf = c(10, 10), start.seed = NULL) {
  
  source("algorithm/model_helper_functions.R")
  
  if(missing(respstart)) {
    
    respstart <- list(fix = list(vel = c(1, 1), acc = c(1, 1), angle = c(0, 2*pi)))
    
    respstart <- lapply(1:length(respstart), function(x) {
      lapply(1:length(respstart[[x]]), function(y) {
        
        if(y < 3) {
          
          out <- sapply(respstart[[x]][[y]], gamma_start, 
                        seed = start.seed[[x]][[y]])
          
        } else {
          
          out <- c(respstart[[x]][[y]][1], gamma_start(respstart[[x]][[y]][2], 
                                                       seed = start.seed[[x]][[y]]))
          
        }
        
        return(out)
      })
    })
  }
  
  if(all(sf > 0)) {
    
    data$vel <- data$vel/sf[1]
    data$acc <- data$acc/sf[2]
    
  }
  
  velresp <- fit(altGamma(data$vel, pstart = respstart[[1]][[1]]))
  accresp <- fit(altGamma(data$acc, pstart = respstart[[1]][[2]]))
  angresp <- fit(unif(data$angle))
  
  ll <- sum(dens(velresp, log = T) + dens(accresp, log = T) + dens(angresp, log = T))
  
  n <- nrow(data)
  
  return(list(LL = ll, N = n, response = list(velresp, accresp, angresp), respstart = respstart))
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
clusterExport(clust, list("gazeHMM", "onestate_HMM", "preprocess", "res", "dim", "dist", "fr", "A2017"))

A2017.fit <- lapply(1:length(A2017[1:2]), function(stim) parLapply(clust, 1:length(A2017[[stim]][1:5]), function(subj) {
  
  data.seed <- stim*1e4 + subj*1e2
  
  df <- A2017[[stim]][[subj]]
  
  if(all(is.nan(df$t) || df$t == 0)) {
    
    df$t <- seq(0, length(df$t)-1)/fr
      
  } else {
    
    df$t <- (df$t - df$t[1])/1e6
    
    df <- df[df$t >= 0,]
    
  }
  
  fit <- list()
  
  start.seed <- lapply(1:3, function(y) return(as.integer(data.seed + 10 + y)))
  
  fit[["1"]] <- try(onestate_HMM(preprocess(x = df$x, y = df$y, t = df$t, unit = "px",
                                            res = res, dim = dim, dist = dist, fr = fr, blink = c(0, 0),
                                            sg.order = 3, sg.length = 5)))
  
  for (k in 4:4) {
    
    start.seed <- lapply(1:k, function(x) lapply(1:3, function(y) return(as.integer(data.seed + x*10 + y))))
    
    fit[[as.character(k)]] <- try(gazeHMM(x = df$x, y = df$y, t = df$t, unit = "px",
                                                     res = res, dim = dim, dist = dist, fr = fr, blink = c(0, 0),
                                                     nstates = k, start.seed = start.seed))
    
  }
  
  return(fit)
}))

stopCluster(clust)

save("A2017.fit", file = "validation/Andersson2017_fitted.Rdata")


# Compare models with different states

schwarz.weights <- function(bic, na.rm = T) {
  
  d.bic <- bic - min(bic, na.rm = na.rm) # eq 2
  
  exp(-0.5 * d.bic)/sum(exp(-0.5 * d.bic), na.rm = na.rm) # eq 4
  
}

A2017.bic <- lapply(A2017.fit, function(stim) {
  out <- lapply(stim, function(subj) {
    out <- lapply(1:length(subj), function(mod) {
      
      if(mod == 1) {
        
        bic <- try(-2*subj[[mod]][["LL"]] + 6*log(subj[[mod]][["N"]]))
        
      } else {
        
        bic <- try(BIC(subj[[mod]]$model))
        
      }
      
      return(ifelse(is.numeric(bic), bic, NA))
    })
    
    return(schwarz.weights(unlist(out)))
  })
  
  df <- as.data.frame(reduce(out, rbind))
  
  names(df) <- paste("model_", 1:length(stim[[1]]), sep = "")
  
  return(df)
})

lapply(A2017.fit, function(x) lapply(x, function(y) try(print(plot_samples_xy(y[["4"]], 0, 1)))))
