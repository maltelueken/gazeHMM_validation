### Validation Andersson et al. (2017) ###

# Author: Malte Lüken
# Date: 17.05.2020

library(tidyverse)
library(depmixS4)
library(parallel)
library(R.matlab)
library(psych)
library(signal)
library(here)

source("algorithm/gazeHMM.R")
source("algorithm/model_helper_functions.R")
source("algorithm/plotting_functions.R")
source("algorithm/preprocessing.R")


# Function to fit 1-state HMM

onestate_HMM <- function(data, respstart, sf = c(10, 10), random.respstart = T, start.seed = NULL) {
  
  source("algorithm/model_helper_functions.R")
  
  if(missing(respstart)) {
    
    respstart <- list(fix = list(vel = c(1, 1), acc = c(1, 1), angle = c(0, 2*pi)))
    
    if(random.respstart) {
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


# Load raw data

A2017 <- list()

for(i in c("dots", "img", "video")) {
  
  filenames.MN <- list.files(path = paste(here("validation/data"), "/", i, "/", sep = ""), pattern = "_MN")
  filenames.RA <- list.files(path = paste(here("validation/data"), "/", i, "/", sep = ""), pattern = "_RA")
  
  data.MN <- lapply(paste(here("validation/data"), "/", i, "/", filenames.MN, sep = ""), function(x) readMat(x)$ETdata[[1]])
  data.RA <- lapply(paste(here("validation/data"), "/", i, "/", filenames.RA, sep = ""), function(x) readMat(x)$ETdata[[1]])
  
  A2017[[i]] <- lapply(1:length(data.MN), function(i, x, y) {

    print(length(x[[i]]) - length(y[[i]]))

    comb <- as.data.frame(x[[i]]) %>% left_join(as.data.frame(y[[i]]), by = names(.)[1:5])

    names(comb) <- c("t", "h_pupil", "v_pupil", "x", "y", "label_MN", "label_RA")

    return(comb)
  }, x = data.MN, y = data.RA)
}


# Load algorithm classification data

events <- readMat(here("validation/data/20150807.mat"))$topResult[1:3,,]

names(events) <- c("img", "dots", "video")


# Remove two data sets from dots condition

A2017[["dots"]] <- A2017[["dots"]][-c(10, 11)]
events[["dots"]] <- array(events[["dots"]][,,-c(10, 11)], c(13, 1, 9))

alg.names <- names(events[[1]][,1,1])


# Merge raw data with algorithm classification

A2017 <- lapply(names(A2017), function(x) {
  
  MN <- lapply(A2017[[x]], function(y) {y$label_MN})
  RA <- lapply(A2017[[x]], function(y) {y$label_RA})
  
  MN.ev <- lapply(1:length(A2017[[x]]), function(y) {
    as.data.frame(events[[x]][,,y] %>% reduce(cbind))[,1]
  })
  RA.ev <- lapply(1:length(A2017[[x]]), function(y) {
    as.data.frame(events[[x]][,,y] %>% reduce(cbind))[,2]
  })
  
  comp <- list(MN, MN.ev, RA, RA.ev)
  
  same <- lapply(1:length(comp[[1]]), function(y) {
    
    all.true <- sapply(1:length(comp[[2]]), function(z) {
      
      all(comp[[1]][[y]] == comp[[2]][[z]])

    })
    
    sum.false <- sapply(1:length(comp[[2]]), function(z) {
      
      sum.false <- sum(comp[[1]][[y]] != comp[[2]][[z]])

    })
    
    if(all.true > 0) {
      which(all.true == TRUE)
    } else {
      which.min(sum.false)
    }
  })
  
  lapply(1:length(A2017[[x]]), function(y) {
    
    toMerge <- as.data.frame(events[[x]][,,same[[y]]] %>% reduce(cbind))
    
    names(toMerge) <- alg.names
    
    out <- cbind(A2017[[x]][[y]][1:length(comp[[2]][[same[[y]]]]),], toMerge)
    
    print(all(out[,6] == out[,8]))
    
    if(all(is.nan(out$t) || out$t == 0)) {
      
      out$t <- seq(0, length(out$t)-1)/500
      
    } else {
      
      out$t <- (out$t - out$t[1])/1e6
      
      out <- out[out$t >= 0,]
      
    }
    
    return(out[-c(6, 7)])
  })
})

save(A2017, file = here("validation/Andersson2017_raw.Rdata"))


# Classify data

res <- c(1024, 768)
dim <- c(380, 300)
dist <- 670
fr <- 500

ncores <- detectCores()
clust <- makeCluster(ncores)
clusterExport(clust, list("gazeHMM", "onestate_HMM", "preprocess", "butter", "filter", "res", "dim", "dist", "fr", "A2017"))

A2017.fit <- lapply(1:length(A2017), function(stim) parLapply(clust, 1:length(A2017[[stim]]), function(subj) {
  
  data.seed <- stim*1e4 + subj*1e2
  
  df <- A2017[[stim]][[subj]]
  
  fit <- list()
  
  start.seed <- lapply(1:3, function(y) return(as.integer(data.seed + 10 + y)))
  
  fit[["1"]] <- try(onestate_HMM(preprocess(x = df$x, y = df$y, t = df$t, unit = "px",
                                            res = res, dim = dim, dist = dist, fr = fr, blink = c(0, 0),
                                            sg.order = 3, sg.length = 5), random.respstart = F))
  
  for (k in 2:5) {
    
    start.seed <- lapply(1:k, function(x) lapply(1:3, function(y) return(as.integer(data.seed + x*10 + y))))
    
    trstart <- matrix(0.1/(k-1), k, k)
    diag(trstart) <- 0.9
    
    fit[[as.character(k)]] <- try(gazeHMM(x = df$x, y = df$y, t = df$t, unit = "px",
                                          res = res, dim = dim, dist = dist, fr = fr, blink = c(0, 0),
                                          nstates = k, sf = c(100, 100), trstart = trstart, random.respstart = F, 
                                          start.seed = start.seed, min.sac = 0.01))
    
  }
  
  return(fit)
}))

stopCluster(clust)

save("A2017.fit", file = here("validation/Andersson2017_fitted.Rdata"))
