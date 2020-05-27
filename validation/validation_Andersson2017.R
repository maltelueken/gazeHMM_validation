### Validation Andersson et al. (2017) ###

# Author: Malte Lüken
# Date:17.05.2020

library(tidyverse)
library(depmixS4)
library(parallel)
library(R.matlab)
library(psych)

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

A2017.fit <- lapply(1:length(A2017), function(stim) parLapply(clust, 1:length(A2017[[stim]]), function(subj) {
  
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
                                            sg.order = 3, sg.length = 5), random.respstart = F))
  
  for (k in 2:5) {
    
    start.seed <- lapply(1:k, function(x) lapply(1:3, function(y) return(as.integer(data.seed + x*10 + y))))
    
    trstart <- matrix(0.1/(k-1), k, k)
    diag(trstart) <- 0.9
    
    fit[[as.character(k)]] <- try(gazeHMM(x = df$x, y = df$y, t = df$t, unit = "px",
                                          res = res, dim = dim, dist = dist, fr = fr, blink = c(0, 0),
                                          nstates = k, trstart = trstart, random.respstart = F, 
                                          start.seed = start.seed, min.sac = 0.01))
    
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

lapply(A2017.fit, function(x) lapply(x, function(y) try(print(plot_samples_time(y[["4"]])))))


# Calculate event descriptives for human coders

k <- 4

A2017.events <- lapply(A2017, function(stim) {
  lapply(1:k, function(e) {
    out <- lapply(stim, function(df) {
      
      if(all(is.nan(df$t) || df$t == 0)) {
        
        df$t <- seq(0, length(df$t)-1)/fr
        
      } else {
        
        df$t <- (df$t - df$t[1])/1e6
        
        df <- df[df$t >= 0,]
        
      }
      
      counter_MN <- 1
      counter_RA <- 1
      
      number_MN <- numeric(nrow(df))
      number_RA <- numeric(nrow(df))
      
      number_MN[1] <- 1
      number_RA[1] <- 1
      
      for (i in 2:nrow(df)) {
        if(df$label_MN[i] != df$label_MN[i-1]) counter_MN <- counter_MN + 1
        if(df$label_RA[i] != df$label_RA[i-1]) counter_RA <- counter_RA + 1
        
        number_MN[i] <- counter_MN
        number_RA[i] <- counter_RA
    
      }
      
      dur_MN <- numeric(max(number_MN))
      event_MN <- numeric(max(number_MN))
      
      for (n in unique(number_MN)) {
        if(n > 1) {
          dur_MN[n] <- max(df$t[number_MN == n], na.rm = T) - max(df$t[number_MN == (n-1)], na.rm = T)
        } else {
          dur_MN[n] <- max(df$t[number_MN == n], na.rm = T) - min(df$t[number_MN == n], na.rm = T)
        }
        
        event_MN[n] <- max(df$label_MN[number_MN == n], na.rm = T)
      }
      
      dur_RA <- numeric(max(number_RA))
      event_RA <- numeric(max(number_RA))
      
      for (n in unique(number_RA)) {
        if(n > 1) {
          dur_RA[n] <- max(df$t[number_RA == n], na.rm = T) - max(df$t[number_RA == (n-1)], na.rm = T)
        } else {
          dur_RA[n] <- max(df$t[number_RA == n], na.rm = T) - min(df$t[number_RA == n], na.rm = T)
        }
        
        event_RA[n] <- max(df$label_RA[number_RA == n], na.rm = T)
      }
      
      return(list(dur_MN[event_MN == e], dur_RA[event_RA == e]))
    })
    
    df <- reduce(out, cbind)
    
    dur_MN <- unlist(df[1,])
    dur_RA <- unlist(df[2,])
    
    return(data.frame(MN = c(mean(dur_MN), sd(dur_MN), length(dur_MN)),
                      RA = c(mean(dur_RA), sd(dur_RA), length(dur_RA))))
  })
})


# Calculate RMSD between event distribution descriptives

ref.max <- list(dot = list(fix = c(0.380, 0.333, 165),
                           sac = c(0.060, 0.026, 93),
                           pso = c(0.024, 0.012, 33)),
                img = list(fix = c(0.399, 0.559, 827),
                           sac = c(0.062, 0.037, 787),
                           pso = c(0.028, 0.013, 319)),
                vid = list(fix = c(0.554, 0.825, 227),
                           sac = c(0.055, 0.053, 1104),
                           pso = c(0.028, 0.013, 97)))

ref.min <- list(dot = list(fix = c(0.060, 0.030, 2),
                           sac = c(0.013, 0.005, 10),
                           pso = c(0.015, 0.005, 17)),
                img = list(fix = c(0.144, 0.136, 251),
                           sac = c(0.017, 0.010, 258),
                           pso = c(0.021, 0.009, 237)),
                vid = list(fix = c(0.202, 0.189, 48),
                           sac = c(0.018, 0.010, 41),
                           pso = c(0.017, 0.008, 78)))

A2017.rmsd <- lapply(1:length(A2017.fit), function(x) { 
  out <- lapply(1:k, function(y) {
    out <- lapply(A2017.fit[[x]], function(z) {
      
      if(class(z[[k]]) == "gazeHMM") z[[k]]$events[[y]]
      
    })
    
    df <- try(reduce(out, rbind))
    
    alg <- try(c(mean(df$dur, na.rm = T), sd(df$dur, na.rm = T), nrow(df)))
    
    if(y == 4) {
      M <- try(cbind(A2017.events[[x]][[y]], alg))
    } else {
      M <- try(cbind(A2017.events[[x]][[y]], alg, ref.max[[x]][[y]], ref.min[[x]][[y]]))
    }
      
    M.norm <- try(apply(M, 1, function(x) {(x - min(x))/(max(x) - min(x))}))
    
    return(try(sum(sqrt((M.norm[3,] - colSums(M.norm[1:2,]/2))^2))))
  })
})


# Calculate agreement with human raters

k <- 3

A2017.acc <- lapply(1:length(A2017.fit), function(x) { 
  out <- lapply(1:length(A2017.fit[[x]]), function(y) {
  
    if(class(A2017.fit[[x]][[y]][[as.character(k)]]) == "gazeHMM") {
      cod.MN <- as.data.frame(cbind(A2017.fit[[x]][[y]][[as.character(k)]]$samples$label, A2017[[x]][[y]]$label_MN))
      cod.RA <- as.data.frame(cbind(A2017.fit[[x]][[y]][[as.character(k)]]$samples$label, A2017[[x]][[y]]$label_RA))
      
      return(rbind(cod.MN, cod.RA))
    }
  })
  
  df <- reduce(out, rbind)
  
  df[df[,1] == 0,1] <- 5
  
  events <- lapply(1:k, function(e) {

    alg <- df[,1] == e
    cod <- df[,2] == e

    mat <- matrix(c(sum(alg & cod, na.rm = T), sum(alg & !cod, na.rm = T), 
                    sum(!alg & cod, na.rm = T), sum(!alg & !cod, na.rm = T)), 2, 2)
    
    over <- mat[2,1]/sum(mat[,1])
    under <- mat[1,2]/sum(mat[,2])
    
    return(list(kappa = cohen.kappa(mat)$kappa, 
                over = over,
                under = under))
  })
  
  ratio <- mean(df[,1] != df[,2], na.rm = T)
  
  conmat <- table(df[,1], df[,2])
  
  return(list(events = events, ratio =  ratio, conmat = conmat))
})

