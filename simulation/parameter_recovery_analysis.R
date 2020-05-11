### Parameter recovery analysis ###

# Author: Malte Lüken
# Date: 11.05.2020


library(tidyverse)
library(depmixS4)


# Load data

for (part in 1:4) {
  
  load(file = paste("simulation/part", part, ".Rdata", sep = ""))
  
}

linkinv <- function(eta,base) {
  linv <- function(eta,base) {
    pp <- numeric(length(eta))
    if(any(is.infinite(eta)) || any(eta > log(.Machine$double.xmax)) || any(eta < log(.Machine$double.xmin))) {
      pp[which(is.infinite(eta))] <- 1
      pp[which(eta > log(.Machine$double.xmax))] <- 1 # change this to something better!
    } else {
      expb <- exp(eta)
      sumb <- sum(expb)
      pp[base] <- 1/sumb
      pp[-base] <- expb[-base]/sumb
    }
    return(pp)
  }
  if(is.matrix(eta)) {
    if(ncol(eta)==1) {
      pp <- as.matrix(apply(eta,1,linv,base=base)) # fixes problem with column matrix eta
    } else pp <- t(apply(eta,1,linv,base=base)) 	
  } else {
    pp <- linv(eta,base)
  }
  return(pp)
}


# Function to transform parameters to normal scale

backtrans <- function(x) {
  
  out <- x
  
  nms <- names(x)
  
  out[str_detect(nms, "(Intercept)")] <- linkinv(x[str_detect(nms, "(Intercept)")], base = 1)
  out[nms %in% c("shape", "scale", "kappa")] <- exp(x[nms %in% c("shape", "scale", "kappa")])
  
  return(out)
}


# Calculate MSE

mse <- list()

for (part in 1:4) {
  
  mse[[part]] <- lapply(get(paste("estimates.", part, sep = "")), function(x) {
    lapply(x, function(y) {
      sqerr <- lapply(y, function(z) {
        
        err <- try((backtrans(z$pars.true) - backtrans(z$pars.est))^2)
        
        if(is.numeric(err)) {
          
          out <- err
          
        } else {
          
          out <- rep(NA, length(z$pars.true))
          
        }
        
        names(out) <- names(z$pars.true)
        
        return(out)
      })
      
      rows <- length(sqerr)

      nms <- names(sqerr[[1]])

      pars <- matrix(unlist(sqerr), nrow = rows, byrow = T)

      msqerr <- apply(pars, 2, mean, na.rm = T)

      names(msqerr) <- nms

      return(msqerr)
    })
  })
}


# Display MSE in data frame

mse.data <- lapply(mse, function(x) lapply(x, as.data.frame))
mse.data <- lapply(mse.data, function(x) lapply(x, function(y) {as.data.frame(t(as.matrix(y)))}))


# Calculate linear regression weights

regw <- list()

for (part in 1:1) {
  
  regw[[part]] <- lapply(get(paste("estimates.", part, sep = "")), function(x) {
    
    index <- 1:length(x)
    
    #print(index)
    
    lapply(index, function(i, y) {
      
      print(i)
      
      lapply(y[[i]], function(z) {
        
        nms <- names(z$pars.true)
        
        pars.tr <- logical(length(z$pars.true))
        pars.resp <- logical(length(z$pars.true))
        
        pars.tr[str_detect(nms, "(Intercept)")] <- T
        pars.resp[nms %in% c("shape", "scale", "kappa")] <- T
        
        #print(sum(pars.varied))
        
        if(i == 1) {
          
          intpar <- try(cbind(backtrans(z$pars.true[pars.tr]), backtrans(z$pars.est[pars.tr])))
          
        } else {
          
          intpar <- try(c(backtrans(z$pars.true[pars.resp][i-1]), backtrans(z$pars.est[pars.resp][i-1])))
        
        }
        
        if(is.numeric(intpar)) {
          return(intpar)
        } else {
          return(rep(NA, 2))
        }
      })
      
      # rows <- length(sqerr)
      # 
      # pars <- matrix(unlist(sqerr), nrow = rows)
      # 
      # return(colMeans(pars, na.rm = T))
    }, y = x)
  })
}
