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
linkfun <- function(p, base) {
  lfun <- function(p, base) {
    p <- p/sum(p)
    beta <- numeric(length(p))
    if (any(p == 1)) 
      beta[which(p == 1)] = Inf
    else beta[-base] <- log(p[-base]/p[base])
    return(beta)
  }
  if (is.matrix(p)) {
    beta <- t(apply(p, 1, lfun, base = base))
  }
  else {
    beta <- lfun(p, base)
  }
  return(beta)
}

# Function to transform parameters to normal scale

backtrans <- function(x) {
  
  out <- x
  
  nms <- names(x)
  
  out[str_detect(nms, "(Intercept)")] <- as.vector(apply(matrix(x[str_detect(nms, "(Intercept)")],
                                                        ncol = sqrt(length(x[str_detect(nms, "(Intercept)")])),
                                                        byrow = T), 1, linkinv, base = 1))
  
  out[nms %in% c("shape", "scale", "kappa")] <- exp(x[nms %in% c("shape", "scale", "kappa")])
  
  return(out)
}
mlogit

#backtrans(estimates.1$`3`[[1]][[3]]$pars.est)


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


# Calculate regression weights for transition probabilities

regw.tr <- list()

regw.tr <- lapply(get("estimates.1"), function(x) {
  
  varpar <- lapply(x[1], function(y) {
    
    lapply(y, function(z) {
      
      nms <- names(z$pars.true)
      
      pars.tr <- logical(length(z$pars.true))
      
      pars.tr[str_detect(nms, "(Intercept)")] <- T
      
      intpar <- try(cbind(backtrans(z$pars.true[pars.tr]), backtrans(z$pars.est[pars.tr])))
      
      if(is.numeric(intpar)) {
        out <- intpar
      } else {
        out <- matrix(NA, nrow = length(z$pars.true[pars.tr]), ncol = 2)
      }
      
      out <- apply(out, 1, list)
      
      return(out)
    })
  })
  
  df <- list()
  
  for (i in 1:length(varpar[[1]][[1]])) {
    
    df[[i]] <- lapply(varpar, function(y, index) {
      
      lapply(y, function(z) {z[[index]]})
      
    }, index = i)
  }
  
  df <- lapply(df, as.data.frame)

  df <- lapply(df, function(z) {
    
    out <- as.data.frame(t(as.matrix(z)))
    
    names(out) <- c("true", "est")
    
    return(out)
  })

  regweights <- lapply(df, function(z) {

    lmfit <- lm(est ~ true, data = z)

    return(lmfit)
  })
})


# Calculate linear regression weights for response parameters

regw.resp <- list()
  
regw.resp <- lapply(get("estimates.1"), function(x) {
  
  index <- 1:length(x[-1])
  
  varpar <- lapply(index, function(i, y) {
    
    lapply(y[[i]], function(z) {
      
      nms <- names(z$pars.true)
      
      #pars.tr <- logical(length(z$pars.true))
      pars.resp <- logical(length(z$pars.true))
      
      #pars.tr[str_detect(nms, "(Intercept)")] <- T
      pars.resp[nms %in% c("shape", "scale", "kappa")] <- T
      
      intpar <- try(c(backtrans(z$pars.true[pars.resp][i]), backtrans(z$pars.est[pars.resp][i])))
      
      if(is.numeric(intpar)) {
        return(intpar)
      } else {
        return(rep(NA, 2))
      }
    })
  }, y = x[-1])
  
  df <- lapply(varpar, as.data.frame)
  df <- lapply(df, function(z) {

    out <- as.data.frame(t(as.matrix(z)))

    names(out) <- c("true", "est")

    return(out)
  })

  regweights <- lapply(df, function(z) {

    lmfit <- lm(est ~ true, data = z)

    return(coef(lmfit))
  })
})
