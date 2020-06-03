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


# Calculate normalized root mean square deviation

rmsd <- list()

for (part in 1:4) {
  
  rmsd[[part]] <- lapply(get(paste("estimates.", part, sep = "")), function(x) {
    lapply(x, function(y) {
      sqerr <- lapply(y, function(z) {
        
        err <- try(((backtrans(z$pars.est) - backtrans(z$pars.true))/
                      ifelse(backtrans(z$pars.true) == 0, 2*pi, backtrans(z$pars.true)))^2)
        
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
      
      msqerr <- apply(pars, 2, median, na.rm = T)
      
      names(msqerr) <- nms
      
      rmsqerr <- sqrt(msqerr)
      
      # parstrue <- lapply(y, function(z) {backtrans(z$pars.true)})
      # 
      # parstruemat <- matrix(unlist(parstrue), nrow = rows, byrow = T)
      # 
      # meantrue <- apply(parstruemat, 2, mean, na.rm = T)
      # 
      # meantrue <- ifelse(meantrue == 0, 2*pi, meantrue)
      # 
      # print(rmsqerr/meantrue)
      
      return(rmsqerr)
    })
  })
}


# Display MSE in data frame

rmsd.data <- lapply(rmsd, function(x) lapply(x, as.data.frame))
rmsd.data <- lapply(rmsd.data, function(x) lapply(x, function(y) {as.data.frame(t(as.matrix(y)))}))
rmsd.data <- lapply(rmsd.data, function(x) lapply(1:length(x), function(y, data) {
  names(data[[y]]) <- names(rmsd[[1]][[y]][[1]])
  return(data[[y]])}, data = x))


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
})


# Summarise accuracy

acc.data <- list()

for (part in 1:4) {
  
  acc.data[[part]] <- lapply(get(paste("estimates.", part, sep = "")), function(x) {
    out <- lapply(x, function(y) { 
      out <- lapply(y, function(z) {
        
        if(is.numeric(z$accuracy)) {
          acc <- z$accuracy
        } else {
          acc <- NA
        }
        
        return(acc)
      })
      
      return(as.vector(reduce(out, cbind)))
    })
    
    return(as.data.frame(t(reduce(out, cbind))))
  })
}


# Plot RMdSPD part 1

plots.rmsd.1 <- lapply(rmsd.data[[1]], function(x) {
  
  names.pars.varied <- list(bquote(a["i=j"]), bquote(alpha["vel"]), bquote(beta["vel"]), bquote(alpha["acc"]),
                            bquote(beta["acc"]), bquote(kappa))
  
  if(ncol(x) == 18) {
    
    trnames <- list(bquote(rho["i"]), bquote(a["i1"]), bquote(a["i2"]))
    
  } else if(ncol(x) == 30) {
    
    trnames <- list(bquote(rho["i"]), bquote(a["i1"]), bquote(a["i2"]), bquote(a["i3"]))
    
  } else {
    
    trnames <- list(bquote(rho["i"]), bquote(a["i1"]), bquote(a["i2"]), bquote(a["i3"]), bquote(a["i4"]))
    
  }
  
  names.pars.est <- append(trnames, list(bquote(alpha["vel"]), bquote(beta["vel"]), bquote(alpha["acc"]),
                         bquote(beta["acc"]), bquote(mu), bquote(kappa)))
  
  x <- as_tibble(x, .name_repair = "unique")
  
  data.long <- x %>%
    mutate(par.varied = c(0, 1, 2, 3, 4, rep(c(1, 2, 3, 4, 5), (nrow(x) %/% 5)-1)),
           state.varied = c(1, rep(1, 4), rep(2:((nrow(x) %/% 5)), each = 5))) %>%
    pivot_longer(names(x), names_to = "par.est", values_to = "RMdSPD", names_repair = "unique") %>%
    mutate(RMdSPD = cut(RMdSPD, breaks = c(0, 0.1, 0.5, Inf), c("good", "moderate", "bad")),
           state.est = rep(c(rep(1:max(state.varied), max(state.varied)+1),
                         rep(1:max(state.varied), each = 6)), nrow(x)),
           par.est = rep(c(rep(1, max(state.varied)), 
                           rep(2:(max(state.varied)+1), each = max(state.varied)), 
                           rep((max(state.varied)+2):(max(state.varied)+7), max(state.varied))), nrow(x))) %>%
    mutate_at(vars(par.varied, state.varied, par.est, state.est), as.factor)

  p <- ggplot(data = data.long, aes(x = par.varied, y = par.est, fill = RMdSPD)) +
    geom_tile() + facet_grid(cols = vars(state.varied), rows = vars(state.est)) +
    scale_x_discrete(labels = names.pars.varied) +
    scale_y_discrete(labels = names.pars.est)
    theme(axis.text = element_text(size = 14))
  
  return(p)
})

print(plots.rmsd.1[[3]])


# Plot RMdSPD parts 2, 3, and 4

plots.rmsd.234 <- lapply(2:4, function(y, data) lapply(data[[y]], function(x) {
  
  if(ncol(x) == 18) {
    
    k <- 2
    
    trnames <- list(bquote(rho["i"]), bquote(a["i1"]), bquote(a["i2"]))
    
  } else if(ncol(x) == 30) {
    
    k <- 3
    
    trnames <- list(bquote(rho["i"]), bquote(a["i1"]), bquote(a["i2"]), bquote(a["i3"]))
    
  } else {
    
    k <-4
    
    trnames <- list(bquote(rho["i"]), bquote(a["i1"]), bquote(a["i2"]), bquote(a["i3"]), bquote(a["i4"]))
    
  }
  
  names.pars.est <- append(trnames, list(bquote(alpha["vel"]), bquote(beta["vel"]), bquote(alpha["acc"]),
                                         bquote(beta["acc"]), bquote(mu), bquote(kappa)))
  
  x <- as_tibble(x, .name_repair = "unique")
  
  data.long <- x %>%
    mutate(cond = 1:3) %>%
    pivot_longer(names(x), names_to = "par.est", values_to = "RMdSPD", names_repair = "unique") %>%
    mutate(state.est = rep(c(rep(1:k, k+1),
                             rep(1:k, each = 6)), nrow(x)),
           par.est = rep(c(rep(1, k), 
                           rep(2:(k+1), each = k), 
                           rep((k+2):(k+7), k)), nrow(x))) %>%
    mutate_at(vars(cond, par.est, state.est), as.factor)
  
  p <- ggplot(data = data.long, aes(x = par.est, y = RMdSPD, color = cond)) +
    geom_point(position = position_dodge(0.25)) + facet_grid(cols = vars(state.est)) +
    scale_x_discrete(labels = names.pars.est) +
    scale_y_continuous(breaks = c(0.1, 0.5, 1, 5, 10)) +
    geom_hline(yintercept = 0.1, linetype = "dashed") +
    geom_hline(yintercept = 0.5, linetype = "dashed")
  
  if(y == 2) {
    names.cond <- c("500", "2500", "10000")
    label.cond <- "N"
  } else if (y == 3) {
    names.cond <- c("1", "2", "3")
    label.cond <- bquote(tau["start"])
  } else {
    names.cond <- c("1", "3", "5")
    label.cond <- "m"
  }
  
  p <- p + scale_color_discrete(name = label.cond, labels = names.cond)
  
  return(p)
}), data = rmsd.data)

print(plots.rmsd.234[[1]][[1]])


# Plot linear regressions for transition probabilities

D <- 100

regw.tr.data <- lapply(regw.tr, function(x) lapply(x, function(y) rbind(y)))
regw.tr.data <- lapply(regw.tr.data, function(x) reduce(x, rbind))

plots.lm.tr <- lapply(1:length(regw.tr.data), function(x, y) {
  
  data <- y[[x]] %>% 
    mutate(par = rep(1:(x+1)^2, each = D),
           from = rep(rep(1:(x+1), each = D), x+1),
           to = rep(1:(x+1), each = D*(x+1)))
  
  p <- ggplot(data, aes(x = true, y = est)) + 
    facet_grid(rows = vars(from), cols = vars(to), labeller = label_both) +
    geom_point() + geom_smooth(method = "lm") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed")
  
  return(p)
}, y = regw.tr.data)

print(plots.lm.tr[[3]])


# Plot linear regressions for response parameters

regw.resp.data <- lapply(regw.resp, function(x) lapply(x, function(y) rbind(y)))
regw.resp.data <- lapply(regw.resp.data, function(x) reduce(x, rbind))

parnames <- list(bquote(alpha["vel;1"]), bquote(beta["vel;1"]), bquote(alpha["acc;1"]), bquote(beta["acc;1"]),
                 bquote(alpha["vel;2"]), bquote(beta["vel;2"]), bquote(alpha["acc;2"]), bquote(beta["acc;2"]), bquote(kappa["2"]),
                 bquote(alpha["vel;3"]), bquote(beta["vel;3"]), bquote(alpha["acc;3"]), bquote(beta["acc;3"]), bquote(kappa["3"]),
                 bquote(alpha["vel;4"]), bquote(beta["vel;4"]), bquote(alpha["acc;4"]), bquote(beta["acc;4"]), bquote(kappa["4"]))

plots.lm.resp <- lapply(regw.resp.data, function(x) {
  
  npar <- nrow(x) %/% D
  
  data <- x %>% 
    mutate(par = rep(1:npar, each = D),
           type = rep(1, 2, 1, 2))
  
  p <- ggplot(data, aes(x = true, y = est)) + 
    facet_wrap(vars(par), scales = "free", labeller = label_bquote(.(parnames[[par]]))) +
    geom_point() + geom_smooth(method = "lm") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed")
  
  return(p)
})

print(plots.lm.resp[[3]])


# Plot accuracy part 1

plots.acc.1 <- lapply(acc.data[[1]], function(x) {
  
  names.pars.varied <- list(bquote(a["i=j"]), bquote(alpha["vel"]), bquote(beta["vel"]), bquote(alpha["acc"]),
                            bquote(beta["acc"]), bquote(kappa))
  
  x <- as_tibble(x, .name_repair = "unique")
  
  data.long <- x %>%
    mutate(par.varied = c(0, 1, 2, 3, 4, rep(c(1, 2, 3, 4, 5), (nrow(x) %/% 5)-1)),
           state.varied = c(1, rep(1, 4), rep(2:((nrow(x) %/% 5)), each = 5))) %>%
    pivot_longer(names(x), names_to = "par.est", values_to = "accuracy") %>%
    mutate_at(vars(par.varied, state.varied), as.factor)
  
  p <- ggplot(data = data.long, aes(x = par.varied, y = accuracy)) +
    geom_boxplot() + facet_grid(cols = vars(state.varied)) +
    scale_x_discrete(labels = names.pars.varied)
  
  return(p)
})

print(plots.acc.1[[2]])


# Plot accuracy parts 2, 3, and 4

acc.data.234 <- lapply(acc.data[2:4], function(x) reduce(x, rbind))

plots.acc.234 <- lapply(1:3, function(y, data) {
  
  x <- as_tibble(data[[y]], .name_repair = "unique")
  
  data.long <- x %>%
    mutate(cond = as.factor(rep(1:3, 3)),
           k = as.factor(rep(2:4, each = 3))) %>%
    pivot_longer(names(x), names_to = "par.est", values_to = "accuracy")
  
  p <- ggplot(data = data.long, aes(x = k, y = accuracy, color = cond)) +
    geom_boxplot() + 
    scale_x_discrete(name = "k (number of states)")
  
  if(y == 1) {
    names.cond <- c("500", "2500", "10000")
    label.cond <- "N"
  } else if (y == 2) {
    names.cond <- c("1", "2", "3")
    label.cond <- bquote(tau["start"])
  } else {
    names.cond <- c("1", "3", "5")
    label.cond <- "m"
  }
  
  p <- p + scale_color_discrete(name = label.cond, labels = names.cond)
  
  return(p)
}, data = acc.data.234)

print(plots.acc.234[[2]])


# Plot accuracy over interval parts 2 and 4

plots.acc.int.24 <- lapply(c(1, 3), function(y, data) {
  
  x <- as_tibble(data[[y]], .name_repair = "unique")
  
  if(y == 1) {
  
    int <- rep(seq(1, 5, length.out = D), 9) 
  
  } else {
    
    int <- rep(floor(seq(1, 200, length.out = D)), 9)
    
  }
  
  data.long <- x %>%
    mutate(cond = as.factor(rep(1:3, 3)),
           k = as.factor(rep(2:4, each = 3))) %>%
    pivot_longer(names(x), names_to = "par.est", values_to = "accuracy") %>%
    mutate(int = int)
  
  p <- ggplot(data = data.long, aes(x = int, y = accuracy, color = cond)) +
    facet_grid(cols = vars(k)) + 
    geom_point() # + geom_smooth(method = "lm")
  
  if(y == 1) {
    names.cond <- c("500", "2500", "10000")
    label.cond <- "N"
    x.name <- bquote(tau["noise"])
  } else {
    names.cond <- c("1", "3", "5")
    label.cond <- "m"
    x.name <- "l"
  }
  
  p <- p + scale_color_discrete(name = label.cond, labels = names.cond) +
    scale_x_continuous(name = x.name)
  
  return(p)
}, data = acc.data.234)

print(plots.acc.int.24[[2]])


# Plot proportion of erroneous models part 1

plots.err.1 <- lapply(acc.data[[1]], function(x) {
  
  names.pars.varied <- list(bquote(a["i=j"]), bquote(alpha["vel"]), bquote(beta["vel"]), bquote(alpha["acc"]),
                            bquote(beta["acc"]), bquote(kappa))
  
  x <- as_tibble(x, .name_repair = "unique")
  
  data.long <- x %>%
    mutate(par.varied = c(0, 1, 2, 3, 4, rep(c(1, 2, 3, 4, 5), (nrow(x) %/% 5)-1)),
           state.varied = c(1, rep(1, 4), rep(2:((nrow(x) %/% 5)), each = 5))) %>%
    pivot_longer(names(x), names_to = "par.est", values_to = "accuracy") %>%
    mutate_at(vars(par.varied, state.varied), as.factor) %>%
    group_by(par.varied, state.varied) %>%
    summarise(error = mean(is.na(accuracy)))
  
  p <- ggplot(data = data.long, aes(x = par.varied, y = error)) +
    geom_point() + facet_grid(cols = vars(state.varied)) +
    scale_x_discrete(labels = names.pars.varied) + 
    scale_y_continuous(name = "proportion erroneous models", limits = c(0, 1))
  
  return(p)
})

print(plots.err.1[[1]])


# Plot proportion of erroneous models parts 2, 3, 4

plots.err.234 <- lapply(1:3, function(y, data) {
  
  x <- as_tibble(data[[y]], .name_repair = "unique")
  
  data.long <- x %>%
    mutate(cond = as.factor(rep(1:3, 3)),
           k = as.factor(rep(2:4, each = 3))) %>%
    pivot_longer(names(x), names_to = "par.est", values_to = "accuracy") %>%
    group_by(cond, k) %>%
    summarise(error = mean(is.na(accuracy)))
  
  p <- ggplot(data = data.long, aes(x = k, y = error, color = cond)) +
    geom_point(position = position_dodge(width = 0.25)) + 
    scale_x_discrete(name = "k (number of states)") +
    scale_y_continuous(name = "proportion erroneous models", limits = c(0, 1))
  
  if(y == 1) {
    names.cond <- c("500", "2500", "10000")
    label.cond <- "N"
  } else if (y == 2) {
    names.cond <- c("1", "2", "3")
    label.cond <- bquote(tau["start"])
  } else {
    names.cond <- c("1", "3", "5")
    label.cond <- "m"
  }
  
  p <- p + scale_color_discrete(name = label.cond, labels = names.cond)
  
  return(p)
}, data = acc.data.234)

print(plots.err.234[[1]])


# Exploratory analysis label switching ------------------------------------

library(psych)

load("simulation/part3_expl.Rdata")

labsw <- lapply(1:length(estimates.3), function(x) {
  lapply(estimates.3[[x]], function(y) {
    lapply(y, function(z) {
      
      kappa <- try(cohen.kappa(z$states)$kappa)
      
      states <- z$states
      
      if(x == 2 & is.numeric(kappa) & kappa < 0.95) {
          
        states$y <- ifelse(z$states$y == 3, 2, ifelse(z$states$y == 2, 3, z$states$y))
          
        kappa <- cohen.kappa(states)$kappa
        
        if(is.numeric(kappa) & kappa < 0.95) {
          
          states$y <- ifelse(z$states$y == 3, 1, ifelse(z$states$y == 1, 3, z$states$y))
          
          kappa <- cohen.kappa(states)$kappa
        }
        
        if(is.numeric(kappa) & kappa < 0.95) {
          
          states$y <- ifelse(z$states$y == 3, 2, ifelse(z$states$y == 1, 3, 1))
          
          kappa <- cohen.kappa(states)$kappa
        }
      }
      
      if(x == 3 & is.numeric(kappa) & kappa < 0.5) {
        
        states$y <- ifelse(z$states$y == 3, 2, ifelse(z$states$y == 2, 3, z$states$y))
        
        kappa <- cohen.kappa(states)$kappa
        
        if(is.numeric(kappa) & kappa < 0.5) {

          states$y <- ifelse(z$states$y == 3, 4, ifelse(z$states$y == 4, 3, z$states$y))

          kappa <- cohen.kappa(states)$kappa
        }

        if(is.numeric(kappa) & kappa < 0.5) {

          states$y <- ifelse(z$states$y == 3, 1, ifelse(z$states$y == 1, 3, z$states$y))

          kappa <- cohen.kappa(states)$kappa
        }
        
        if(is.numeric(kappa) & kappa < 0.5) {
          
          kappa <- cohen.kappa(z$states)$kappa
        }
      }
      
      return(kappa)
    })
  })
})

