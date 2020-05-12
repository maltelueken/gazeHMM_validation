### Model functions ###

# Author: Malte Lüken
# Date: 01.04.2020

library(depmixS4)
library(CircStats)


# Create response model for uniform distribution --------------------------

setClass("unif", contains = "response")

setGeneric("unif", function(y, ...) standardGeneric("unif"))


# Define the method that creates the response class

setMethod("unif", 
          signature(y = "ANY"), 
          function(y, fixed = NULL, ...) {
            
            y <- matrix(y, length(y))
            x <- matrix(1)
            parameters <- list()
            npar <- 2 # min and max
            
            if(is.null(fixed)) fixed <- as.logical(rep(0,npar))
            
            parameters$min <- min(y)
            parameters$max <- max(y)
            
            mod <- new("unif", parameters = parameters, fixed = fixed, x = x, y = y, npar = npar)
            
            return(mod)
          }
)


# Define method that prints parameters

setMethod("show", "unif",
          function(object) {
            
            cat("Model of type uniform\n")
            cat("Parameters: \n")
            cat("min: ", object@parameters$min, "\n")
            cat("max: ", object@parameters$max, "\n")
          }
)


# Define method that evaluates the density

setMethod("dens", "unif",
          function(object, log = FALSE) {
            
            dens <- ifelse(is.na(object@y), 1, dunif(object@y, min = object@parameters$min, max = object@parameters$max, log = log))
            
          }
)


# Define method that retrieves parameters

setMethod("getpars", "response",
          function(object, which = "pars", ...) {
            
            switch(which,
                   "pars" = {
                     parameters <- numeric()
                     parameters <- unlist(object@parameters)
                     pars <- parameters
                   },
                   "fixed" = {
                     pars <- object@fixed
                   }
            )
            return(pars)
          }
)


# Define method that sets parameters

setMethod("setpars", "unif",
          function(object, values, which = "pars", ...) {
            
            npar <- npar(object)
            
            if(length(values) != npar) stop("length of 'values' must be", npar)
            
            nms <- names(object@parameters)
            
            switch(which,
                   "pars"= {
                     object@parameters$min <- values[1]
                     object@parameters$max <- values[2]
                   },
                   "fixed" = {
                     object@fixed <- as.logical(values)
                   }
            )
            names(object@parameters) <- nms
            
            return(object)
          }
)


# Define method that estimates the parameters

setMethod("fit", "unif",
          function(object, w) {
            
            if(missing(w)) w <- NULL
            
            y <- object@y
            pars <- c(min(y[!is.na(y)]), max(y[!is.na(y)]))
            object <- setpars(object, pars)
            
            return(object)
          }
)


# Define method that predicts observations

setMethod("predict", "unif", 
          function(object) {
            
            ret <- sample(object@y, 1)
            
            return(ret)
          }
)


# Define method to simulate

setMethod("simulate", "unif",
          function(object, nsim = 1, seed) {
            
            if(!is.null(seed)) set.seed(seed)
            
            nt <- nrow(object@y)
            
            sim <- runif(n = nt*nsim, min = object@parameters$min, max = object@parameters$max)
            
            return(as.matrix(sim))
          }
)


# Create response model for von-Mises distribution ------------------------

setClass("vMF", contains = "response")

setGeneric("vMF", function(y, pstart = NULL, fixed = NULL, ...) standardGeneric("vMF"))


# Define the method that creates the response class

setMethod("vMF", 
          signature(y = "ANY"), 
          function(y, pstart=NULL, fixed=NULL, ...) {
            
            y <- matrix(y, length(y))
            x <- matrix(1)
            parameters <- list()
            npar <- 2 # mu and kappa
            
            if(is.null(fixed)) fixed <- as.logical(rep(0, npar))
            
            if(!is.null(pstart)) {
              
              if(length(pstart) != npar) stop("length of 'pstart' must be ", npar)
              
              parameters$mu <- pstart[1]
              parameters$kappa <- log(pstart[2])
            }
            
            mod <- new("vMF", parameters = parameters, fixed = fixed, x = x, y = y, npar = npar)
            
            return(mod)
          }
)


# Define method that prints the parameters

setMethod("show", "vMF",
          function(object) {
            cat("Model of type vMF\n")
            cat("Parameters: \n")
            cat("mu: ", object@parameters$mu, "\n")
            cat("kappa: ", exp(object@parameters$kappa), "\n")
          }
)


# Define method that evaluates the density

setMethod("dens","vMF",
          function(object) {
            
            dens <- ifelse(is.na(object@y), 1, dvm(object@y, mu = object@parameters$mu, kappa = exp(object@parameters$kappa)))
            
          }
)


# Define method that retrieves the parameters

setMethod("getpars","response",
          function(object,which="pars",...) {
            
            switch(which,
                   "pars" = {
                     parameters <- numeric()
                     parameters <- unlist(object@parameters)
                     pars <- parameters
                   },
                   "fixed" = {
                     pars <- object@fixed
                   }
            )
            
            return(pars)
          }
)


# Define method that sets the parameters

setMethod("setpars","vMF",
          function(object, values, which="pars", ...) {
            
            npar <- npar(object)
            
            if(length(values)!=npar) stop("length of 'values' must be",npar)
            
            nms <- names(object@parameters)
            
            switch(which,
                   "pars"= {
                     object@parameters$mu <- values[1]
                     object@parameters$kappa <- values[2]
                   },
                   "fixed" = {
                     object@fixed <- as.logical(values)
                   }
            )
            
            names(object@parameters) <- nms
            
            return(object)
          }
)


# Define method that estimates the parameters

setMethod("fit","vMF",
          function(object,w) {
            
            if(missing(w)) w <- NULL
            
            y <- object@y
            
            
            # Create weighted log-likelihood function for von-Mises distribution
            
            ll_vMF <- function(par, y, w) {
              
              if(is.null(w)) w <- 1 # no weights: all values weighted equally (= 1)
              
              miss <- is.na(y)
              
              dens <- ifelse(miss, 1, dvm(y, mu = par[1], exp(kappa = par[2])))
              
              return(-sum(w*log(dens)))
            }
            
            init <- c(object@parameters$mu, object@parameters$kappa) # start values
            
            
            # Optimize weighted ll for von-Mises distribution
            
            fit <- BB::BBoptim(par = init, fn = ll_vMF, y = y, w = w,
                               lower = c(-Inf, -Inf), upper = c(Inf, Inf), quiet = T)
            
            pars <- fit$par
            object <- setpars(object,pars)
            
            return(object)
          }
)


# Define method that predicts observations

setMethod("predict","vMF", 
          function(object) {
            
            ret <- object@parameters$mu
            
            return(ret)
          }
)


# Define method to simulate

setMethod("simulate", signature(object = "vMF"),
          function(object, nsim = 1, seed) {
            
            if(!is.null(seed)) set.seed(seed)
            
            nt <- nrow(object@y)
            
            sim <- rvm(nt*nsim, mean = object@parameters$mu, k = exp(object@parameters$kappa))
            
            return(as.matrix(sim))
          }
)


# Create alternative response model for gamma distribution ----------------

setClass("altGamma", contains = "response")

setGeneric("altGamma", function(y, pstart = NULL, fixed = NULL, ...) standardGeneric("altGamma"))


# Define the method that creates the response class

setMethod("altGamma", 
          signature(y="ANY"), 
          function(y,pstart=NULL,fixed=NULL, ...) {
            y <- matrix(y,length(y))
            x <- matrix(1)
            parameters <- list()
            npar <- 2
            if(is.null(fixed)) fixed <- as.logical(rep(0,npar))
            if(!is.null(pstart)) {
              if(length(pstart)!=npar) stop("length of 'pstart' must be ",npar)
              parameters$shape <- log(pstart[1])
              parameters$scale <- log(pstart[2])
            }
            mod <- new("altGamma",parameters=parameters,fixed=fixed,x=x,y=y,npar=npar)
            mod
          }
)

setMethod("show","altGamma",
          function(object) {
            cat("Model of type altGamma (see ?altGammamlss for details) \n")
            cat("Parameters: \n")
            cat("shape: ", exp(object@parameters$shape), "\n")
            cat("scale: ", exp(object@parameters$scale), "\n")
          }
)

setMethod("dens","altGamma",
          function(object,log=FALSE) {
            
            dens <- ifelse(is.na(object@y), 1, dgamma(object@y, shape = exp(object@parameters$shape), 
                                                      scale = exp(object@parameters$scale),log = log))
            
          }
)

setMethod("getpars","response",
          function(object,which="pars",...) {
            switch(which,
                   "pars" = {
                     parameters <- numeric()
                     parameters <- unlist(object@parameters)
                     pars <- parameters
                   },
                   "fixed" = {
                     pars <- object@fixed
                   }
            )
            return(pars)
          }
)

setMethod("setpars","altGamma",
          function(object, values, which="pars", ...) {
            npar <- npar(object)
            if(length(values)!=npar) stop("length of 'values' must be",npar)
            # determine whether parameters or fixed constraints are being set
            nms <- names(object@parameters)
            switch(which,
                   "pars"= {
                     object@parameters$shape <- values[1]
                     object@parameters$scale <- values[2]
                   },
                   "fixed" = {
                     object@fixed <- as.logical(values)
                   }
            )
            names(object@parameters) <- nms
            return(object)
          }
)

setMethod("fit","altGamma",
          function(object,w) {
            
            if(missing(w)) w <- NULL
            
            y <- object@y
          
            
            # Create weighted log-likelihood function for gamma distribution
            
            ll_gamma <- function(par, y, w) {
              
              if(is.null(w)) w <- 1 # no weights: all values weighted equally (= 1)
              
              # w[w == 0] <- 1e-5
              
              miss <- is.na(y)
              
              dens <- ifelse(miss, log(1), dgamma(y, shape = exp(par[1]), scale = exp(par[2]), log = T))
              
              return(-sum(w*dens))
            }
            
            
            init <- c(object@parameters$shape, object@parameters$scale) # start values
            
            
            # Optimize weighted ll for von-Mises distribution
            
            fit <- BB::BBoptim(par = init, fn = ll_gamma, y = y, w = w,
                         lower = c(-Inf, -Inf), upper = c(Inf, Inf), quiet = T)
            
            pars <- fit$par
            object <- setpars(object, pars)
            
            return(object)
          }
)

setMethod("predict","altGamma", 
          function(object) {
            ret <- object@parameters$shape*object@parameters$scale
            return(ret)
          }
)


# Define method to simulate

setMethod("simulate", signature(object = "altGamma"),
          function(object, nsim = 1, seed) {
            
            if(!is.null(seed)) set.seed(seed)
            
            nt <- nrow(object@y)
            
            sim <- rgamma(nt*nsim, shape = exp(object@parameters$shape), scale = exp(object@parameters$scale))
            
            return(as.matrix(sim))
          }
)