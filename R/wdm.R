#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# wdm.R
#
# (c) 2016 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2016-11-04
# last mod 2016-12-20 18:39 DW
#

bwdm <- function(data, yvar=c("q", "resp"), alpha=NULL, tau=NULL, beta=NULL, delta=NULL,
               xvar=NULL, start=NULL, fixed=0, jagscmd=NULL) {
  # save original function call
  cl <- match.call()

  # prepare passed arguments
  #RWiener::verifydata(data)
  if (is.numeric(data) & is.null(xvar)) data <- RWiener::revamp(data, yvar=yvar)
  else if (length(yvar)==1) {
    cbind(RWiener::revamp(data[,yvar]), data)
    yvar <- c("q", "resp")
  }
  fpar <- c("alpha"=unname(alpha), "tau"=unname(tau), 
    "beta"=unname(beta), "delta"=unname(delta))

  if(is.null(data$id)) data$id <- factor(1)

  # estimate parameters
  if (!is.null(xvar)) {
    if(length(xvar)==1) {
      if(class(data[,xvar]) == "factor"){
        res <- list()
        res$coefficients <- fpar
        for (l in levels(data[,xvar])) {
          est <- mcmce(data[data[,xvar]==l,c(yvar, "id")], fpar, start, jagscmd)
          est$coefficients <- est$coefficients[!(names(est$coefficients) %in% names(fpar))]
          names(est$coefficients) <- paste(l,names(est$coefficients), sep=":")
          res$coefficients <- append(res$coefficients, est$coefficients)
          res$counts <- append(res$counts, est$counts)
          res$algorithm <- append(res$algorithm, list(est$algorithm))
          res$convergence <- append(res$convergence, est$convergence)
          res$message <- append(res$message, est$message)
          res$hessian <- append(res$hessian, list(est$hessian))
          res$loglik <- sum(res$loglik, est$loglik)
        }
      } else stop("xvar has to be a factor")
    }
  }
  else
    res <- mcmce(data[,c(yvar, "id")], fpar, start, jagscmd)

  # prepare return object
  res$nobs <- length(data[,1])
  res$npar <- length(res$coefficients)-fixed
  res$data <- data
  res$yvar <- yvar
  res$estpar <- c("alpha"=is.null(alpha), "tau"=is.null(tau),
                   "beta"=is.null(beta), "delta"=is.null(delta))
  res$call <- cl
  class(res) <- c("bwdm", "wdm")
  return(res)
}

## internal function
mcmce <- function(data, fpar, start=NULL, jagscmd=NULL) {
  rjags::load.module("wiener")
  if(is.null(jagscmd)) jagscmd <- list(n.chains=1, n.adapt=1000, n.iter=10000, thin=10)

  data$y <- RWiener::revamp(data)

  data$id <- factor(data$id)
  I <- length(levels(data$id))
  J <- max(aggregate(y~id, data, length)$y)
  yvar <- matrix(NA, I, J)
  for(i in 1:I) {
    for(j in 1:length(data$y[data$id==levels(data$id)[i]])) yvar[i,j] <- data$y[data$id==i][j]
  }

  if(is.null(start)) {
    jinits <- list()
    for(ch in 1:jagscmd$n.chains) {
      inits <- list()
      if(!("alpha" %in% names(fpar))) inits <- append(inits, list(alpha=runif(1,0.01,2)))
      if(!("tau" %in% names(fpar))) inits <- append(inits, list(tau=min(data[,1])/3))
      if(!("beta" %in% names(fpar))) inits <- append(inits, list(beta=runif(1,.2,.8)))
      if(!("delta" %in% names(fpar))) inits <- append(inits, list(delta=runif(1,-1,1)))
      jinits[[ch]] <- inits
    }
  }

  jdata <- list(y=yvar, I=I, J=J)

  modelfile <- textConnection("wdmodel", "w")
  cat("model { # wiener diffusion model\n
  ", file=modelfile)

  if (is.null(fpar)) {
    cat("alpha ~ dunif(0.01,100)\n  ", file=modelfile)
    cat("beta ~ dunif(0.01,1)\n  ", file=modelfile)
    cat("tau ~ dunif(0.01,1)\n  ", file=modelfile)
    cat("delta ~ dunif(-100,100)\n  ", file=modelfile)
  }
  else {
    if("alpha" %in% names(fpar))
      cat(paste0("alpha <- ", fpar["alpha"], "\n  "), file=modelfile)
    else {
      cat("alpha ~ dunif(0.01,100)\n  ", file=modelfile)
    }
    if("tau" %in% names(fpar)) 
      cat(paste0("tau <- ", fpar["tau"], "\n  "), file=modelfile)
    else {
      cat("tau ~ dunif(0.01,1)\n  ", file=modelfile)
    }
    if("beta" %in% names(fpar)) 
      cat(paste0("beta <- ", fpar["beta"], "\n  "), file=modelfile)
    else {
      cat("beta ~ dunif(0.01,1)\n  ", file=modelfile)
    }
    if("delta" %in% names(fpar)) 
      cat(paste0("delta <- ", fpar["delta"], "\n  "), file=modelfile)
    else {
      cat("delta ~ dunif(-100,100)\n  ", file=modelfile)
    }
  }

  cat("
	for(i in 1:I) {              # I <- number of subjects
    for(j in 1:J) {            # J <- number of trials
      y[i,j] ~ dwiener(alpha, tau, beta, delta)
    }
	}
  \n}", file=modelfile)
  close(modelfile)

  jmodel <- rjags::jags.model(textConnection(wdmodel), jdata, jinits, jagscmd$n.chains, jagscmd$n.adapt, TRUE)
  jsamples <- rjags::coda.samples(jmodel, c("alpha", "tau", "beta", "delta"), #names(jinits[[1]]), 
                                  jagscmd$n.iter, jagscmd$thin)

  res <- list(
    coefficients = apply(coda::as.array.mcmc.list(jsamples), 2, mean),
    sd = apply(coda::as.array.mcmc.list(jsamples), 2, sd),
    loglik = NA,
    algorithm = list(type="MCMC sampling with JAGS", samples=jsamples)
  )

  return(res)
}
