#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# wdm.R
#
# (c) 2016 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2016-11-04
# last mod 2016-12-22 16:32 DW
#

bwdm <- function(data, yvar=c("q", "resp"), alpha=NULL, tau=NULL, beta=NULL, delta=NULL,
               xvar=NULL, start=NULL, fixed=0, jagscmd=NULL, idvar="id") {
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

  if(!(idvar == "id")) data$id <- factor(data[[idvar]])
  else if(is.null(data$id)) data$id <- factor(1)
  else data$id <- factor(data$id)

  # estimate parameters
  if (!is.null(xvar)) {
    if(length(xvar)==1) {
      if(class(data[,xvar]) == "factor"){
        res <- mcmce(data[,c(yvar, "id", xvar)], fpar, start, jagscmd)
      } else stop("xvar has to be a factor")
    }
  }
  else
    res <- mcmce(data[,c(yvar, "id")], fpar, start, jagscmd)

  # prepare return object
  res$nobs <- length(data[,1])
  res$data <- data
  res$yvar <- yvar
  res$estpar <- c("alpha"=is.null(alpha), "tau"=is.null(tau),
                   "beta"=is.null(beta), "delta"=is.null(delta))
  res$coefficients <- append(fpar, res$coefficients)
  res$npar <- length(res$coefficients)-fixed
  res$call <- cl
  class(res) <- c("bwdm", "wdm")
  return(res)
}

## internal function
mcmce <- function(data, fpar, start=NULL, jagscmd=NULL) {
  rjags::load.module("wiener", quiet=TRUE)
  if(is.null(jagscmd)) jagscmd <- list(n.chains=1, n.adapt=1000, n.iter=10000, thin=10)

  data$y <- RWiener::revamp(data)

  if (length(data) == 4) {
    K <- 1
    I <- length(levels(data$id))
    J <- max(aggregate(y~id, data, length)$y)
    yvar <- matrix(NA, I, J)
    for(i in 1:I) {
      for(j in 1:length(data$y[data$id==levels(data$id)[i]])) yvar[i,j] <- data$y[data$id==levels(data$id)[i]][j]
    }
    yvar <- array(yvar, c(I,J,K))
  }
  else if (length(data) == 5)
  {
    names(data)[4] <- "xvar" # the 4th variable should always be xvar
    K <- length(levels(data$xvar))
    I <- length(levels(data$id)) 
    J <- max(aggregate(y~id+xvar, data, length)$y)
    yvar <- array(NA, c(I,J,K))
    for(k in 1:K) {
      for(i in 1:I) {
        for(j in 1:length(data$y[data$id==levels(data$id)[i] & data$xvar==levels(data$xvar)[k]])) 
          yvar[i,j,k] <- 
          data$y[data$id==levels(data$id)[i] & data$xvar==levels(data$xvar)[k]][j]
      }
    }
    yvar <- array(yvar, c(I,J,K))
  }

  if(is.null(start)) {
    jinits <- list()
    for(ch in 1:jagscmd$n.chains) {
      inits <- list()
      if(!("alpha" %in% names(fpar))) inits <- append(inits, list(alpha=runif(K,0.5,2)))
      if(!("tau" %in% names(fpar))) inits <- append(inits, list(tau=rep(min(data[,1])/3,K)))
      if(!("beta" %in% names(fpar))) inits <- append(inits, list(beta=runif(K,.2,.8)))
      if(!("delta" %in% names(fpar))) inits <- append(inits, list(delta=runif(K,-1,1)))
      jinits[[ch]] <- inits
    }
  }

  jdata <- list(y=yvar, I=I, J=J, K=K)

  modelfile <- textConnection("wdmodel", "w")
  cat("model { # wiener diffusion model\n
  ", file=modelfile)

  if (is.null(fpar)) {
    mon <- c("alpha", "tau", "beta", "delta")
    cat("for(k in 1:K) {\n    ", file=modelfile)
    cat("alpha[k] ~ dunif(0.01,100)\n    ", file=modelfile)
    cat("beta[k] ~ dunif(0.01,1)\n    ", file=modelfile)
    cat("tau[k] ~ dunif(0.01,1)\n    ", file=modelfile)
    cat("delta[k] ~ dunif(-100,100)\n  ", file=modelfile)
    cat("}\n  ", file=modelfile)
  }
  else {
    mon <- vector()
    cat("for(k in 1:K) {\n    ", file=modelfile)
    if("alpha" %in% names(fpar))
      cat(paste0("alpha[k] <- ", fpar["alpha"], "\n  "), file=modelfile)
    else {
      mon <- c(mon, "alpha")
      cat("alpha[k] ~ dunif(0.01,100)\n  ", file=modelfile)
    }
    if("tau" %in% names(fpar)) 
      cat(paste0("tau[k] <- ", fpar["tau"], "\n  "), file=modelfile)
    else {
      mon <- c(mon, "tau")
      cat("tau[k] ~ dunif(0.01,1)\n  ", file=modelfile)
    }
    if("beta" %in% names(fpar)) 
      cat(paste0("beta[k] <- ", fpar["beta"], "\n  "), file=modelfile)
    else {
      mon <- c(mon, "beta")
      cat("beta[k] ~ dunif(0.01,1)\n  ", file=modelfile)
    }
    if("delta" %in% names(fpar)) 
      cat(paste0("delta[k] <- ", fpar["delta"], "\n  "), file=modelfile)
    else {
      mon <- c(mon, "delta")
      cat("delta[k] ~ dunif(-100,100)\n  ", file=modelfile)
    }
    cat("}\n  ", file=modelfile)
  }

  cat("
	for(i in 1:I) {              # I <- number of subjects
    for(j in 1:J) {            # J <- number of trials
      for(k in 1:K) {          # K <- number of conditions (xvar)
        y[i,j,k] ~ dwiener(alpha[k], tau[k], beta[k], delta[k])
      }
    }
	}
  \n}", file=modelfile)
  close(modelfile)

  jmodel <- rjags::jags.model(textConnection(wdmodel), jdata, jinits, jagscmd$n.chains, jagscmd$n.adapt, TRUE)
  jsamples <- rjags::coda.samples(jmodel, mon, #names(jinits[[1]]), 
                                  jagscmd$n.iter, jagscmd$thin)

  res <- list(
    coefficients = apply(coda::as.array.mcmc.list(jsamples), 2, mean),
    sd = apply(coda::as.array.mcmc.list(jsamples), 2, sd),
    loglik = NA,
    algorithm = list(type="MCMC sampling with JAGS", samples=jsamples)
  )

  rjags::unload.module("wiener", quiet=TRUE)
  return(res)
}
