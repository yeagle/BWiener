#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# hwdm.R
#
# (c) 2016 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2016-10-25
# last mod 2016-12-21 18:18 DW
#

hwdm <- function(data, jagscmd=NULL, model="hwdm", idvar="id", condvar=NULL) {
  rjags::load.module("wiener")

  if(!(idvar == "id")) data$id <- data[[idvar]]
  if (is.null(data$id)) data$id <- factor(1) # assume one is if there is no id column
  if (is.null(data$y)) data$y <- revamp(data) # add univariate variable to data if missing

  if (model=="hwdm") res <- mlehwdm(data, jagscmd)
  else if (model=="rdm") res <- mlerc(data, jagscmd)
  else if (model=="wdmhd") res <- mlehd(data, jagscmd)
  else if (model=="wdmhdc") {
    if(!is.null(condvar)) {
      if(!(condvar == "cond")) data$cond <- data[["condvar"]]
      res <- mlehdc(data, jagscmd)
    }
    else stop("'condvar' missing for wdmhdc model")
  }
  else stop("model type not supported")
  res$nobs <- length(data[,1])
  res$npar <- length(res$coefficients)
  res$data <- data
  res$call <- match.call()
  class(res) <- c(class(res), "wdm")
  return(res)
}

mlehwdm <- function(data, jagscmd=NULL) {
  if(is.null(jagscmd)) jagscmd <- list(n.chains=1, n.adapt=1000, n.iter=10000, thin=10)

  data$id <- factor(data$id)
  I <- length(levels(data$id))
  J <- max(aggregate(y~id, data, length)$y)
  yvar <- matrix(NA, I, J)
  for(i in 1:I) {
    for(j in 1:length(data$y[data$id==levels(data$id)[i]])) yvar[i,j] <- data$y[data$id==i][j]
  }

  jinits <- list()
  for(ch in 1:jagscmd$n.chains) {
    mAlpha <- runif(1,0.01,2); mTau <- min(data[,1])/3; mBeta <- runif(1,.2,.8)
    uAlpha <- runif(1,0.01,1); uTau <- runif(1,0,min(data[,1])/3); uBeta <- runif(1,0.01,0.1)
    mDelta <- runif(1,-1,1); sDelta <- runif(1,0,1)
    jinits <- append(jinits, list(mAlpha=mAlpha,mTau=mTau,mBeta=mBeta,
                                  uAlpha=uAlpha,uTau=uTau,uBeta=uBeta,
                                  mDelta=mDelta,sDelta=sDelta))
  }

  jdata <- list(y=yvar, I=I, J=J)

  hwdm.jmodel <- textConnection("model { # hierarchical wiener diffusion model

  mAlpha ~ dunif(0.01,100)      # k
  uAlpha ~ dunif(0.0001,100)    # k
  mBeta ~ dunif(0.01,1)         # k
  uBeta ~ dunif(0.0001,1)       # k
  mTau ~ dunif(0.01,1)          # k
  uTau ~ dunif(0.0001,1)        # k

  mDelta ~ dunif(-100,100)      # k
  sDelta ~ dunif(0.0001,100)    # k
  prDelta <- pow(sDelta, -2)

	for(i in 1:I) {               # I <- number of subjects

    alpha[i] ~ dunif(mAlpha-uAlpha/2, mAlpha+uAlpha/2)
    tau[i] ~ dunif(mTau-uTau/2, mTau+uTau/2)
    beta[i] ~ dunif(mBeta-uBeta/2, mBeta+uBeta/2)
    delta[i] ~ dnorm(mDelta, prDelta)

    for(j in 1:J) {             # J <- number of trials
      y[i,j] ~ dwiener(alpha[i], tau[i], beta[i], delta[i])
    }
	}

}")

  jmodel <- rjags::jags.model(hwdm.jmodel, jdata, jinits, jagscmd$n.chains, jagscmd$n.adapt, TRUE)
  jsamples <- rjags::coda.samples(jmodel, c("mAlpha", "mTau", "mBeta", "uAlpha", "uTau", "uBeta", "mDelta", "sDelta"), 
                                  jagscmd$n.iter, jagscmd$thin)

  res <- list(
    coefficients = apply(coda::as.array.mcmc.list(jsamples), 2, mean),
    sd = apply(coda::as.array.mcmc.list(jsamples), 2, sd),
    loglik = NA,
    algorithm = list(type="MCMC sampling with JAGS", samples=jsamples)
  )

  rjags::unload.module("wiener")
  return(res)
}

mlehd <- function(data, jagscmd=NULL) {
  if(is.null(jagscmd)) jagscmd <- list(n.chains=1, n.adapt=1000, n.iter=10000, thin=10)

  data$id <- factor(data$id)
  I <- length(levels(data$id))
  J <- max(aggregate(y~id, data, length)$y)
  yvar <- matrix(NA, I, J)
  for(i in 1:I) {
    for(j in 1:length(data$y[data$id==levels(data$id)[i]])) yvar[i,j] <- data$y[data$id==i][j]
  }

  jinits <- list()
  for(ch in 1:jagscmd$n.chains) {
    alpha <- runif(1,0.01,2); tau <- min(data[,1])/3; beta <- runif(1,.2,.8)
    mDelta <- runif(1,-1,1); sDelta <- runif(1,0,1)
    jinits <- append(jinits, list(alpha=alpha,tau=tau,beta=beta,
                                  mDelta=mDelta,sDelta=sDelta))
  }

  jdata <- list(y=yvar, I=I, J=J)

  wdmhd.jmodel <- textConnection("model { # hierarchical wiener diffusion model

  alpha ~ dunif(0.01,100)     # k
  tau ~ dunif(0.01,1)         # k
  beta ~ dunif(0.01,1)        # k

  mDelta ~ dunif(-100,100)    # k
  sDelta ~ dunif(0.0001,100)  # k
  prDelta <- pow(sDelta, -2)

	for(i in 1:I) {              # I <- number of subjects

    delta[i] ~ dnorm(mDelta, prDelta)

    for(j in 1:J) {            # J <- number of trials
      y[i,j] ~ dwiener(alpha, tau, beta, delta[i])
    }
	}

}")

  jmodel <- rjags::jags.model(wdmhd.jmodel, jdata, jinits, jagscmd$n.chains, jagscmd$n.adapt, TRUE)
  jsamples <- rjags::coda.samples(jmodel, c("alpha", "tau", "beta", "mDelta", "sDelta"), 
                                  jagscmd$n.iter, jagscmd$thin)

  res <- list(
    coefficients = apply(coda::as.array.mcmc.list(jsamples), 2, mean),
    sd = apply(coda::as.array.mcmc.list(jsamples), 2, sd),
    loglik = NA,
    algorithm = list(type="MCMC sampling with JAGS", samples=jsamples)
  )

  rjags::unload.module("wiener")
  return(res)
}

mlehdc <- function(data, jagscmd=NULL) {
  if(is.null(jagscmd)) jagscmd <- list(n.chains=1, n.adapt=1000, n.iter=10000, thin=10)

  data$id <- factor(data$id)
  data$cond <- factor(data$cond)
  I <- length(levels(data$id))
  C <- length(levels(data$cond))
  J <- max(aggregate(y~id, data, length)$y)
  yvar <- array(NA, c(I,C,J))
  for(i in 1:I) {
    for(c in 1:C) {
      for(j in 1:length(data$y[data$id==levels(data$id)[i]])) 
        yvar[i,c,j] <- data$y[data$id==levels(data$id)[i] & data$cond==levels(data$cond)[c]][j]
    }
  }

  jinits <- list()
  for(ch in 1:jagscmd$n.chains) {
    alpha <- runif(1,0.01,2); tau <- min(data[,1])/3; beta <- runif(1,.2,.8)
    sDelta <- runif(1,0,1)
    mDelta <- vector()
    for(c in 1:C) {
      mDelta[c] <- runif(1,-1,1)
    }
    jinits <- append(jinits, list(alpha=alpha,tau=tau,beta=beta,
                                  mDelta=mDelta,sDelta=sDelta))
  }

  jdata <- list(y=yvar, I=I, J=J, C=C)

  wdmhd.jmodel <- textConnection("model { # hierarchical wiener diffusion model

  alpha ~ dunif(0.01,100)          # k
  tau ~ dunif(0.01,1)              # k
  beta ~ dunif(0.01,1)             # k
  sDelta ~ dunif(0.0001,100)       # k
  prDelta <- pow(sDelta, -2)

  for (c in 1:C) {                 # C <- number of conditions

    mDelta[c] ~ dunif(-100,100)    # k[c]

    for(i in 1:I) {                 # I <- number of subjects
      delta[i,c] ~ dnorm(mDelta[c], prDelta)
      for(j in 1:J) {               # J <- number of trials
        y[i,c,j] ~ dwiener(alpha, tau, beta, delta[i,c])
      }
    }
	}

}")

  jmodel <- rjags::jags.model(wdmhd.jmodel, jdata, jinits, jagscmd$n.chains, jagscmd$n.adapt, TRUE)
  jsamples <- rjags::coda.samples(jmodel, c("alpha", "tau", "beta", "mDelta", "sDelta"), 
                                  jagscmd$n.iter, jagscmd$thin)#, na.rm=FALSE)

  res <- list(
    coefficients = apply(coda::as.array.mcmc.list(jsamples), 2, mean),
    sd = apply(coda::as.array.mcmc.list(jsamples), 2, sd),
    loglik = NA,
    algorithm = list(type="MCMC sampling with JAGS", samples=jsamples)
  )

  rjags::unload.module("wiener")
  return(res)
}

mlerc <- function(data, jagscmd=NULL) {
  if(is.null(jagscmd)) jagscmd <- list(n.chains=1, n.adapt=1000, n.iter=10000, thin=10)

  data$id <- factor(data$id)
  I <- length(levels(data$id))
  J <- max(aggregate(y~id, data, length)$y)
  yvar <- matrix(NA, I, J)
  for(i in 1:I) {
    for(j in 1:length(data$y[data$id==levels(data$id)[i]])) yvar[i,j] <- data$y[data$id==i][j]
  }

  jinits <- list()
  for(ch in 1:jagscmd$n.chains) {
    alpha <- vector()
    mTau <- vector(); mBeta <- vector(); mDelta <- vector()
    uTau <- vector(); uBeta <- vector(); sDelta <- vector()
    for(i in 1:I) {
      alpha[i] <- runif(1,1,2)
      mTau[i] <- min(data[,1])/3; mBeta[i] <- runif(1,.2,.8); mDelta[i] <- runif(1,-1,1)
      uTau[i] <- runif(1,0,min(data[,1])/3); uBeta[i] <- runif(1,0,0.1); sDelta[i] <- runif(1,0,1)
    }
    jinits <- append(jinits, list(alpha=alpha,mTau=mTau,mBeta=mBeta,mDelta=mDelta, 
                 uTau=uTau,uBeta=uBeta,sDelta=sDelta))
  }

  jdata <- list(y=yvar, I=I, J=J)

  ratcliff.jmodel <- textConnection("model { # ratcliff diffusion model

	for(i in 1:I) {                  # I <- number of subjects

    alpha[i] ~ dunif(0.01,100)     # k=1
    mDelta[i] ~ dunif(-100,100)    # k=2
    sDelta[i] ~ dunif(0.0001,100)  # k=3
    prDelta[i] <- pow(sDelta[i], -2)
    mBeta[i] ~ dunif(0.01,1)       # k=4
    uBeta[i] ~ dunif(0.0001,1)     # k=5
    mTau[i] ~ dunif(0.01,1)        # k=6
    uTau[i] ~ dunif(0.0001,1)      # k=7

    mDrift[i] ~ dnorm(mDelta[i], prDelta[i])
    # sDrift <- 1 # fixed to 1 in dwiener code

    for(j in 1:J) {                # J <- number of trials
      tau[i,j] ~ dunif(mTau[i]-uTau[i]/2, mTau[i]+uTau[i]/2)
      beta[i,j] ~ dunif(mBeta[i]-uBeta[i]/2, mBeta[i]+uBeta[i]/2)
      y[i,j] ~ dwiener(alpha[i], tau[i,j], beta[i,j], mDrift[i])
    }
	}

}")

  jmodel <- rjags::jags.model(ratcliff.jmodel, jdata, jinits, jagscmd$n.chains, jagscmd$n.adapt, TRUE)
  jsamples <- rjags::coda.samples(jmodel, c("alpha", "mTau", "uTau", "mBeta", "uBeta", "mDelta", "sDelta"), 
                                  jagscmd$n.iter, jagscmd$thin)

  res <- list(
    coefficients = apply(coda::as.array.mcmc.list(jsamples), 2, mean),
    sd = apply(coda::as.array.mcmc.list(jsamples), 2, sd),
    loglik = NA,
    algorithm = list(type="MCMC sampling with JAGS", samples=jsamples)
  )

  rjags::unload.module("wiener")
  return(res)
}
