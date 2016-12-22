#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# rrdm.R
#
# (c) 2016 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2016-12-22
# last mod 2016-12-22 18:02 DW
#

rrdm <- function(n, alpha=2, mDelta=0, sDelta=1, mBeta=0.5, uBeta=0.1, mTau=0.3, uTau=.1) {
  rjags::load.module("wiener", quiet=TRUE)

  ratcliff.jmodel <- textConnection("model { # ratcliff diffusion model

	for(i in 1:I) {                  # I <- number of subjects

    prDelta[i] <- pow(sDelta[i], -2)

    mDrift[i] ~ dnorm(mDelta[i], prDelta[i])
    # sDrift <- 1 # fixed to 1 in dwiener code

    for(j in 1:J) {                # J <- number of trials
      tau[i,j] ~ dunif(mTau[i]-uTau[i]/2, mTau[i]+uTau[i]/2)
      beta[i,j] ~ dunif(mBeta[i]-uBeta[i]/2, mBeta[i]+uBeta[i]/2)
      y ~ dwiener(alpha[i], tau[i,j], beta[i,j], mDrift[i])
    }
	}

}")

  jdata <- list(I=1, J=1, # I and J are fixed to 1 for sampling
    alpha=alpha,
    mDelta=mDelta, sDelta=sDelta,
    mBeta=mBeta, uBeta=uBeta,
    mTau=mTau, uTau=uTau)

  jmodel <- rjags::jags.model(ratcliff.jmodel, jdata, n.chains=1, n.adapt=0, quiet=TRUE)
  jsamples <- rjags::coda.samples(jmodel, c("y"), n, 1)

  res <- unname(coda::as.array.mcmc.list(jsamples))
  res <- RWiener::revamp.numdata.wiener(res)

  rjags::unload.module("wiener", quiet=TRUE)
  return(res)
}
