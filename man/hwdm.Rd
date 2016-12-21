\name{hwdm}
\alias{hwdm}
\title{Hierarchical Bayesian Wiener diffusion model fit functions}
\description{
  \code{hwdm} creates parameter estimates for the parameters
  of the underlying hierarchical Wiener diffusion model used.
}
\usage{
  hwdm(data, jagscmd=NULL, model="hwdm", idvar="id", condvar=NULL)
}
\arguments{
  \item{data}{is the data object containing data coming from a
  (hypothetical) Wiener Diffusion process. For further details on the data
  object, see \link{is.wiener}.}
  \item{jagscmd}{a list with options for the jags call. Defaults to
    \code{NULL} which is equal to: \code{list(n.chains=1, n.adapt=1000,
    n.iter=10000, thin=10)}.}
  \item{model}{the model to use. Options are the following:
    \code{hwdm} for a simple hierarchical extension to the standard
    model, with a normal distribution assumed for the \code{delta}
    parameter and a uniform distribution for the remaining 3 parameters;
    \code{wdmhd} for a hierarchical extension only on the \code{delta}
    parameter by assuming that this parameter stems from a normal
    distribution;
    \code{wdmhdc} same as the \code{wdmhd} extension, but allowing for
    different means of the \code{delta} parameter depending on
    groups/conditions;
    \code{rdm} provides a hierarchical extension in the form of the
    Ratcliff diffusion model with 7 parameters: 
    \code{alpha}, \code{mDelta}, \code{sDelta},
    \code{mBeta}, \code{uBeta}, \code{mTau}, \code{uTau}.}
  \item{idvar}{variable/factor name defining the individual/person in data
    in data (typically column name in data.frame), defaults to \code{id}.}
  \item{condvar}{variable/factor name defining the condition in data
    (typically column name in data.frame), if necessary for model.}
}
\details{
  The \code{hwdm} function calls an mcmc based estimation
  routine (JAGS with the JAGS Wiener module), 
  to estimate the model parameters of the specified model. 
  }
%\references{} 
%\note{}
%\author{}
\examples{

\dontrun{
## generate random data
dat <- rbind(cbind(rwiener(100, 2,.3,.5,1), group="A"),
        cbind(rwiener(100,2,.3,.5,-1), group="B"))

## hierarchical models
hwdm1 <- hwdm(dat, model="hwdm")
wdmhd1 <- hwdm(dat, model="wdmhd")
}

}
\keyword{hwdm}
