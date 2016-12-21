\name{bwdm}
\alias{bwdm}
\title{Bayesian Wiener diffusion model fit functions}
\description{
  \code{bwdm} creates parameter estimates for the four parameters
  of the Wiener diffusion model.
}
\usage{
  bwdm(data, yvar=c("q", "resp"), alpha=NULL, tau=NULL, beta=NULL, delta=NULL,
               xvar=NULL, start=NULL, fixed=0, jagscmd=NULL, idvar="id")
}
\arguments{
  \item{data}{is the data object containing data coming from a
  (hypothetical) Wiener Diffusion process. For further details on the data
  object, see \link{is.wiener}.}
  \item{yvar}{represents an optional vector, that can be used to define the
  names of the reaction time column. For further details on the data
  object, see \link{is.wiener}.}
  \item{alpha}{optional, can be used to fix the alpha parameter to the
  given value.}
  \item{tau}{optional, can be used to fix the tau parameter to the given
  value.}
  \item{beta}{optional, can be used to fix the beta parameter to the given
  value.}
  \item{delta}{optional, can be used to fix the delta parameter to the
  given value.}
  \item{xvar}{optional: group factor variable to estimate all unfixed
  parameters separate for the given groups.}
  \item{start}{an optional vector with the four starting parameter values
  in the following order: alpha, tau, beta, delta.}  
  \item{fixed}{a number indicating how many of the parameters are fixed
    (not free). This number will be subtracted from the number of free
    parameters. Defaults to 0.}  
  \item{jagscmd}{a list with options for the jags call. Defaults to
    \code{NULL} which is equal to: \code{list(n.chains=1, n.adapt=1000,
    n.iter=10000, thin=10)}.}
  \item{idvar}{variable/factor name defining the individual/person in data
    in data (typically column name in data.frame), defaults to \code{id}.
    Can be omitted if irrelevant.}
}
\details{
  The \code{bwdm} function calls an mcmc based estimation
  routine (JAGS with the JAGS Wiener module), 
  to estimate the model parameters. 
  }
%\references{} 
%\note{}
%\author{}
\examples{

\dontrun{
## generate random data
dat <- rbind(cbind(rwiener(100, 2,.3,.5,1), group="A"),
        cbind(rwiener(100,2,.3,.5,-1), group="B"))

## fit wdm
bwdm1 <- bwdm(dat)

## extract parameters
coef(bwdm1)

## further models
bwdm2 <- bwdm(dat, beta=.5)
bwdm3 <- bwdm(dat, alpha=bwdm1$coefficients["alpha"], 
              tau=bwdm1$coefficients["tau"],
              beta=bwdm1$coefficients["beta"], 
              xvar="group")
}

}
\keyword{bwdm}
