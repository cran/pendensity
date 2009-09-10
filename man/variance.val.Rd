\name{variance.val}
\alias{variance.val}
\title{Calculating variance and standard deviance of each observation.}
\description{
 Calculating the variance and standard deviance of each observation. Therefore we use the variance of the parameter set beta, called 'var.par'.
}
\usage{
variance.val(base.den, var.par, weight, K, x, list.len, Z, x.factor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{base.den}{base values}
  \item{var.par}{variance of the parameter set beta}
  \item{weight}{weights ck}
  \item{K}{number of knots}
  \item{x}{covariates}
  \item{list.len}{number of covariate combinations}
  \item{Z}{covariate matrix}
  \item{x.factor}{list of covariate combinations}
}

\value{
Returning a vector with the standard deviance of each observation.
}

\keyword{nonparametric}
