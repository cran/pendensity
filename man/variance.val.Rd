\name{variance.val}
\alias{variance.val}
\title{Calculating variance and standard deviance of each observation.}
\description{
 Calculating the variance and standard deviance of each observation. Therefore we use the variance of the parameter set beta, called 'var.par'.
}
\usage{
variance.val(base.den, var.par, weight, K, x, list.len, Z, x.factor, y.val=NULL)
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
  \item{y.val}{optimal values for calculating the variance in any point
    yi in the case of a factorial density}
}

\value{
Returning a vector with the standard deviance of each observation.
}

\references{Density Estimation with a Penalized Mixture Approach, Kauermann G. and Schellhase C. (2009), to appear.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld>}

\keyword{nonparametric}
