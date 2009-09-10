\name{marg.likelihood}
\alias{marg.likelihood}

\title{Calculating the marginal likelihood}

\description{
Calculating the marginal likelihood.
}

\usage{
marg.likelihood(penden.env,pen.likelihood)
}

\arguments{
  \item{penden.env}{Containing all information, environment of pendensity()}
  \item{pen.likelihood}{penalized log likelihood}
}
\details{
Calculating is done using a Laplace approximation for the integral of the marginal likelihood.

The needed values are saved in the environment.
}
\value{
Returning the marginal likelihood.
}
\references{Penalized Density Estimation, Kauermann G. and Schellhase C. (2009), to appear.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\keyword{math}

