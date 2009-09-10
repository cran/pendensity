\name{L.mat}
\alias{L.mat}
\title{Calculates the difference matrix of order m}
\description{
Calculating the differences matrix 'L' of order 'm', depending on the number of knots 'k'.
}
\usage{
L.mat(penden.env)
}

\arguments{
  \item{penden.env}{Environment of pendensity()}
}
\value{
Returns the difference matrix of order 'm' for given number of knots 'K'.

The needed values are saved in the environment.
}
\references{Penalized Density Estimation, Kauermann G. and Schellhase C. (2009), to appear.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\note{
Right now, the difference matrix is implemented for m=1,2,3,4.
}
\keyword{nonparametric}

