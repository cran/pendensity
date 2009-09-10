\name{test.equal}
\alias{test.equal}

\title{Testing pairwise equality of densities}
\description{
Every group of factor combination is tested pairwise for equality to all other groups.
}
\usage{
test.equal(obj)
}

\arguments{
  \item{obj}{object of class pendensity}
}
\details{
We consider the distribution of the integrated B-Spline density base. This is saved in the program in the object named 'mat1'. Moreover, we use the variance 'var.par' of the estimation, the weights and some matrices 'C' of the two compared densities to construct the matrix 'W'. We simulate the distribution of the test statistic using a spectral decomposition of W.
}
\value{
Returning a list of p-values for testing pairwise for equality.
}
\references{Penalized Density Estimation, Kauermann G. and Schellhase C. (2009), to appear.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld>}
\keyword{nonparametric}
