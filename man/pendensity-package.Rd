\name{pendensity-package}
\alias{pendensity-package}
\docType{package}
\title{
The package 'pendensity' offers routines for estimating penalized unconditional and conditional (on factor groups) densities.
}
\description{
The package 'pendensity' offers routines for estimating penalized
unconditional and conditional (on factor groups) densities. For details
see the description of the function pendensity.
}
\details{
\tabular{ll}{
Package: \tab pendensity\cr
Type: \tab Package\cr
Version: \tab 0.1.1\cr
Date: \tab 2009-09-10\cr
License: GPL (>= 2)
LazyLoad: \tab yes\cr
}
The packages contributes the function 'pendensity()' for estimating densities using penalized splines techniques. 
}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\references{
Penalized Density Estimation, Kauermann G. and Schellhase C. (2009), to appear.
}
\keyword{nonparametric}
\examples{
y <- rnorm(100)
test <- pendensity(y~1)
plot(test)

}
