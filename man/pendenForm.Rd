\name{pendenForm}
\alias{pendenForm}
\alias{string.help}
\title{Formula interpretation and data transfer}
\description{
Function 'pendenForm' interprets the input 'form' of the function pendensity(),transfers the response and covariates data back to the main program and categorize the values to groupings.
}
\usage{
pendenForm(form)
string.help(string, star = " ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{form}{formula from the input of pendensity()}
  \item{string}{string of the formula}
  \item{star}{separating letter}
}

\value{
Returning the values and the structure of response and covariates.
}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\keyword{IO}
