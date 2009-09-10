\name{pendensity}
\alias{pendensity}

\title{Calculating penalized density}
\description{
Main program for estimation penalized densities. The estimation can be done for response with or without any covariates. The covariates have to be factors. The response is called 'y', the covariates 'x'. We estimate densities using penalized splines. This done by using a number of knots and a penalty parameter, which are sufficient large. We penalize the m-order differences of the beta-coefficients to estimate the weights 'ck' of the used base functions.
}
\usage{
pendensity(form, base = "bspline", no.base = NULL, max.iter = 20,
 lambda0 = 50000, q = 3, plot.bsp = FALSE, sort = TRUE, with.border = NULL, m = q)
}

\arguments{
  \item{form}{formula describing the density, the formula is\cr
	\deqn{y \sim x_1 + x_2 + \dots x_n}{%
	\eqn{y ~ x1 + x2 + \dots + xn}}
   where 'x' have to be factors.}
  \item{base}{supported bases are "bspline" or "gaussian"}
  \item{no.base}{how many knots 'K', following the approach to use 2'no.base'+1 knots, if 'no.base' is NULL, default is K=41.}
  \item{max.iter}{maximum number of iteration, the default is max.iter=20.}
  \item{lambda0}{start penalty parameter, the default is lambda0=50000}
  \item{q}{order of B-Spline base, the default is 'q=3'}
  \item{plot.bsp}{TRUE or FALSE, if TRUE the used B-Spline base should be plotted}
  \item{sort}{TRUE or FALSE, if TRUE the response and the covariates should be sorted}
  \item{with.border}{determining the number of additional knots on the left and the right of the support of the response. The number of knots 'no.base' is not influenced by this parameter. The amount of knots 'no.base' are placed on the support of the response. The amount of knots determined in 'with.border' is placed outside the support and reduce the amount of knots on the support about its value.}
  \item{m}{m-th order difference for penalization. Default is m=3.}
}
\details{
pendensity() begins with setting the parameters for the estimation. Checking the formula and transferring the data into the program, setting the knots and creating the base, depending on the chosen parameter 'base'. Moreover the penalty matrix is constructed. At the beginning of the first iteration the beta parameter are set equal to zero. With this setup, the first log likelihood is calculated and is used for the first iteration for a new beta parameter.\cr
The iteration for a new beta parameter is done with a Newton-Raphson-Iteration and implemented in the function 'new.beta.val'. We calculate the direction of the Newton Raphson step for the known \eqn{beta_t} and iterate a step size bisection to control the maximizing of the penalized likelihood \deqn{l_p(\beta_t,\lambda_0)}{%
\eqn{l(beta,lambda0)}}. This means we set \deqn{\beta_{t+1}=\beta_t - 2^{-v} \{s_p(\beta_t,\lambda_0) \cdot (-J_p(\beta_t,\lambda_0))^{-1}\}}{%
\eqn{beta[t+1]=beta[t]-(2/v)*sp(beta,lambda0)*(-Jp(beta[t],lambda0))^-1}} with \eqn{s_p} as penalized first order derivative and \eqn{J_p} as penalized second order derivative. We begin with \eqn{v=0}. Not yielding a new maximum for a current v, we increase v step by step respectively bisect the step size. We terminate the iteration, if the step size is smaller than some reference value epsilon (eps=1e-3) without yielding a new maximum. We iterate for new parameter beta until the new log likelihood depending on the new estimated parameter beta differ less than 0.1 log-likelihood points from the log likelihood estimated before.\cr
After reaching the new parameter beta, we iterate for a new penalty parameter lambda. This iteration is done by the function 'new.lambda'. The iteration formula is \deqn{ {\lambda}^{-1}=\frac{\hat{\beta}^T D_m \hat{\beta}}{df(\hat{\lambda})-p(m-1)}}{%
\eqn{lambda^-1=beta^T Dm beta / (df(lambda)-p(m-1))}.}\cr
The iteration for the new lambda is terminated, if the approximate degree of freedom minus p*(m-1) is smaller than some epsilon2 (eps2=0.01). Moreover, we terminate the iteration if the new lambda is approximatively converted, i.e. the new lambda differs only 0.001*old lambda (*) from the old lambda. If these both criteria doesn't fit, the lambda iteration is terminated after eleven iterations. \cr
We begin a new iteration with the new lambda, restarting with parameter
  beta setting equal to zero again. This procedure is repeated until
  convergence of lambda, i.e. that the new lambda fulfills the criteria
  (*). If this criteria is not fulfilled after 20 iterations, the total iteration terminates.\cr
After terminating all iterations, the final AIC, ck and beta are saved in the output.\cr
For speediness, all values, matrices, vectors etc. are saved in an environment called 'penden.env'. Most of the used programs get only this environment as input.

}
\value{
Returning an object of class pendensity.
The class pendensity consists of the "call" and three main groups "values", "splines" and "results".

call: the formula prompted for calculation of the penalized density.

#####

\$values contains:
y: the values of the response variable
x: the values of the covariate(s)
sort: TRUE/FALSE if TRUE the response (and covariates) have been sorted in increasing order of the response

\$values\$covariate contains 
Z: matrix Z
levels: existing levels of each covariate
how.levels: number of existing levels of all covariates
how.combi: number of combination of levels
x.factor: list of all combination of levels

#####

\$splines contains 
K: number of knots
N: number of coefficients estimated for each base, depending on the number of covariates.
MeanW: values of the knots used for splines and means of the Gaussian densities
Stand.abw: values of the standard deviance of the Gaussian densities
h: distance between the equidistant knots
m: used difference order for penalization
q: used order of the B-Spline base
stand.num: calculated values for standardization getting B-Spline densities
base: used kind of base, "bspline" or "gaussian"
base.den: values of the base of order q created with knots=knots.val$val
base.den2: values of the base of order q+1 created with knots=knots.val$val. Used for calculating the distribution function(s).
L: used difference matrix
Dm: used penalty matrix, depending on lambda0, L (,Z) and n=number of observations
help.degree: additional degree(s) depending on the number of knots K and the used order q

\$splines\$knots.val contains:
val: list of the used knots in the support of the response
all: list of the used knots extended with additional knots used for constructing

#####

results contains:
ck: final calculated weights ck
beta.val: final calculated parameter beta
lambda0: final calculated lambda0
f.hat.val: final fitted density values of the response y
variance.par: final variances of the parameter beta
bias.par: final bias of the parameter beta

results\$AIC contains:
my.AIC: final AIC value
my.trace: trace component of the final AIC

results\$Derv contains:
Derv2.pen: final penalized second order derivation 
Derv2.cal: final non-penalized second order derivation 
Derv1.cal: final non-penalized first order derivation 

results\$iterations contains:
list.opt.results: list of the final results of each iteration of new beta + new lambda
all.lists: list of lists. Each list contains the result of one iteration

results\$likelihood contains:
pen.Likelihood: final penalized log likelihood
opt.Likelihood: final log likelihood
marg.Likelihood: final marginal likelihood
}

\references{Penalized Density Estimation, Kauermann G. and Schellhase C. (2009), to appear.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\examples{
#first simple example

set.seed(27)

y <- rnorm(100)
test <- pendensity(y~1)

#plotting the estimated density
plot(test)

#expand the support at the boundary
test2 <- pendensity(y~1,with.border=8)
plot(test2)

#expand the support at the boundary and enlarge the number of knots to get the same number of knots in the support
test3 <- pendensity(y~1,with.border=8,no.base=28)
plot(test3)

test4 <- pendensity(y~1,with.border=10,no.base=35)
plot(test4)

#################

#second simple example
#with covariate

x <- rep(c(0,1),200)
y <- rnorm(400,x*0.2,1)
test <- pendensity(y~as.factor(x))
plot(test)


#################

#density-example of the stock exchange Allianz in 2006
data(Allianz)

form<-'\%d.\%m.\%y'

time.Allianz <- strptime(Allianz[,1],form)

#looking for all dates in 2006
data.Allianz <- Allianz[which(time.Allianz$year==106),2]

#building differences of first order
d.Allianz <- diff(data.Allianz)

#estimating the density
density.Allianz <- pendensity(d.Allianz~1)
plot(density.Allianz)

#################

#density-example of the stock exchange Allianz in 2006 and 2007
data(Allianz)

form<-'\%d.\%m.\%y'

time.Allianz <- strptime(Allianz[,1],form)

#looking for all dates in 2006
data.Allianz <- Allianz[which(time.Allianz$year==106|time.Allianz$year==107),2]

#building differences of first order
d.Allianz <- diff(data.Allianz)

#estimating the density
density.Allianz <- pendensity(d.Allianz~as.factor(time.Allianz$year))
plot(density.Allianz,legend.txt=c("2006","2007"))
}

\keyword{nonparametric}