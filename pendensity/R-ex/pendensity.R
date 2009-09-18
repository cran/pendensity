### Name: pendensity
### Title: Calculating penalized density
### Aliases: pendensity
### Keywords: nonparametric

### ** Examples

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

form<-'%d.%m.%y %H:%M'

time.Allianz <- strptime(Allianz[,1],form)

#looking for all dates in 2006
data.Allianz <- Allianz[which(time.Allianz$year==106),2]

#building differences of first order
Allianz1 <- c()
for(i in 2:length(data.Allianz)) Allianz1[i-1] <- data.Allianz[i]-data.Allianz[i-1]

#estimating the density
density.Allianz <- pendensity(Allianz1~1)
plot(density.Allianz)

#################

#density-example of the stock exchange Allianz in 2006 and 2007
data(Allianz)

form<-'%d.%m.%y %H:%M'

time.Allianz <- strptime(Allianz[,1],form)

#looking for all dates in 2006
data.Allianz <- Allianz[which(time.Allianz$year==106|time.Allianz$year==107),2]

#building differences of first order
Allianz1 <- c()
for(i in 2:length(data.Allianz)) Allianz1[i-1] <- data.Allianz[i]-data.Allianz[i-1]

#estimating the density
density.Allianz <- pendensity(Allianz1~as.factor(time.Allianz$year))
plot(density.Allianz,legend.txt=c("2006","2007"))



