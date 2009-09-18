### Name: Allianz
### Title: Daily final prices (DAX) of the German stock Allianz in the
###   years 2006 and 2007
### Aliases: Allianz
### Keywords: datasets

### ** Examples

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



