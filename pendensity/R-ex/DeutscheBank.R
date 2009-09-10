### Name: DeutscheBank
### Title: Daily final prices (DAX) of the German stock Deutsche Bank in
###   the years 2006 and 2007
### Aliases: DeutscheBank
### Keywords: datasets

### ** Examples

data(DeutscheBank)

form<-'%d.%m.%y %H:%M'

time.DeutscheBank <- strptime(DeutscheBank[,1],form)

#looking for all dates in 2006
data.DeutscheBank <- DeutscheBank[which(time.DeutscheBank$year==106),2]

#building differences of first order
DeutscheBank1 <- c()
for(i in 2:length(data.DeutscheBank)) DeutscheBank1[i-1] <- data.DeutscheBank[i]-data.DeutscheBank[i-1]

#estimating the density
density.DeutscheBank <- pendensity(DeutscheBank1~1)
plot(density.DeutscheBank)



