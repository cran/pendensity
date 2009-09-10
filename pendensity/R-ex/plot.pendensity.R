### Name: plot.pendensity
### Title: Plotting estimated penalized densities
### Aliases: plot.pendensity plot
### Keywords: aplot

### ** Examples

y <- rnorm(100)
test <- pendensity(y~1)
plot(test)

#empirical distribution
plot(test,plot.val=2)



