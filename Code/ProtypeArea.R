wSize <- sum(community[nrow(community),2:length(community)] > 10^-15)
r <- list()

for(i in 1:10) {
  z <- c()
  for(j in 1:length(d)){
    z <- c(z, lengths(d[[j]]["Living"], use.names=FALSE))
  }
  z <- data.frame(z)
  colnames(z) <- "y"
  z["x"] <- c(1:nrow(z))
  stepSize <- 1
  z <- z[seq(1, nrow(z), stepSize), ]
}
stepPlot <- ggplot(data=z, aes(x=x, y=y)) +
  geom_step() +
  geom_point() +
  xlab("Step Number") +
  ylab("Nisle") +
  expand_limits(y=c(0,wSize))