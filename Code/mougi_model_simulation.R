library(deSolve)
library(lattice)


  ## make some random, cascade, and niche food webs
  S = 40     ## set species richness
  C = 0.2   ## set connectance
  N = 1     ## set the number of replicate webs to make
  L = S^2*C  ## calculate number of links from S and C
  
  
  xxx <- Cascade.model(S, L, N)
  
  
  # par(mfrow = c(1,2))
  # Plot.matrix(xxx)
  # box()
  
  # Number of Species
  n <- S
  r <- runif(n, -1,1)
  s <- runif(n, 1,1)
  g <- runif(n)
  # a <- matrix(runif(n*n, 0, 0.1),nrow=n)
  a <- xxx * matrix(runif(n*n, 0,1),nrow=n)
  diag(a) <- rep(0,n)
  
  # init.x <- rep(1, n)
  init.x <- runif(n)
  
  mougi_model <- function(t,x,parms){
    dx <- x * (r - s*x + g * (a %*% x) - (t(a) %*% x))
    list(dx)
  }
  
  n.integrate <- function(time=time, init.x= init.x, model=model){
    t.out <- seq(time$start,time$end,length=time$steps)
    as.data.frame(lsoda(init.x, t.out, model, parms = parms))
  }
  
  
  
  
  
  # Integration window
  time <- list(start = 0, end = 100, steps = 100)
  # dummy variable for lvm() function defined above
  parms <- c(0) ### dummy variable (can have any numerical value)
  
  
  out <- n.integrate(time, init.x, model = mougi_model)
  
  
  plot(out[,1], out[,2], type='l', ylim = c(0, max(out[,2:ncol(out)])), lwd = 2)
  for (i in 3:ncol(out)){
    points(out[,1], out[,i], type='l', col = i-1, lwd = 2)
  }

  
  
  sum(out[nrow(out),2:n] > 10^-5) / n
  
  
  