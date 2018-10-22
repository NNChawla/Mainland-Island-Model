library(deSolve)
library(lattice)
library(tidyverse)

CvNs <- function(S, C, step) {

  ## make some random, cascade, and niche food webs
  ## S is species richness
  ## C is connectance
  N <- 1     ## set the number of replicate webs to make
  matrixSize <- round(C/step)
  replicates <- 10
  
  communities <- rep(list(rep(list(vector("list", matrixSize)), matrixSize)), replicates)
  persistences <- rep(list(rep(list(vector("list", matrixSize)), matrixSize)), replicates)
  means <- rep(list(rep(list(0), matrixSize)), matrixSize)
  stdDevs <- rep(list(rep(list(vector("integer", replicates)), matrixSize)), matrixSize)

  for (CvN in 1:replicates) {
    C_step <- step
    S_step <- S/(C/step)
    
    for(rowC in 1:matrixSize) {
      
      # if (C_step > 0.5*(1-(1/S))) {
      #     C_step <- C_step + step
      #     next
      # }
      
      for(colN in 1:matrixSize) {
        #print(c(C_step, S_step))
        
        L <- round(S_step^2*C_step)  ## calculate number of links from S and C
        
        # skipping community generation because formula is less than 0 which throws error
        if (((S_step^2 - S_step)/2 - L) < 0) {
          #print(c(S_step, C_step, "error"))
          S_step <- S_step + S/(C/step)
          next
        }
        
        #generating model based on parameters
        #print(c(S_step, C_step))
        xxx <- Cascade.model(S_step, L, N)
        
        # Number of Species
        n <- S_step
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
        communities[[CvN]][[rowC]][[colN]] <- out
        persistences[[CvN]][[rowC]][[colN]] <- mean(out[nrow(out),2:n+1] > 10^-5)
        
        S_step <- S_step + S/(C/step)
      }
      
      S_step <- S/(C/step)
      C_step <- C_step + step
      
    }
  
  }
  
  for(CvN in 1:replicates) {
      
    for(rowC in 1:matrixSize) {
      
      for(colN in 1:matrixSize) {
        
        #print(c(CvN, rowC, colN))
        means[[rowC]][[colN]] <- means[[rowC]][[colN]] + persistences[[CvN]][[rowC]][[colN]]
        if(!is.null(persistences[[CvN]][[rowC]][[colN]]))
          stdDevs[[rowC]][[colN]][[CvN]] <- persistences[[CvN]][[rowC]][[colN]]
      }
    }
  }
  
  for(rowC in 1:matrixSize) {
    for(colN in 1:matrixSize) {
      stdDevs[[rowC]][[colN]] <- sd(stdDevs[[rowC]][[colN]])
      #print(c(rowC, colN, stdDevs[[rowC]][[colN]]))
    }
  }
  
  speciesLabels <- c(seq((S/(C/step)), S, (S/(C/step))))
  connectanceLabels <- c(seq(step, C, step))
  stdDevs <- as.data.frame(matrix(unlist(stdDevs), nrow=length(unlist(stdDevs[1]))))
  rownames(stdDevs) <- connectanceLabels
  colnames(stdDevs) <- speciesLabels
  
  matricies <- list("communities" = communities, "persistences" = persistences, "mean" = means, "stdDev" = stdDevs)
  return(matricies)

}

plotGraph <- function(dataFrame, graphType) {
    ggplot(melt(as.matrix(dataFrame), id.vars = c("X1", "X2", "value"))) +
    scale_fill_gradient(low = "steelblue", high = "white") +
    ylab("Connectance") +
    xlab("Species") +
    geom_tile(mapping = aes(x = X2, y = X1, fill = value)) +
    labs(fill = graphType)
}