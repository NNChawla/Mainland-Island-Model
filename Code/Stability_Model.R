library(deSolve)
library(lattice)
library(tidyverse)
library(reshape)

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
      means[[rowC]][[colN]] <- means[[rowC]][[colN]]/replicates
      stdDevs[[rowC]][[colN]] <- sd(stdDevs[[rowC]][[colN]])
    }
  }
  
  speciesLabels <- c(seq((S/(C/step)), S, (S/(C/step))))
  connectanceLabels <- c(seq(step, C, step))
  stdDevs <- as.data.frame(matrix(unlist(stdDevs), nrow=length(unlist(stdDevs[1]))))
  rownames(stdDevs) <- connectanceLabels
  colnames(stdDevs) <- speciesLabels
  
  means <- as.data.frame(matrix(unlist(means), nrow=length(unlist(means[1]))))
  rownames(means) <- connectanceLabels
  colnames(means) <- speciesLabels[1:length(means)]
  
  #removing null columns from sd dataframe
  stdDevs <- stdDevs[1:length(stdDevs), 1:length(means)]
  
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

subsetPath <- function(community, numSpecies) {
  randSpecies <- sample(2:length(community), numSpecies, replace=FALSE)
  subsetData <- data.frame(matrix(nrow=100, ncol=numSpecies))
  for(i in 1:numSpecies) {
    subsetData[i] <- community[[randSpecies[i]]]
  }
  return(list(randSpecies, subsetData))
}

nstarGraph <- function(meanData, Nstar, interval=0.5) {
  nstarMatrix <- data.frame(matrix(nrow=nrow(meanData), ncol=ncol(meanData)))
  for(i in 1:ncol(nstarMatrix)){
    nstarMatrix[i] <- meanData[,i]*as.numeric(colnames(meanData)[i])
  }
  #Use indicies to get CVN pairs along matrix and sample them using island assembly
  indicies <- which(nstarMatrix > Nstar-interval & nstarMatrix < Nstar+interval, arr.ind = TRUE)
  
  
}