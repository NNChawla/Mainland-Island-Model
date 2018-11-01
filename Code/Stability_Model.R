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
  
  #initializing containers to be returned
  communities <- rep(list(rep(list(vector("list", matrixSize)), matrixSize)), replicates)
  persistences <- rep(list(rep(list(vector("list", matrixSize)), matrixSize)), replicates)
  means <- rep(list(rep(list(0), matrixSize)), matrixSize)
  stdDevs <- rep(list(rep(list(vector("integer", replicates)), matrixSize)), matrixSize)

  #creating a list of lists of lists of dimensions N x C x S
  #C x S throws errors if it's not square
  for (CvN in 1:replicates) {
    C_step <- step
    S_step <- S/(C/step)
    
    for(rowC in 1:matrixSize) {
      
      for(colN in 1:matrixSize) {
        
        L <- round(S_step^2*C_step)  ## calculate number of links from S and C
        
        # skipping community generation because formula is less than 0 which throws error
        if (((S_step^2 - S_step)/2 - L) < 0) {
          S_step <- S_step + S/(C/step)
          next
        }
        
        #generating model based on parameters
        xxx <- Cascade.model(S_step, L, N)
        
        # Number of Species
        n <- S_step
        r <- runif(n, -1,1)
        s <- runif(n, 1,1)
        g <- runif(n)
        a <- xxx * matrix(runif(n*n, 0,1),nrow=n)
        diag(a) <- rep(0,n)
        
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
  
  #calculating mean and sd of persistences across all N matricies of C x S
  for(CvN in 1:replicates) {
      
    for(rowC in 1:matrixSize) {
      
      for(colN in 1:matrixSize) {
        
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
  
  #assigning labels to means and sd data.frames for easier reading
  speciesLabels <- c(seq((S/(C/step)), S, (S/(C/step))))
  connectanceLabels <- c(seq(step, C, step))
  stdDevs <- as.data.frame(matrix(unlist(stdDevs), nrow=length(unlist(stdDevs[1]))))
  rownames(stdDevs) <- connectanceLabels
  colnames(stdDevs) <- speciesLabels
  
  means <- as.data.frame(matrix(unlist(means), nrow=length(unlist(means[1]))))
  rownames(means) <- connectanceLabels
  colnames(means) <- speciesLabels[1:length(means)]
  
  #removing null columns (CvN pairs past max L) from sd dataframe
  stdDevs <- stdDevs[1:length(stdDevs), 1:length(means)]
  
  #returning containers
  matricies <- list("communities" = communities, "persistences" = persistences, "mean" = means, "stdDev" = stdDevs)
  return(matricies)

}

#graphs mean or sd as heatmap
plotGraph <- function(dataFrame, graphType, xax="Species", yax="Connectance") {
  ggplot(melt(as.matrix(dataFrame), id.vars = c("X1", "X2", "value"))) +
  scale_fill_gradient(low = "steelblue", high = "white") +
  ylab(yax) +
  xlab(xax) +
  geom_tile(mapping = aes(x = X2, y = X1, fill = value)) +
  labs(fill = graphType)
}

#returns subset of community integrated through time using original final densities
subsetPath <- function(community, numSpecies, C) {
  
  #selecting a subset of species from the ones that survived initial integration
  speciesSurvived <- community[nrow(community),2:length(community)] > 10^-5
  persistingSpecies <- c(1:sum(speciesSurvived == TRUE))
  counter = 1
  for(i in 1:length(speciesSurvived)){
    if(speciesSurvived[i]) {
      persistingSpecies[counter] <- i+1
      counter <- counter + 1
    }
  }
  if(numSpecies > length(persistingSpecies))
    return(print("That many species didn't survive"))
  
  randSpecies <- sample(persistingSpecies, numSpecies, replace=FALSE)
  speciesIndicies <- c(1:numSpecies+1)
  speciesIndicies[1] <- "time"
  for(i in 1:numSpecies) {
    speciesIndicies[i+1] <- randSpecies[i]-1
    randSpecies[i] <- community[1, randSpecies[i]]
  }
  
  #integration function
  L <- round(numSpecies^2*C)  ## calculate number of links from S and C
  
  xxx <- Cascade.model(numSpecies, L, N)
  
  n <- numSpecies
  r <- runif(n, -1,1)
  s <- runif(n, 1,1)
  g <- runif(n)
  a <- xxx * matrix(runif(n*n, 0,1),nrow=n)
  diag(a) <- rep(0,n)
  
  init.x <- randSpecies
  
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
  
  
  subsetData <- n.integrate(time, init.x, model = mougi_model)
  
  colnames(subsetData) <- speciesIndicies
  return(subsetData)
}

nStarGraph <- function(container, Nstar, interval=0.5, nI = 5) {
  #creating matrix of NStar for all CvNs
  meanData <- container$mean
  community <- container$communities
  nstarMatrix <- data.frame(matrix(nrow=nrow(meanData), ncol=ncol(meanData)))
  for(i in 1:ncol(nstarMatrix)){
    nstarMatrix[i] <- meanData[,i]*as.numeric(colnames(meanData)[i])
  }
  colnames(nstarMatrix) <- colnames(meanData)
  rownames(nstarMatrix) <- rownames(meanData)
  print(plotGraph(nstarMatrix, "N*", xax="Species", yax="Connectance"))
  
  #finding the communities that have the desired Nstar within the set interval
  indicies <- which(nstarMatrix > Nstar-interval & nstarMatrix < Nstar+interval, arr.ind = TRUE)
  indicies <- as.data.frame(indicies)
  colnames(indicies) <- c("Connectance", "Species")
  colnames(nstarMatrix) <- colnames(meanData)
  rownames(nstarMatrix) <- rownames(meanData)
  indicies <- indicies[sample(nrow(indicies), 1), ]
  connectance <- as.numeric(rownames(meanData)[[indicies[[1]]]])
  species <- as.numeric(colnames(meanData)[[indicies[[2]]]])
  
  #selecting a random community from the communities that satisfy Nstar
  community <- a$communities[[sample(length(community), 1)]][[indicies[[1]]]][[indicies[[2]]]]
  
  #plotting a subset of the species in the Nstar community, integrated through time
  subset <- subsetPath(community, nI, connectance)
  subsetMelt <- melt(subset, id.vars = "time")
  colnames(subsetMelt) <- c("time", "Species", "Density")
  densityPlot <- ggplot(subsetMelt) +
    ylab("Species Density") +
    xlab("Time") +
    geom_line(mapping = aes(x=time, y=Density, color = Species))
  print(densityPlot)
  return(c(connectance, round(species*0.25)))
}