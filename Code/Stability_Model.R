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
  
  #remove row a$mean[-c(2), ] where 2 is the row number and a$mean is the dataframe
  remove <- which(as.numeric(rownames(means)) > (0.5*(1-1/S)))
  if(length(remove) > 0){
    means <- means[-c(remove), ]
    stdDevs <- stdDevs[-c(remove), ]
  }
  
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
  
  #Getting the species that survived Mainland integration
  speciesSurvived <- community[nrow(community),2:length(community)] > 10^-5
  persistingSpecies <- c(1:sum(speciesSurvived == TRUE))
  counter = 1
  for(i in 1:length(speciesSurvived)){
    if(speciesSurvived[i]) {
      persistingSpecies[counter] <- i+1
      counter <- counter + 1
    }
  }
  nStar <- length(persistingSpecies)
  if(numSpecies > nStar)
    return(print("nI is greater than N*"))
  
  #the returned data set should be all of the integrations throughout time
  numSteps <- ceiling(nStar/numSpecies)
  allIslandPersistences <- c()
  allIslandIndicies <- c()
  endIslandIndicies <- c()
  endIslandPersistences <- c()
  subsets <- list()
  
  for(step in 1:numSteps) {
    
    if(numSpecies > remainingMainlandSpecies)
      numSpecies <- remainingMainlandSpecies
    #selecting random subset of nI species from the persisting species
    randSpecies <- sample(persistingSpecies, numSpecies, replace=FALSE)
    persistingSpecies <- setdiff(persistingSpecies, randSpecies)
    speciesIndicies <- c(1:numSpecies+1)
    speciesIndicies[1] <- "time"
    for(i in 1:numSpecies) {
      speciesIndicies[i+1] <- randSpecies[i]-1
      randSpecies[i] <- community[1, randSpecies[i]]
    }
    
    allIslandIndicies <- c(allIslandIndicies, speciesIndicies[2:numSpecies+1])
    preIntegration <- c(endIslandPersistences, randSpecies)
    
    #integration function
    L <- round(numSpecies^2*C)  ## calculate number of links from S and C
  
    xxx <- Cascade.model(numSpecies, L, N)
    
    n <- numSpecies
    r <- runif(n, -1,1)
    s <- runif(n, 1,1)
    g <- runif(n)
    a <- xxx * matrix(runif(n*n, 0,1),nrow=n)
    diag(a) <- rep(0,n)
    
    #what should this be?
    init.x <- preIntegration
    
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
    
    dataset <- n.integrate(time, init.x, model = mougi_model)
    subsets[step] <- dataset
    
    postIntegration <- unlist(dataset[nrow(dataset), 2:length(dataset)], use.names=FALSE)
    allIslandPersistences <- c(allIslandPersistences, postIntegration)
    
    speciesSurvived <- postIntegration > 10^-5
    survivors <- c(1:sum(speciesSurvived == TRUE))
    #counter = 1
    #for(i in 1:length(speciesSurvived)){
    #  if(speciesSurvived[i]) {
    #    survivors[counter] <- i+1
    #    counter <- counter + 1
    #  }
    #}
    remainingMainlandSpecies <- length(survivors)
    
    endIslandIndicies <- c(endIslandIndicies, survivors)
    for(i in endIslandIndicies){
      tmp <- match(i, postIntegration)
      endIslandPersistences <- c(endIslandPersistences, postIntegration[tmp])
    }
  }
  
  colnames(subsetData) <- speciesIndicies
  return(subsetData)
}

nStarGraph <- function(container, Nstar, interval=0.5, nI = 5) {
  #creating matrix of NStar for all CvNs
  meanData <- container$mean
  community <- container$communities
  nStarMatrix <- data.frame(matrix(nrow=nrow(meanData), ncol=ncol(meanData)))
  for(i in 1:ncol(nStarMatrix)){
    nStarMatrix[i] <- meanData[,i]*as.numeric(colnames(meanData)[i])
  }
  
  colnames(nStarMatrix) <- colnames(meanData)
  rownames(nStarMatrix) <- rownames(meanData)
  print(plotGraph(nStarMatrix, "N*", xax="Species", yax="Connectance"))
  print("N* Matrix:")
  print(nStarMatrix)
  
  #finding the communities that have the desired Nstar within the set interval
  indicies <- which(nStarMatrix > Nstar-interval & nStarMatrix < Nstar+interval, arr.ind = TRUE)
  indicies <- as.data.frame(indicies)
  if(nrow(indicies)==0)
    return("No communities with that persistence.")
  colnames(indicies) <- c("Connectance", "Species")
  colnames(nStarMatrix) <- colnames(meanData)
  rownames(nStarMatrix) <- rownames(meanData)
  
  #Deprecated check for invalid CvN Pairs
  # for(i in 1:nrow(indicies)) {
  #   R <- indicies[i,][[1]]
  #   C <- indicies[i,][[2]]
  #   S <- as.numeric(rownames(meanData)[R])
  #   C <- as.numeric(colnames(meanData)[C])
  #   if(round(S^2*C)==-1)
  #     print(c(S,C))
  # }
  
  #Select random community from indicies of acceptable communities
  indicies <- indicies[sample(nrow(indicies), 1), ]
  connectance <- as.numeric(rownames(meanData)[[indicies[[1]]]])
  startingSpecies <- as.numeric(colnames(meanData)[[indicies[[2]]]])

  #selecting a random community from the communities that satisfy Nstar
  L <- round(startingSpecies^2*connectance)
  if(((startingSpecies^2 - startingSpecies)/2 - L) == -1)
    return("L value below 0")
  community <- a$communities[[sample(length(community), 1)]][[indicies[[1]]]][[indicies[[2]]]]
  if(is.null(community))
    return("There are no communities that satisfy Nstar. Increase the search interval or pick a new value.")

  #plotting a subset of the species in the Nstar community, integrated through time
  subset <- subsetPath(community, nI, connectance)
  subsetMelt <- melt(subset, id.vars = "time")
  colnames(subsetMelt) <- c("time", "Species", "Density")
  densityPlot <- ggplot(subsetMelt) +
    ylab("Species Density") +
    xlab("Time") +
    geom_line(mapping = aes(x=time, y=Density, color = Species))
  return(densityPlot)
}