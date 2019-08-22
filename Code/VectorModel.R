#!/usr/bin/env Rscript

library(deSolve)
library(lattice)
library(tidyverse)
library(reshape)
library(plyr)
library(doParallel)

VectorCvNs <- function(S, C, step, p.m, p.e, s = 1, alpha = 0.5, xo = 10, r = 1, K = 100, h = 20, delta = 0.001, tl = 5000, e = 0.01, replicates = 10) {
  ## make some random, cascade, and niche food webs
  ## S is species richness
  ## C is connectance
  N <- 1     ## set the number of replicate webs to make
  matrixSize <- round(C/step)
  
  reps <- 1:replicates
  cSteps <- seq(step, C, step)
  x <- S/(C/step)
  sSteps <- seq(x, sum(rep(x, matrixSize)), x)
  
  #creating a list of lists of lists of dimensions N x C x S
  #C x S throws errors if it's not square
  
  numCores <- detectCores()
  #parCluster <- makeCluster(numCores, type="FORK")
  parCluster <- makeCluster(numCores, type="PSOCK")
  registerDoParallel(parCluster)
  
  out <- foreach(i = reps) %:%
    foreach(C_step = cSteps) %:%
    foreach(S_step = sSteps, .export = c("fullSim", "QianMatrix"), .packages = c("deSolve", "lattice")) %dopar% {
      fullSim(S_step, C_step, p.m, p.e, s, alpha, xo, r, K, h, delta, tl, e)
    }
  
  stopCluster(parCluster)
  
  communities <- foreach(i = out) %:%
    foreach(j = i) %:%
    foreach(k = j) %do% k[[1]]
  
  persistences <- foreach(i = out) %:%
    foreach(j = i) %:%
    foreach(k = j) %do% k[[2]]
  
  interactions <- foreach(i = out) %:%
    foreach(j = i) %:%
    foreach(k = j) %do% k[[3]]
  
  rm(out)
  
  # for(i in 1:length(interactions)){
  #   for(j in 1:length(interactions[[i]])) {
  #     for(k in 1:length(interactions[[i]][[j]])){
  #       interactions[[i]][[j]][[k]][[1]] <- as.big.matrix(interactions[[i]][[j]][[k]][[1]])
  #       interactions[[i]][[j]][[k]][[2]] <- as.big.matrix(interactions[[i]][[j]][[k]][[2]])
  #       interactions[[i]][[j]][[k]][[3]] <- as.big.matrix(interactions[[i]][[j]][[k]][[3]])
  #       interactions[[i]][[j]][[k]][[4]] <- as.big.matrix(interactions[[i]][[j]][[k]][[4]])
  #     }
  #   }
  # }
  
  #initializing containers to be returned
  means <- rep(list(rep(list(0), matrixSize)), matrixSize)
  stdDevs <- rep(list(rep(list(vector("integer", replicates)), matrixSize)), matrixSize)
  
  for(CvN in 1:replicates) {
    
    for(rowC in 1:matrixSize) {
      
      for(colN in 1:matrixSize) {
        
        if(!length(persistences[[CvN]][[rowC]])==0){
          means[[colN]][[rowC]] <- means[[colN]][[rowC]] + persistences[[CvN]][[rowC]][[colN]]
          stdDevs[[colN]][[rowC]][[CvN]] <- persistences[[CvN]][[rowC]][[colN]]
        }
      }
    }
  }
  
  for(rowC in 1:matrixSize) {
    for(colN in 1:matrixSize) {
      means[[colN]][[rowC]] <- means[[colN]][[rowC]]/replicates
      stdDevs[[colN]][[rowC]] <- sd(stdDevs[[colN]][[rowC]])
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
  
  #remove row a$mean[-c(2), ] where 2 is the row number and a$mean is the dataframe
  remove <- which(rowSums(means)==0)
  if(length(remove) > 0){
    means <- means[-c(remove), ]
    stdDevs <- stdDevs[-c(remove), ]
  }
  #returning containers
  matrices <- list("communities" = communities, "persistences" = persistences, "mean" = means, "stdDev" = stdDevs, "interactions" = interactions, "parms" = c(S, C, step, p.m, p.e, 1-(p.m+p.e)))
  return(matrices)
  
}

#graphs mean or sd as heatmap
heatMap <- function(dataFrame, graphType, xax="Species", yax="Connectance") {
  ggplot(melt(as.matrix(dataFrame))) +
  scale_fill_gradient(low = "steelblue", high = "white") +
  ylab(yax) +
  xlab(xax) +
  geom_tile(mapping = aes(x = X1, y = X2, fill = value)) +
  labs(fill = graphType)
}

#returns subset of community integrated through time using original final densities
VectorPath <- function(community, interactions, numSpecies, C, e, k) {
  
  livingAndDead <- community[nrow(community),2:ncol(community)] > e #getting the status of each species
  w <- list()
  counter = 1
  for(i in 1:length(livingAndDead)){
    if(livingAndDead[i]) {
      w[paste("Species", i+1)] <- i+1
      counter <- counter + 1
    }
  }
  
  #Exit if the sample size is greater than the number of persisting species in the community
  nStar <- length(w)
  if(numSpecies > nStar){
    print("nI is greater than N*")
    return(NA)
  }
  
  xo <- NULL
  y <- NULL
  BAL <- list()
    
  step <- 1
      
  path <- pathSim(w, numSpecies, community, interactions)
  BAL[step] <- list(path[[1]])
  while((length(path[[1]]$Living) < nStar) && (step < k)) {
    path <- pathSim(w, numSpecies, community, interactions, path[[2]])
    step <- step + 1
    BAL[step] <- list(path[[1]])
  }
  
  return(BAL)
}

nStarGraph <- function(container, Nstar, tolerance = 0.5, nI = 5,
                       graphStep = 1, replicates = 10, e = 0.01, k = 5000) {
  #creating matrix of NStar for all CvNs
  meanData <- container$mean
  communities <- container$communities
  interactions <- container$interactions
  nStarMatrix <- data.frame(matrix(nrow=nrow(meanData), ncol=ncol(meanData)))
  for(i in 1:ncol(nStarMatrix)){
    nStarMatrix[i] <- meanData[,i]*as.numeric(colnames(meanData)[i])
  }
  
  colnames(nStarMatrix) <- colnames(meanData)
  rownames(nStarMatrix) <- rownames(meanData)
  print(heatMap(nStarMatrix, "N*", xax="Species", yax="Connectance"))
  print("N* Matrix:")
  print(nStarMatrix)
  
  #finding the communities that have the desired Nstar within the set interval
  indicies <- which(nStarMatrix > Nstar-tolerance & nStarMatrix < Nstar+tolerance, arr.ind = TRUE)
  indicies <- as.data.frame(indicies)
  if(nrow(indicies)==0)
    return("No communities with that persistence.")
  colnames(indicies) <- c("Connectance", "Species")
  colnames(nStarMatrix) <- colnames(meanData)
  rownames(nStarMatrix) <- rownames(meanData)
  
  #Select random community from indicies of acceptable communities
  indicies <- indicies[sample(nrow(indicies), 1), ]
  connectance <- as.numeric(rownames(meanData)[[indicies[[1]]]])
  
  #selecting a random community from the communities that satisfy Nstar
  rpNum <- sample(length(communities), 1)
  community <- communities[[rpNum]][[indicies[[1]]]][[indicies[[2]]]]
  interaction <- interactions[[rpNum]][[indicies[[1]]]][[indicies[[2]]]]
  
  if(is.null(community))
    return("There are no communities that satisfy Nstar. Increase the search interval or pick a new value.")
  
  print(c(indicies[[1]], indicies[[2]]))

  #plotting a subset of the species in the Nstar community, integrated through time
  wSize <- sum(community[nrow(community),2:ncol(community)] > e)
  stepSize <- graphStep
  subset <- list()
  frames <- list()
  for(i in 1:replicates) {
    subset[i] <- list(VectorPath(community, interaction, nI, connectance, e, k))
    z <- c()
    for(j in 1:length(subset[[i]])){
      z <- c(z, lengths(subset[[i]][[j]]["Living"], use.names=FALSE))
    }
    z <- data.frame(z)
    colnames(z) <- i
    z["Step"] <- c(1:nrow(z))
    frames[[i]] <- z[seq(1, nrow(z), stepSize), ]
  }
  frames <- Reduce(function(x, y) merge(x=x, y=y, by="Step"), frames)
  frames["Mean"] <- rowMeans(frames[2:length(frames)])
  yield <- frames
  frames <- melt(frames, id.var="Step")
  colnames(frames) <- c("Step", "Replicates", "value")
  
  stepPlot <- ggplot(data=frames, aes(x=Step, y=value, col=Replicates)) +
    geom_line() +
    geom_point() +
    ggtitle(paste("Archipelago Migration Simulation: Qian Model")) +
    xlab("Step Number") +
    ylab("Nisle") +
    expand_limits(y=c(0,wSize))
  
  print(stepPlot)
  return(yield)
}