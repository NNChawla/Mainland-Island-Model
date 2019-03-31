library(deSolve)
library(lattice)
library(tidyverse)
library(reshape)
library(plyr)

CvNs <- function(S, C, step, intTime = 100) {

  ## make some random, cascade, and niche food webs
  ## S is species richness
  ## C is connectance
  N <- 1     ## set the number of replicate webs to make
  matrixSize <- round(C/step)
  replicates <- 10
  
  #initializing containers to be returned
  communities <- rep(list(rep(list(vector("list", matrixSize)), matrixSize)), replicates)
  persistences <- rep(list(rep(list(vector("list", matrixSize)), matrixSize)), replicates)
  interactions <- rep(list(rep(list(vector("list", matrixSize)), matrixSize)), replicates)
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
        if (((S_step^2 - S_step)/2 - L) < -0.5) {
          S_step <- S_step + S/(C/step)
          next
        }
        
        #generating model based on parameters
        xxx <- Cascade.model(S_step, L, N)
        
        # Number of Species
        n <- S_step
        r <- runif(n, -1,1) #growth rate
        s <- runif(n, 1,1) #desnity dependent inhibition
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
        time <- list(start = 0, end = intTime, steps = intTime)
        # dummy variable for lvm() function defined above
        parms <- c(0) ### dummy variable (can have any numerical value)
        
        
        out <- n.integrate(time, init.x, model = mougi_model)
        communities[[CvN]][[rowC]][[colN]] <- out
        persistences[[CvN]][[rowC]][[colN]] <- mean(out[nrow(out),2:n+1] > 10^-15)
        interactions[[CvN]][[rowC]][[colN]] <- list(r, s, g, a, init.x)
        
        S_step <- S_step + S/(C/step)
      }
      
      S_step <- S/(C/step)
      C_step <- C_step + step
      
    }
  
  }
  print("Hello")
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
  print("I am here")
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
  print("here now")
  #returning containers
  matricies <- list("communities" = communities, "persistences" = persistences, "mean" = means, "stdDev" = stdDevs, "interactions" = interactions)
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
subsetPath <- function(community, numSpecies, C, replace_sp, stepTime = 100) {
  
  livingAndDead <- community[nrow(community),2:length(community)] > 10^-15 #getting the status of each species
  w <- list()
  counter = 1
  for(i in 1:length(livingAndDead)){
    if(livingAndDead[i]) {
      w[paste("Species", i+1)] <- i+1
      counter <- counter + 1
    }
  }
  
  #print("2")
  
  #Mandatory Check
  nStar <- length(w)
  if(numSpecies > nStar){
    print("nI is greater than N*")
    return(NA)
  }
  
  xo <- NULL
  y <- NULL
  BAL <- list()
  
  wSize <- 10000000 #arbitrarily large value
  if(!replace_sp) {
    numSteps <- ceiling(nStar/numSpecies)
  }
  else
    numSteps <- 100 ################## Eventually, we will do away with this and dynamically decide when to stop at runtime
  # We will do this by identifying when the island has reached equilibrium and cut it off there
  
  #print("3")
  
  for(step in 1:numSteps) {
    
    if(!replace_sp & (numSpecies > wSize))
      numSpecies <- wSize
    
    #print("4")
    
    if(step == 1) {
      
      tmpBAL <- list()
      
      #Sample W
      xo <- sample(w, numSpecies)
      if(!replace_sp)
        w <- setdiff(w, xo)
      
      #print("5")
      ############################################################################### Since this is a list it removes the names when you setdiff; needs to re-name each element as "Species Element"
      
      #Retrieving the appropriate persistences
      for(i in unlist(xo, use.names=FALSE)){
        xo[paste("Species", i)] <- community[nrow(community), i]
      }
      
      
      tmp <- list()
      for(i in 1:length(xo)){
        if(is.null(tmp[[names(xo[i])]]))
          tmp[[names(xo[i])]] <- xo[[i]]
        else
          tmp[[names(xo[i])]] <- tmp[[names(xo[i])]] + xo[[i]]
      }
      xo <- tmp
      
      tmpBAL["Before"] <- list(xo)
      
      #print("7")
      
      L <- round(numSpecies^2*C)
      N <- 1
      
      xxx <- Cascade.model(numSpecies, L, N)
      n <- numSpecies
      r <- runif(n, -1,1)
      s <- runif(n, 1,1)
      g <- runif(n)
      a <- xxx * matrix(runif(n*n, 0,1),nrow=n)
      diag(a) <- rep(0,n)
      
      init.x <- unlist(xo, use.names=FALSE)
      
      mougi_model <- function(t,x,parms){
        dx <- x * (r - s*x + g * (a %*% x) - (t(a) %*% x))
        list(dx)
      }
      
      n.integrate <- function(time=time, init.x= init.x, model=model){
        t.out <- seq(time$start,time$end,length=time$steps)
        as.data.frame(lsoda(init.x, t.out, model, parms = parms))
      }
      
      time <- list(start = 0, end = stepTime, steps = stepTime)
      parms <- c(0)
      tmp <- n.integrate(time, init.x, model = mougi_model)
      
      tmpBAL["Between"] <- list(tmp)
      #print("8")
      
      tmp <- unlist(tmp[nrow(tmp), 2:length(tmp)], use.names=FALSE)
      for(i in 1:length(tmp))
      {
        xo[i] <- tmp[i]
      }
      
      tmpBAL["After"] <- list(xo)
      
      #print("9")
      
      count <- length(xo)
      i <- 1
      while(!(i>count)) {
        if(xo[i] < 10^-15){
          xo[i] <- NULL
          count <- length(xo)
        }
        else{
          i <- i + 1
        }
      }
      
      tmpBAL["Living"] <- list(xo)
      
      #print("10")
      
      if(!replace_sp)
        wSize <- length(w)
      
      BAL[step] <- list(tmpBAL)
      
      #print("11")
    }
    
    else 
    {
      tmpBAL <- list()
      
      #Sample W
      y <- sample(w, numSpecies)
      if(!replace_sp)
        w <- setdiff(w, y)
      
      #print("5")
      ############################################################################### Since this is a list it removes the names when you setdiff; needs to re-name each element as "Species Element"
      
      #Retrieving the appropriate persistences
      for(i in unlist(y, use.names=FALSE))
        y[paste("Species", i)] <- community[nrow(community), i]
      
      #print("6")
      y <- c(y, xo)
      
      tmp <- list()
      for(i in 1:length(y)){
        if(is.null(tmp[[names(y[i])]]))
          tmp[[names(y[i])]] <- y[[i]]
        else
          tmp[[names(y[i])]] <- tmp[[names(y[i])]] + y[[i]]
      }
      y <- tmp
      
      tmpBAL["Before"] <- list(y)
      
      ySize <- length(y)
      
      #print("7")
      
      L <- round(ySize^2*C)
      N <- 1
      
      xxx <- Cascade.model(ySize, L, N)
      n <- ySize
      r <- runif(n, -1,1)
      s <- runif(n, 1,1)
      g <- runif(n)
      a <- xxx * matrix(runif(n*n, 0,1),nrow=n)
      diag(a) <- rep(0,n)
      
      init.x <- unlist(y, use.names=FALSE)
      
      mougi_model <- function(t,x,parms){
        dx <- x * (r - s*x + g * (a %*% x) - (t(a) %*% x))
        list(dx)
      }
      
      n.integrate <- function(time=time, init.x= init.x, model=model){
        t.out <- seq(time$start,time$end,length=time$steps)
        as.data.frame(lsoda(init.x, t.out, model, parms = parms))
      }
      
      time <- list(start = 0, end = stepTime, steps = stepTime)
      parms <- c(0)
      
      tmp <- n.integrate(time, init.x, model = mougi_model)
      
      tmpBAL["Between"] <- list(tmp)
      #print("8")
      
      tmp <- unlist(tmp[nrow(tmp), 2:length(tmp)], use.names=FALSE)
      for(i in 1:length(tmp))
      {
        y[i] <- tmp[i]
      }
      
      tmpBAL["After"] <- list(y)
      
      #print("9")
      
      count <- length(y)
      i <- 1
      while(!(i>count)) {
        if(y[i] < 10^-15){
          y[i] <- NULL
          count <- length(y)
        }
        else{
          i <- i + 1
        }
      }
      
      tmpBAL["Living"] <- list(y)
      
      #print("10")
      
      xo <- y
      
      if(!replace_sp)
        wSize <- length(w)
      
      BAL[step] <- list(tmpBAL)
      #print("11")
    }
  }
  
  return(BAL)
}

nStarGraph <- function(container, Nstar, tolerance = 0.5, nI = 5, replace_sp = TRUE, graphStep = 1, replicates = 10, stepTime = 100) {
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
  indicies <- which(nStarMatrix > Nstar-tolerance & nStarMatrix < Nstar+tolerance, arr.ind = TRUE)
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
  L <- round(nI^2*connectance)
  if(((nI^2 - nI)/2 - L) < -0.5)
    return("L value below 0")
  community <- container$communities[[sample(length(community), 1)]][[indicies[[1]]]][[indicies[[2]]]]
  if(is.null(community))
    return("There are no communities that satisfy Nstar. Increase the search interval or pick a new value.")
  
  print(c(indicies[[1]], indicies[[2]]))

  #plotting a subset of the species in the Nstar community, integrated through time
  wSize <- sum(community[nrow(community),2:length(community)] > 10^-15)
  stepSize <- graphStep
  subset <- list()
  frames <- list()
  for(i in 1:replicates) {
    subset[i] <- list(subsetPath(community, nI, connectance, replace_sp, stepTime))
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
    ggtitle("Archipelago Migration Simulation") +
    xlab("Step Number") +
    ylab("Nisle") +
    expand_limits(y=c(0,wSize))
  
  print(stepPlot)
  return(yield)
}