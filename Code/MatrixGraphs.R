matrixGraph <- function(container, nI = 5, replace_sp = TRUE, graphStep = 1, replicates = 10, stepTime = 100) {
  #creating matrix of NStar for all CvNs
  meanData <- container$mean
  mainlands <- container$communities
  nStarMatrix <- data.frame(matrix(nrow=nrow(meanData), ncol=ncol(meanData)))
  for(i in 1:ncol(nStarMatrix)){
    nStarMatrix[i] <- meanData[,i]*as.numeric(colnames(meanData)[i])
  }
  
  colnames(nStarMatrix) <- colnames(meanData)
  rownames(nStarMatrix) <- rownames(meanData)
  meanMatrix <- list()
  
  #Data.Frame [r, c] List [[r]][[c]]
  print("N* Matrix:")
  print(nStarMatrix)
  
  for(i in 1:nrow(nStarMatrix)){
    meanMatrix[[i]] <- list()
    for(j in 1:ncol(nStarMatrix)) {
      
      connectance <- as.numeric(rownames(meanData)[[i]])
      startingSpecies <- as.numeric(colnames(meanData)[[j]])
      
      #selecting a random community from the communities that satisfy Nstar
      L <- round(nI^2*connectance)
      if(((nI^2 - nI)/2 - L) < -0.5) {
        meanMatrix[[i]][[j]] <- NA
        next
      }
      
      rpNum <- sample(length(mainlands), 1)
      #print(rpNum)
      community <- container$communities[[rpNum]][[i]][[j]]
      if(is.null(community)) {
        meanMatrix[[i]][[j]] <- NA
        next
      }
      print(c(connectance, startingSpecies))
      
      #plotting a subset of the species in the Nstar community, integrated through time
      stepSize <- graphStep
      subset <- list()
      frames <- list()
      for(k in 1:replicates) {
        subset[k] <- list(subsetPath(community, nI, connectance, replace_sp, stepTime))
        
        z <- c()
        for(l in 1:length(subset[[k]])){
          z <- c(z, lengths(subset[[k]][[l]]["Living"], use.names=FALSE))
        }
        z <- data.frame(z)
        colnames(z) <- k
        z["Step"] <- c(1:nrow(z))
        frames[[k]] <- z[seq(1, nrow(z), stepSize), ]
      }
      frames <- Reduce(function(x, y) merge(x=x, y=y, by="Step"), frames)
      
      frames <- rowMeans(frames[2:length(frames)])
      mat <- data.frame(matrix(nrow=length(frames), ncol=2))
      colnames(mat) <- c("Step", paste(i, j))
      mat["Step"] <- c(1:length(frames))
      mat[paste(i, j)] <- frames
      #print(c(i, j))
      meanMatrix[[i]][[j]] <- list(mat)
     }
  }
  return(meanMatrix)
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