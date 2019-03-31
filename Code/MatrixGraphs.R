meanMatrix <- function(container, nI = 5, replace_sp = TRUE, graphStep = 1, replicates = 10, stepTime = 100, modelType = "Cascade") {
  #creating matrix of NStar for all CvNs
  meanData <- container$mean
  mainlands <- container$communities
  interactions <- container$interactions
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
    #if(i!=4)
      #next
    for(j in 1:ncol(nStarMatrix)) {
      connectance <- as.numeric(rownames(meanData)[[i]])
      startingSpecies <- as.numeric(colnames(meanData)[[j]])
      
      L <- round(nI^2*connectance)
      if(identical(modelType, "Cascade")) {
        if(((nI^2 - nI)/2 - L) < -0.5) {
          meanMatrix[[i]][[j]] <- NULL
          next
        }
      }
      rpNum <- sample(length(mainlands), 1)
      #print(rpNum)
      community <- mainlands[[rpNum]][[i]][[j]]
      interaction <- interactions[[rpNum]][[i]][[j]]
      if(is.null(community)) {
        meanMatrix[[i]][[j]] <- NULL
        next
      }
      print(c(connectance, startingSpecies))
      
      stepSize <- graphStep
      subset <- list()
      frames <- list()
      for(k in 1:replicates) {
        subset[k] <- list(subsetPath(community, interaction, nI, connectance, replace_sp, stepTime))
        z <- c()
        for(l in 1:length(subset[[k]])){
          z <- c(z, lengths(subset[[k]][[l]]["Living"], use.names=FALSE))
        }
        z <- data.frame(z)
        colnames(z) <- k
        z["Step"] <- c(1:nrow(z))
        if(is.na(subset[k]))
          z[[1]] <- NA
        frames[[k]] <- z[seq(1, nrow(z), stepSize), ]
      }
      #if(j==ncol(nStarMatrix))
        #return(frames)
      frames <- Reduce(function(x, y) merge(x=x, y=y, by="Step", all.y = TRUE), frames)
      frames <- rowMeans(frames[2:length(frames)], na.rm = TRUE)
      mat <- data.frame(matrix(nrow=length(frames), ncol=2))
      colnames(mat) <- c("Step", paste(i, j))
      mat["Step"] <- c(1:length(frames))
      mat[paste(i, j)] <- frames
      meanMatrix[[i]][[j]] <- mat
    }
  }
  
  frames <- list()
  for(i in 1:length(meanMatrix)){
    frames[[i]] <- Reduce(function(x, y) merge(x=x, y=y, by="Step", all.y = TRUE), meanMatrix[[i]])
  }
  frames <- Reduce(function(x, y) merge(x=x, y=y, by="Step", all.y = TRUE), frames)
  frames["Mean"] <- rowMeans(frames[2:length(frames)], na.rm = TRUE)
  return(frames)
}

immStep <- c(5, 10)
timeStep <- c(10, 50, 100, 200)
massMatrix <- function(container, imms, times){
  plotMeans <- list()
  for(i in 1:length(imms)){
    plotMeans[[i]] <- list()
    for(j in 1:length(times)){
      plotMeans[[i]][[j]] <- meanMatrix(container, nI=imms[[i]], stepTime=times[[j]])
    }
  }
  return(plotMeans)
  
  for(i in 1:length(plotMeans)){
    for(j in 1:length(plotMeans[[i]])) {
      plotMeans[[i]][[j]] <- plotMeans[[i]][[j]]["Mean"]
    }
  }
}

matrixGraph <- function(massMat, paths = c("Mean"), graphMean = FALSE, include = TRUE) {
  if(sum(is.element(paths, names(massMat)), na.rm=TRUE) == length(paths)) {
    mat <- massMat['Step']
    for(i in paths){
      mat[i] <- massMat[i]
    }
  }
  else if(identical(paths, "All")) {
    mat <- massMat
  }
  else{
    return(print("Enter valid values for the desired path"))
  }
  
  if(graphMean && include) {
    mat["Graph Mean"] <- rowMeans(mat[2:length(mat)])
  }
  else if(graphMean) {
    tmpMat <- mat
    mat <- tmpMat['Step']
    mat["Graph Mean"] <- rowMeans(tmpMat[2:length(tmpMat)])
  }
  mat <- melt(mat, id.var="Step")
  colnames(mat) <- c("Step", "Replicates", "value")

  stepPlot <- ggplot(data=mat, aes(x=Step, y=value, col=Replicates)) +
    geom_line() +
    geom_point() +
    ggtitle("Archipelago Migration Simulation") +
    xlab("Step Number") +
    ylab("Nisle") +
    expand_limits(y=c(0,max(massMat[nrow(massMat), 2:ncol(massMat)])+5))
  
  print(stepPlot)
}