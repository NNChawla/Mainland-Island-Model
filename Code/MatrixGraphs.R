meanMatrix <- function(container, nI = 5, replace_sp = TRUE, graphStep = 1, replicates = 10, stepTime = 100) {
  #creating matrix of NStar for all CvNs
  modelType <- container$model
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
      
      mainLiving <- sum(community[nrow(community), 2:length(community)] > 10^-15)
      stepSize <- graphStep
      subset <- list()
      frames <- list()
      for(k in 1:replicates) {
        subset[k] <- list(subsetPath(community, interaction, nI, connectance, replace_sp, stepTime))
        z <- c()
        for(l in 1:length(subset[[k]])){
          z <- c(z, lengths(subset[[k]][[l]]["Living"], use.names=FALSE))
        }
        z <- z/mainLiving
        z <- data.frame(z)
        colnames(z) <- k
        z["Step"] <- c(1:nrow(z))
        if(is.na(subset[k]))
          z[[1]] <- NA
        frames[[k]] <- z[seq(1, nrow(z), stepSize), ]
      }
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
  return(list(frames, c(modelType, nI, stepTime)))
}

immStep <- c(2, 5, 10)
timeStep <- c(10, 50, 100, 200)
massMatrix <- function(containers, imms, times){
  plotMeans <- list()
  for(i in 1:length(containers)) {
    plotMeans[[i]] <- list()
    for(j in 1:length(imms)){
      plotMeans[[i]][[j]] <- list()
      for(k in 1:length(times)){
        plotMeans[[i]][[j]][[k]] <- meanMatrix(containers[[i]], nI=imms[[j]], stepTime=times[[k]])
      }
    }
  }
  return(plotMeans)
}

timeMatrix <- function(massMat) {
  timeMat <- data.frame(matrix(nrow=3, ncol=ncol(massMat)-1))
  colnames(timeMat) <- colnames(massMat)[2:length(massMat)]
  rownames(timeMat) <- rownames(c("nFinal", "nHalfI", "nHalfM"))
  for(i in 2:length(massMat)){
    nFinal <- massMat[[i]][[length(massMat[[i]])]]
    nHalfI <- which.min(abs(massMat[[i]]-nFinal/2))
    nHalfM <- which.min(abs(massMat[[i]]-0.5))
    timeMat[[i-1]] <- c(nFinal, nHalfI, nHalfM)
  }
  #timeMat[[length(timeMat)]] <- rowMeans(timeMat[1:(length(massMat)-1)])
  return(timeMat)
}

massGraph <- function(massMat, paths = c("Mean"), graphMean = FALSE, include = TRUE) {
  modelType <- massMat[[2]][[1]]
  nI <- massMat[[2]][[2]]
  sT <- massMat[[2]][[3]]
  massMat <- massMat[[1]]
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
  
  yLimit <- c(0, 1.0)

  stepPlot <- ggplot(data=mat, aes(x=Step, y=value, col=Replicates)) +
    geom_line() +
    geom_point() +
    ggtitle(paste("Archipelago Migration Simulation:", modelType, "Model - nI", nI, "- timeStep", sT)) +
    xlab("Step Number") +
    ylab("Nisle") +
    expand_limits(y=yLimit)
  
  return(stepPlot)
}

multiGraph <- function(multiMat){
  graphs <- list()
  count <- 1
  for(i in 1:length(multiMat)){
    for(j in 1:length(multiMat[[i]])){
      for(k in 1:length(multiMat[[i]][[j]])){
        graphs[[count]] <- massGraph(multiMat[[i]][[j]][[k]][[1]], multiMat[[i]][[j]][[k]][[2]], paths=c("All"))
        count <- count + 1
      }
    }
  }
  return(graphs)
}