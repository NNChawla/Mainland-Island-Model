VectorMetaMatrix <- function(container, nI = 5, replace_sp = TRUE, graphStep = 1, replicates = 10, stepTime = 100, stepCount = 100) {
  mainS <- container$parms[[1]]
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
  
  rownums <- 1:nrow(nStarMatrix)
  colnums <- 1:ncol(nStarMatrix)
  
  cond <- function(c, s){
    if(identical(modelType, "Cascade")) {
      # skipping Connectance and Species pairing because number of links is less than 0 which throws error
      return(ifelse((((s^2 - s)/2 - round(s^2*c)) < -0.5), FALSE, TRUE))
    }
    else if(identical(modelType, "Niche")) {
      return(TRUE)
    }
  }
  
  cs <- as.numeric(rownames(meanData))[cond(as.numeric(rownames(meanData)), nI)]
  ss <- as.numeric(colnames(meanData))
  
  numCores <- detectCores()
  parCluster <- makeCluster(numCores, type="PSOCK")
  registerDoParallel(parCluster)
  
  out <- foreach(i = cs) %:%
    foreach(j = ss, .export = c("subsetPath"), .packages = c("deSolve", "lattice")) %dopar% {
      rpNum <- sample(length(mainlands), 1)
      c <- which(cs == i)
      s <- which(ss == j)
      
      community <- mainlands[[rpNum]][[c]][[s]]
      interaction <- interactions[[rpNum]][[c]][[s]]
      if(is.null(community)) {
        list("NULL", "NULL")
      }
      else {
        
        mainLiving <- sum(community[nrow(community), 2:length(community)] > 10^-15)
        stepSize <- graphStep
        subset <- list()
        frames <- list()
        pldRs <- list()
        for(k in 1:replicates) {
          pldRs[[k]] <- subsetPath(community, interaction, nI, i, replace_sp, stepTime, stepCount)
          subset[k] <- list(pldRs[[k]])
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
        colnames(frames) <- c("Step", paste(i, j, 1:replicates))
        list(frames, pldRs)
      }
    }
  
  stopCluster(parCluster)
  
  meanMat <- list()
  pathLDRatios <- list()
  
  for(i in 1:length(out)){
    meanMat[[i]] <- list()
    pathLDRatios[[i]] <- list()
    for(j in 1:length(out[[i]])){
      meanMat[[i]][[j]] <- out[[i]][[j]][[1]]
      pathLDRatios[[i]][[j]] <- out[[i]][[j]][[2]]
    }
  }
  
  frames <- list()
  for(i in 1:length(meanMat)){
    frames[[i]] <- Reduce(function(x, y) merge(x=x, y=y, by="Step", all.y = TRUE), meanMat[[i]])
  }
  frames <- Reduce(function(x, y) merge(x=x, y=y, by="Step", all.y = TRUE), frames)
  
  for(i in 1:length(pathLDRatios)){
    for(j in 1:length(pathLDRatios[[i]])){
      if(sum(is.na(pathLDRatios[[i]][[j]][[1]])) == 1) {
        pathLDRatios[[i]][[j]] <- NA
        next
      }
      for(k in 1:replicates){
        z <- c()
        for(l in 1:stepCount){
          live <- length(pathLDRatios[[i]][[j]][[k]][[l]][["Living"]])
          dead <- length(pathLDRatios[[i]][[j]][[k]][[l]][["Before"]])
          z <- c(z, live/dead)
        }
        pathLDRatios[[i]][[j]][[k]] <- z
      }
    }
  }
  
  return(list(frames, c(modelType, nI, stepTime, mainS), timeMatrix(frames), pathLDRatios))
}

multiMatrix <- function(containers, imms, times){
  plotMeans <- list()
  for(i in 1:length(containers)) {
    plotMeans[[i]] <- list()
    for(j in 1:length(imms)){
      plotMeans[[i]][[j]] <- list()
      for(k in 1:length(times)){
        plotMeans[[i]][[j]][[k]] <- tryCatch(VectorMetaMatrix(containers[[i]], nI=imms[[j]], stepTime=times[[k]]), error=function(e) NA)
        print(c(i, j, k))
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

multiGraph <- function(multiMat, containers){
  massGraphs <- list()
  tGraphs <- list()
  pGraphs <- list()
  count <- 1
  for(i in 1:length(multiMat)){
    for(j in 1:length(multiMat[[i]])){
      for(k in 1:length(multiMat[[i]][[j]])){
        if(is.na(multiMat[[i]][[j]][[k]])){
          next
        }
        massMat <- multiMat[[i]][[j]][[k]]
        meanData <- containers[[i]]$mean
        S <- containers[[i]]$parms[[1]]
        fName <- paste(S, "-", massMat[[2]][[1]], "-", massMat[[2]][[2]], "-", massMat[[2]][[3]], ".png", sep="")
        massGraphs[[count]] <- massGraph(massMat, paths=c("All"))
        tGraphs[[count]] <- t50graph(massMat, rownames(meanData), colnames(meanData))
        pGraphs[[count]] <- pathGraph(massMat)
        ggsave(filename=paste("mass-", fName, sep=""), plot=massGraphs[[count]], dpi=320, width=20, height=10)
        ggsave(filename=paste("t50-", fName, sep=""), plot=tGraphs[[count]], dpi=320, width=20, height=10)
        ggsave(filename=paste("path-", fName, sep=""), plot=pGraphs[[count]], dpi=320, width=20, height=10)
        
        count <- count + 1
      }
    }
  }
  return(list(massGraphs, tGraphs, pGraphs))
}

matGraphs <- function(massMat, container) {
  meanData <- container$mean
  fName <- paste(massMat[[2]][[1]], "-", massMat[[2]][[2]], "-", massMat[[2]][[3]], ".png", sep="")
  mGraph <- massGraph(massMat, paths=c("All"))
  tGraph <- t50graph(massMat, rownames(meanData), colnames(meanData))
  pGraph <- pathGraph(massMat)
  ggsave(filename=paste("mass-", fName, sep=""), plot=mGraph, dpi=320, width=20, height=10)
  ggsave(filename=paste("t50-", fName, sep=""), plot=tGraph, dpi=320, width=20, height=10)
  ggsave(filename=paste("path-", fName, sep=""), plot=pGraph, dpi=320, width=20, height=10)
  return(list(mGraph, tGraph, pGraph))
}

pathGraph <- function(massMat, paths=c("All")){
  pathMat <- massMat[[4]]
  C <- rep(seq_along(pathMat), lengths(pathMat))
  N <- sequence(lengths(pathMat))
  step <- lengths(unlist(pathMat, rec=FALSE))
  graph <- data.frame(
    CvN = paste(rep(C, step), rep(N, step)),
    step = sequence(step),
    nIM = unlist(pathMat)
  )
  #return(graph)
  if(!identical(paths, "All"))
    graph <- filter(graph, CvN %in% paths)
  pathPlot <- ggplot(data=graph, aes(x=step, y=nIM, col=CvN)) +
  geom_line() +
  geom_point()
  
  return(pathPlot)
}
#pathGraph(pathLDRatios, paste(1, 1:9))

t50graph <- function(massMat, rNames, cNames){
  timeMat <- massMat[[3]]
  dims <- c(length(massMat[[4]]), max(lengths(massMat[[4]])))
  sC <- max(lengths(massMat[[4]][[1]]))
  
  x <- dims[[1]]
  y <- dims[[2]]
  vectors <- data.frame(matrix(nrow=x, ncol=y))
  vectors[, length(vectors)+1] <- as.double(rNames)[1:x]
  for(i in 1:y){
    vectors[, i] <- unlist(select(timeMat, ends_with(toString(i)))[2,])
  }
  rownames(vectors) <- rNames[1:x]
  colnames(vectors) <- c(cNames[1:y], "C")
  vectors <- melt(vectors, id.var="C")
  colnames(vectors) <- c("C", "N", "t50")
  yLimit <- c(0, sC)
  stepPlot <- ggplot(data=vectors, aes(x=C, y=t50, col=N)) +
    geom_line() +
    geom_point() +
    ggtitle("Archipelago t50s") +
    xlab("C") +
    ylab("t50")
    #expand_limits(y=yLimit)
  return(stepPlot)
}

dataTable <- function(ds, S, nI, sT) {
  headers <- c("mainS","nI", "nST", "C", "N", "Rep", "t50I")
  dtable <- as.tibble(matrix(0,0,length(headers)))
  
  MS = NIS = NSTS = CS = NS = RS = TS = 0
  for(i in 1:length(S)) {
    for(j in 1:length(nI)){
      for(k in 1:length(sT)){
        if(length(ds[[i]][[j]][[k]])==1) #if entry is NA
          next
        entry <- ds[[i]][[j]][[k]]
        MS <- as.integer(entry[[2]][[4]])
        NIS <- as.integer(entry[[2]][[2]])
        NSTS <- as.integer(entry[[2]][[3]])
        
        for(x in 1:length(entry[[3]])){
          cnr <- unlist(strsplit(names(entry[[3]])[[x]], " "))
          CS <- as.double(cnr[[1]])
          NS <- as.integer(cnr[[2]])
          RS <- as.integer(cnr[[3]])
          TS <- as.double(entry[[3]][[x]][[2]])
          
          dtable <- rbind(dtable, c(MS, NIS, NSTS, CS, NS, RS, TS))
        }
        
      }
    }
  }
  names(dtable) <- headers
  return(dtable)
}