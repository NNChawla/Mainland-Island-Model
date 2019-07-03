#original return values
#list(frames, c(modelType, nI, stepTime, mainS), timeMatrix(frames), pathLDRatios)
#new return values
#list(frames, c(nI, mainS), timeMatrix(frames), pathLDRatios)

library(zoo)

pltTbFrames <- function(mats, single = FALSE, connectance = 0.1, richness = 10){
  if(!single) {
    plt <- ggplot(mats, aes(x=Step, y=value, color = factor(R)))
    plt <- plt + facet_wrap(vars(factor(C)))
    plt <- plt + geom_point() + geom_line() + theme(legend.position = "None") 
  }
  else{
    mats <- mats %>% filter(C == connectance & N == richness)
    plt <- ggplot(mats, aes(x=Step, y=value, color=factor(R)))
    plt <- plt + geom_point() + geom_line() + xlab("Step Number") + ylab("Nisle")
    plt <- plt + ggtitle(paste("Connectance", connectance, "Species #", richness))
  }
  return(plt)
}

vPathGraph <- function(path, nStar) {
  plt <- list()
  
  for(i in 1:length(path)) {
    plt[[i]] <- length(path[[i]][["Living"]])/nStar
  }
  
  plt <- melt(plt)
  colnames(plt) <- c("value", "Step")
  
  yLimit <- c(0, 1.0)
  
  stepPlot <- ggplot(data=plt, aes(x=Step, y=value)) +
    geom_line() +
    geom_point() +
    xlab("Step Number") +
    ylab("Nisle") +
    expand_limits(y=yLimit)
  
  return(stepPlot)
}

massGraph <- function(massMat, paths = c("All"), graphMean = FALSE, include = TRUE) {
  nI <- massMat[[2]][[1]]
  massMat <- massMat[[1]]
  
  massMat <- na.locf(massMat)
  
  for(i in paths)
    mat <- matrix(nrow=nrow(massMat), ncol=ncol(massMat))
  mat[i] <- unlist(select(massMat, ends_with(toString(i)))[1:nrow(massMat),])
  browser()
  mat <- melt(mat, id.var="Step")
  colnames(mat) <- c("Step", "Replicates", "value")
  
  yLimit <- c(0, 1.0)
  
  stepPlot <- ggplot(data=mat, aes(x=Step, y=value, col=Replicates)) +
    geom_line() +
    geom_point() +
    ggtitle(paste("Archipelago Migration Simulation: Qian Model - nI", nI)) +
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

timeMatrixGraph <- function(massMat, rNames, cNames, t100){
  if(t100) {
    ylabel <- "t100"
    rN <- 2
  }
  else {
    ylabel <- "p100"
    rN <- 1
  }
  
  timeMat <- massMat[[3]]
  x <- 100#length(massMat[[4]])
  y <- max(lengths(massMat[[4]]))
  sC <- max(lengths(massMat[[4]][[1]]))
  
  vectors <- data.frame(matrix(nrow=x, ncol=y))
  vectors[, length(vectors)+1] <- as.double(rNames)[1:x]
  for(i in 1:y){
    vectors[, i] <- unlist(select(timeMat, ends_with(toString(i)))[rN,])
  }
  browser()
  rownames(vectors) <- rNames[1:x]
  colnames(vectors) <- c(cNames[1:y], "C")
  vectors <- melt(vectors, id.var="C")
  colnames(vectors) <- c("C", "N", "factor")
  yLimit <- c(0, sC)
  stepPlot <- ggplot(data=vectors, aes(x=C, y=factor, col=N)) +
    geom_line() +
    geom_point() +
    ggtitle("Archipelago t50s") +
    xlab("C") +
    ylab(ylabel)
  #expand_limits(y=yLimit)
  return(stepPlot)
}