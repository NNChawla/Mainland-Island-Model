#!/usr/bin/env Rscript

VectorMetaMatrix <- function(container, nI = 5, graphStep = 1, replicates = 10, e = 0.01, kl = 5000) {
  
  mainS <- container$parms[[1]]
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
  
  cs <- as.numeric(rownames(meanData))
  ss <- as.numeric(colnames(meanData))

  numCores <- detectCores()
  parCluster <- makeCluster(numCores, type="PSOCK")
  registerDoParallel(parCluster)
  out <- foreach(i = cs) %:%
    foreach(j = ss, .export = c("VectorPath", "pathSim", "QianMatrix"), .packages = c("deSolve", "lattice")) %dopar% {
      rpNum <- sample(length(mainlands), 1)
      c <- which(cs == i)
      s <- which(ss == j)
      
      community <- mainlands[[rpNum]][[c]][[s]]
      interaction <- interactions[[rpNum]][[c]][[s]]
      if(is.null(community)) {
        list("NULL", "NULL", "NULL")
      }
      else {
        
        mainLiving <- sum(community[nrow(community), 2:ncol(community)] > e)
        if(mainLiving > nI) {
          stepSize <- graphStep
          subset <- list()
          frames <- list()
          pldRs <- list()
          steps <- list()
          
          for(k in 1:replicates) {
            
            pldRs[[k]] <- VectorPath(community, interaction, nI, i, e, kl)
            steps[[k]] <- sapply(pldRs[[k]], FUN = function(x) x$Steps)
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
          frames <- Reduce(function(x, y) merge(x=x, y=y, by="Step", all.y = TRUE, all.x = TRUE), frames)
          colnames(frames) <- c("Step", paste(i, j, 1:replicates))
          list(frames, pldRs, steps)
        }
        else {
          frames <- as.data.frame(rbind(c(1, rep(0, replicates))))
          colnames(frames) <- c("Step", paste(i, j, 1:replicates))
          pldRs <- list()
          for(k in 1:replicates){
            pldRs[[k]] <- NA
          }
          list(frames, pldRs, NA)
        }
      }
    }
  
  stopCluster(parCluster)
  
  rm(interactions)
  rm(mainlands)
  
  meanMat <- list()
  pathLDRatios <- list()
  eqLengths <- list()
  
  for(i in 1:length(out)){
    meanMat[[i]] <- list()
    pathLDRatios[[i]] <- list()
    eqLengths[[i]] <- list()
    for(j in 1:length(out[[i]])){
      meanMat[[i]][[j]] <- out[[i]][[j]][[1]]
      pathLDRatios[[i]][[j]] <- out[[i]][[j]][[2]]
      eqLengths[[i]][[j]] <- out[[i]][[j]][[3]]
    }
  }
  
  rm(out)
  
  frames <- list()
  for(i in 1:length(meanMat)){
    frames[[i]] <- Reduce(function(x, y) merge(x=x, y=y, by="Step", all.y = TRUE, all.x = TRUE), meanMat[[i]])
  }
  frames <- Reduce(function(x, y) merge(x=x, y=y, by="Step", all.y = TRUE, all.x = TRUE), frames)
  
  for(i in 1:length(pathLDRatios)){
    for(j in 1:length(pathLDRatios[[i]])){
      if(sum(is.na(pathLDRatios[[i]][[j]][[1]])) == 1) {
        pathLDRatios[[i]][[j]] <- NA
        next
      }
      for(k in 1:replicates){
        z <- c()
        for(l in 1:length(pathLDRatios[[i]][[j]][[k]])){
          live <- length(pathLDRatios[[i]][[j]][[k]][[l]][["Living"]])
          dead <- length(pathLDRatios[[i]][[j]][[k]][[l]][["Before"]])
          z <- c(z, live/dead)
        }
        pathLDRatios[[i]][[j]][[k]] <- z
      }
    }
  }
  
  return(list("Frames" = frames, "Parms" = container$parms, "TimeMatrix" = "TBA", "PLDRS" = pathLDRatios, "Steps" = eqLengths, "Data"="Data"))
}

timeMatrix <- function(massMat) {
  timeMat <- data.frame(matrix(nrow=3, ncol=ncol(massMat)-1))
  colnames(timeMat) <- colnames(massMat)[2:length(massMat)]
  rownames(timeMat) <- rownames(c("nFinalP", "nFinalS", "nHalfMS"))
  for(i in 2:length(massMat)){
    path <- massMat[[i]][!is.na(massMat[[i]])] #removing NAs from path so functions will work
    nFinalP <- path[[length(path)]] #Final Persistence Reached
    nFinalS <- length(path) #Step at which community reached mainland's final persistence
    nHalfM <- which.min(abs(path-0.5)) #Step at which community reched half of the mainland's final persistence
    timeMat[[i-1]] <- c(nFinalP, nFinalS, nHalfM)
  }
  return(timeMat)
}

tbTimeMatrix <- function(timeMat) { # for a single metaMatrix run, takes the timeMat and reformats it
  headers <- c("C", "N", "R", "p100", "t100", "t50")
  dtable <- as.tibble(matrix(0,0,length(headers)))
  for(x in 1:length(timeMat)){
    cnr <- unlist(strsplit(names(timeMat)[[x]], " "))
    CS <- as.double(cnr[[1]])
    NS <- as.integer(cnr[[2]])
    RS <- as.integer(cnr[[3]])
    NFP <- as.double(timeMat[[x]][[1]])
    NFS <- as.double(timeMat[[x]][[2]])
    THalf <- as.double(timeMat[[x]][[3]])
    
    dtable <- rbind(dtable, c(CS, NS, RS, NFP, NFS, THalf))
  }
  names(dtable) <- headers
  return(dtable)
}

tbFrames <- function(frames, replicates){
  
  plotMat <- frames[,2:ncol(frames)]
  mats <- list()
  count <- 1
  
  na.pad <- function(x,len){
    x[1:len]
  }
  
  makePaddedDataFrame <- function(l,...){
    maxlen <- max(sapply(l,length))
    frame <- data.frame(lapply(l,na.pad,len=maxlen),...)
    colnames(frame) <- names(l)
    return(frame)
  }
  
  for(i in seq(1, ncol(plotMat), replicates)){
    path <- plotMat[, i:(i+4)]
    if(all(is.na(path[1,ncol(path)]))) path[1,1:ncol(path)] <- rep(0, ncol(path))
    vals <- makePaddedDataFrame(lapply(path, FUN = function(x) rle(x[!is.na(x)])$values))
    rps <- makePaddedDataFrame(lapply(path, FUN = function(x) rle(x[!is.na(x)])$lengths))
    rps <- melt(rps, id.vars = NULL)[[2]]
    path <- melt(vals, id.vars = NULL) %>% as_tibble() %>% separate(variable, c("C", "N", "R"), sep = "([\\ ])")
    path <- path %>% add_column(RPS = rps) %>% drop_na
    path <- path %>% add_column(Step = unlist(lapply(rle(path$R)$lengths, FUN = function(y) seq(1:y))))
    mats[[count]] <- path
    count <- count + 1
  }
  mats <- bind_rows(mats)
  
  return(mats)
}

tbSteps <- function(eqLengths, tableFrame){
  
  mats <- list()
  count <- 1
  
  na.pad <- function(x,len){
    x[1:len]
  }
  
  paddedFrame <- function(l,...){
    maxlen <- max(sapply(l,length))
    frame <- data.frame(lapply(l,na.pad,len=maxlen),...)
    colnames(frame) <- seq(1, ncol(frame))
    return(frame)
  }
  
  for(i in 1:length(unique(tableFrame$C))){
    tmpFrame <- paddedFrame(sapply(eqLengths[[i]], paddedFrame))
    tmpFrame[1,] <- replace(tmpFrame[1,], is.na(tmpFrame[1,]), 0)
    colnames(tmpFrame) <- unlist(lapply(tableFrame$C[[i]], sapply(unique(tableFrame$N), 1:length(unique(tableFrame$R)), FUN = paste), FUN = paste))
    tmpFrame <- melt(tmpFrame, id.vars = NULL) %>% as_tibble() %>% separate(variable, c("C", "N", "R"), sep = "([\\ ])")
    tmpFrame <- tmpFrame %>% drop_na
    colnames(tmpFrame) <- c("C", "N", "R", "eqSteps")
    tmpFrame <- tmpFrame %>% add_column(Step = unlist(lapply(rle(tmpFrame$R)$lengths, FUN = function(y) seq(1:y))))
    mats[[count]] <- tmpFrame
    count <- count + 1
  }
  
  eqs <- bind_rows(mats)$eqSteps
  rps <- tableFrame$RPS
  mats <- tableFrame %>% add_column(eqSteps = seq(1, nrow(tableFrame)))
  
  neqs <- rep(NA, length(rps))
  i <- 1
  j <- 1
  
  while(i < length(rps)) {
    if(rps[[i]]!=1){
      neqs[[i]] <- mean(eqs[j:(j+rps[[i]]-1)])
      j <- j + rps[[i]]
      i <- i + 1
    }
    else{
      neqs[[i]] <- eqs[[j]]
      i <- i + 1
      j <- j + 1
    }
  }
  mats$eqSteps <- round(neqs)
  
  return(mats)
}