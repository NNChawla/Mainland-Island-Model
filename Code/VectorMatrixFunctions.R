VectorMetaMatrix <- function(container, nI = 5, graphStep = 1, replicates = 10, e = 0.01, kl = 1000) {
  tic()
  
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
        list("NULL", "NULL")
      }
      else {
        
        mainLiving <- sum(community[nrow(community), 2:length(community)] > e)
        stepSize <- graphStep
        subset <- list()
        frames <- list()
        pldRs <- list()
        
        for(k in 1:replicates) {
          
          pldRs[[k]] <- VectorPath(community, interaction, nI, i, e, kl)
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
  
  # for(i in 1:length(pathLDRatios)){
  #   for(j in 1:length(pathLDRatios[[i]])){
  #     if(sum(is.na(pathLDRatios[[i]][[j]][[1]])) == 1) {
  #       pathLDRatios[[i]][[j]] <- NA
  #       next
  #     }
  #     for(k in 1:replicates){
  #       z <- c()
  #       for(l in 1:length(pathLDRatios[[i]][[j]][[k]])){
  #         live <- length(pathLDRatios[[i]][[j]][[k]][[l]][["Living"]])
  #         dead <- length(pathLDRatios[[i]][[j]][[k]][[l]][["Before"]])
  #         z <- c(z, live/dead)
  #       }
  #       pathLDRatios[[i]][[j]][[k]] <- z
  #     }
  #   }
  # }
  
  toc()
  return(list(frames, pathLDRatios, c(nI, mainS)))
  #return(list(frames, c(nI, mainS), timeMatrix(frames), pathLDRatios))
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