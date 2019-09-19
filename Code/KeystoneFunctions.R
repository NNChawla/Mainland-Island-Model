library("qgraph")

trophicIndexes <- function(m, numEdges) {
  
  m <- t(m)
  size <- nrow(m)
  nodes <- 1:size
  k <- numEdges+1
  mus <- matrix(nrow=size, ncol=size)
  epsilons <- matrix(nrow=size, ncol=size)
  aijs <- matrix(nrow=size, ncol=size)
  
  aCalc <- function(row, col) {
    mu <- mus[row,col]
    numerator <- epsilons[row, col]/mu
    denominator <- sum(epsilons[row,]/mu, na.rm = TRUE)
    return(numerator/denominator)
  }
  
  for (i in nodes) {
    for (j in nodes) {
      
      if(i == j){
        mus[i,j] <- NA
        epsilons[i,j] <- NA
        next
      }
      
      interactionType <- c(m[i,j] < 0, m[j,i] < 0)
      epsilons[i,j] <- abs(m[i,j]) + abs(m[j,i])
      
      if(interactionType[[1]] & !interactionType[[2]]){
        mus[i,j] <- sum(abs(m[,i]), na.rm = TRUE)
        
      } else if (!interactionType[[1]] & interactionType[[2]]) {
        mus[i,j] <- sum(abs(m[i,]), na.rm = TRUE)
        
      } else if (!interactionType[[1]] & !interactionType[[2]]) {
        mus[i,j] <- sum(abs(m[i,]), na.rm = TRUE)
        
      } else if (interactionType[[1]] & interactionType[[2]]){
        mus[i,j] <- sum(abs(m[,i]), na.rm = TRUE)
        
      } else {
        print("Not True or False")
      }
      
    }
  }
  
  for (i in nodes) {
    for (j in nodes) {
      aijs[i, j] <- aCalc(i, j)
    }
  }
  
  tis <- array(0, dim=c(size,size, numEdges))
  path <- c()
  
  travel <- function(node, path) {
    path <- c(path, node)
    
    if(length(path)==k){
      pathSum <- 1
      for (pair in 2:k) {
        pathSum <- pathSum * aijs[path[pair-1], path[pair]]
      }
      tis[path[1],path[k],k-1] <<- tis[path[1],path[k],k-1] + pathSum
      path <- path[-length(path)]
      return(path)
    }
    else if (length(path) < k && length(path) > 1){
      len <- length(path)
      pathSum <- 1
      for (pair in 2:len) {
        pathSum <- pathSum * aijs[path[pair-1], path[pair]]
      }
      tis[path[1],path[len],len-1] <<- tis[path[1],path[len],len-1] + pathSum
    }
    
    links <- setdiff(which(m[node,] != 0), path)
    for(link in links){
      path <- travel(link,path)
    }
    path <- path[-length(path)]
    return(path)
  }
  
  tmp <- lapply(nodes, travel, path = path)
  
  tiIndexes <- rep(0, size)
  for (i in nodes) {
    for (e in 1:numEdges){
      tiIndexes[i] <- tiIndexes[i] + sum(tis[i,,e], na.rm = TRUE)
    }
  }
  tiIndexes <- tiIndexes/numEdges
  
  return(list(tiIndexes, tis))
}

impactsIndex <- function(tis) {
  rows <- dim(tis)[[1]]
  cols <- dim(tis)[[2]]
  edges <- dim(tis)[[3]]
  tiM <- matrix(nrow=rows, ncol=cols)
  for(i in 1:rows){
    for(j in 1:cols) {
      tiM[i,j] = sum(tis[i,j,])/edges
    }
  }
  
  Q <- tiM - t(tiM)
  I <- diag(rows)
  M <- solve(I-Q) - I
  IMA <- rep(0, rows)
  for(i in 1:rows){
    IMA[i] <- sum(abs(M[i,]))
  }
  
  return(IMA)
}

KeystoneIndexes <- function(m, walk_length) {
  
  g <- qgraph(t(abs(m)), directed = TRUE, DoNotPlot = TRUE)
  m_centrality <- centrality(g, weighted = TRUE)
  values <- trophicIndexes(m, walk_length)
  impacts <- impactsIndex(values[[2]])
  
  return(list(m_centrality, values[[1]], impacts))
}

KeystoneSimulate <- function(interactions, community, walk_length) {
  
  mCommunity <- community[nrow(community),2:ncol(community)]
  m <- Reduce('+', interactions[1:4])
  
  delta <- interactions[[10]]
  tl <- interactions[[11]]
  e <- interactions[[12]]
  
  extinctSpecies <- which(mCommunity[nrow(mCommunity),] < e)
  
  A.m <- interactions[[1]][-extinctSpecies, -extinctSpecies]
  A.ep <- interactions[[2]][-extinctSpecies, -extinctSpecies]
  A.p <- A.m + A.ep
  A.en <- interactions[[3]][-extinctSpecies, -extinctSpecies]
  A.c <- interactions[[4]][-extinctSpecies, -extinctSpecies]
  
  h <- interactions[[5]]
  K <- interactions[[6]][-extinctSpecies]
  r <- interactions[[7]][-extinctSpecies]
  s <- interactions[[8]][-extinctSpecies]
  
  
  m <- m[-extinctSpecies, -extinctSpecies]
  mCommunity <- mCommunity[-extinctSpecies]
  mCommunity <- unlist(mCommunity, use.names = FALSE)
  
  predictions <- KeystoneIndexes(m, walk_length)
  
  simulate <- function(i) {
    
    init.x <- mCommunity[-i]
    A.m <- A.m[-i, -i]
    A.ep <- A.ep[-i, -i]
    A.p <- A.m + A.ep
    A.en <- A.en[-i, -i]
    A.c <- A.c[-i, -i]
    K <- K[-i]
    r <- r[-i]
    s <- s[-i]
    sK <- s/K
    
    rootfn <- function(t,x,parms){
      dx <- unlist(qian_model(t,x,parms))
      sum(abs(dx) > delta) - 0
    }
    
    qian_model <- function(t,x,parms){
      x <- pmax(x, 0)
      dx <- x * (r + (sK*x) + ((A.p %*% (x/(h+x))) + (A.c %*% x) + (A.en %*% x)/(h+x)))
      list(dx)
    }
    
    n.integrate <- function(time=time, init.x= init.x, model=model){
      t.out <- seq(time$start,time$end)
      as.data.frame(lsodar(init.x, t.out, model, parms = parms, rootfunc = rootfn))
    }
    
    time <- list(start = 0, end = tl)
    parms <- c(0)
    out <- n.integrate(time, init.x, model = qian_model)
    
    return(out[nrow(out),])
  }
  
  numStillAlive <- c(1:nrow(m))
  
  for (i in 1:nrow(m)) {
    #print(i)
    run <- simulate(i)
    numStillAlive[[i]] <- sum(run[2:ncol(run)] > e)
  }
  
  #print(numStillAlive)
  return(list(predictions, numStillAlive))
}
