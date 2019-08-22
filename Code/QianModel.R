#!/usr/bin/env Rscript

QianMatrix <- function(n, C, p.m, p.e, s, alpha, K) {
  ##number of species
  #n <- 10
  
  ## connectance
  #C <- 0.5
  
  ##proportion of each interaciton type
  #p.m <- 0.3
  #p.e <- 0.3
  p.c <- 1 - (p.m + p.e)
  
  ##self regulation
  #s <- 1
  s <- rep(s, n)
  
  ##iniate matrix, and set diagonal to self regulation
  a <- array(0, dim = c(n,n))
  diag(a) <- s
  
  ## population off diagonals with interactions with probability C
  connected <- rbinom(sum(upper.tri(a, diag = F)), 1, C)
  a[upper.tri(a, diag = F)] <- connected
  
  ## make symmetrical (copy upper diagonal to lower diagonal)
  ## this is the binary interaction matrix
  A <- a + t(a) - diag(diag(a))
  
  ## find matrix indices for interactions (row and column ids)
  pairs_index <- which(upper.tri(A) == 1, arr.ind = TRUE)
  interaction_index <- array(dim = c(sum(A[which(upper.tri(A))]), 2))
  count <- 1
  for (i in 1:nrow(pairs_index)){
    if (A[pairs_index[i,1], pairs_index[i,2]] == 1){
      interaction_index[count,] <- pairs_index[i,]
      count <- count+1
    }
  }
  
  
  ## draw multinomial for each interaction to determine interaction type
  interaction_type <- rmultinom(nrow(interaction_index), 1, c(p.m, p.e, p.c))
  
  #parameter for half normal distributions (see below)
  #alpha <- 0.5
  
  ## Loop through all interactions and fill appropriately
  
  for (i in 1:nrow(interaction_index)){
    
    ## Mutualism
    if(interaction_type[1,i] == 1){
      
      A[interaction_index[i,1], interaction_index[i,2]] <- abs(rnorm(1,0, alpha))
      A[interaction_index[i,2], interaction_index[i,1]] <- abs(rnorm(1,0, alpha))
      
    }
    
    ## Exploitation
    if(interaction_type[2,i] == 1){
      
      x <- sample(c(-1,1), 1) ## randomly assign one species as exploiter or victim
      
      A[interaction_index[i,1], interaction_index[i,2]] <- x * abs(rnorm(1,0, alpha))
      A[interaction_index[i,2], interaction_index[i,1]] <- -x * abs(rnorm(1,0, alpha))
      
    }
    
    ## Competition
    if(interaction_type[3,i] == 1){
      
      A[interaction_index[i,1], interaction_index[i,2]] <- -abs(rnorm(1,0, alpha/K))
      A[interaction_index[i,2], interaction_index[i,1]] <- -abs(rnorm(1,0, alpha/K))
      
    }  
  }
  
  A.m <- matrix(rep(0, n*n),nrow=n)
  A.ep <- matrix(rep(0, n*n),nrow=n)
  A.en <- matrix(rep(0, n*n),nrow=n)
  A.e <- matrix(rep(0, n*n),nrow=n)
  A.c <- matrix(rep(0, n*n),nrow=n)
  
  for(i in 1:nrow(interaction_index)){
    
    indicies <- interaction_index[i,]
    int_type <- interaction_type[,i]
    
    if(interaction_type[1,i]) {
      
      A.m[indicies[1], indicies[2]] <- A[indicies[1], indicies[2]]
      A.m[indicies[2], indicies[1]] <- A[indicies[2], indicies[1]]
      
    } else if (interaction_type[2,i]) {
      
      A.e[indicies[1], indicies[2]] <- A[indicies[1], indicies[2]]
      A.e[indicies[2], indicies[1]] <- A[indicies[2], indicies[1]]
      
      if(A[indicies[1], indicies[2]] > 0) {
        A.ep[indicies[1], indicies[2]] <- A[indicies[1], indicies[2]]
        A.en[indicies[2], indicies[1]] <- A[indicies[2], indicies[1]]
      }
      else {
        A.ep[indicies[2], indicies[1]] <- A[indicies[2], indicies[1]]
        A.en[indicies[1], indicies[2]] <- A[indicies[1], indicies[2]]
      }
      
    } else {
      
      A.c[indicies[1], indicies[2]] <- A[indicies[1], indicies[2]]
      A.c[indicies[2], indicies[1]] <- A[indicies[2], indicies[1]]
      
    }
    
  }
  
  return(list(A.m, A.ep, A.en, A.c, A.e, interaction_index, interaction_type))
}

fullSim <- function(n, C, p.m, p.e, s, alpha, xo, r, K, h, delta, tl, e) {
  
  s <- -s
  
  interaction_matrices <- QianMatrix(n, C, p.m, p.e, s, alpha, K)
  A.m <- interaction_matrices[[1]]
  A.ep <- interaction_matrices[[2]]
  A.p <- A.m + A.ep
  A.en <- interaction_matrices[[3]]
  A.c <- interaction_matrices[[4]]
  interaction_index <- interaction_matrices[[6]]
  interaction_type <- interaction_matrices[[7]]
  
  s <- rep(s, n)
  r <- rep(r, n)
  K <- rep(K, n)
  sK <- s/K
  
  init.x <- runif(n) * xo
  
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
  
  #This is to preserve memory; remove if no longer an issue
  zout <- out[nrow(out),]
  
  steps <- nrow(out)
  return(list(zout, mean(out[nrow(out),2:ncol(out)] > e), list(A.m, A.ep, A.en, A.c, h, K, r, s, steps, delta, tl, e, interaction_index, interaction_type)))
}

pathSim <- function(w, numSpecies, community, interactions, y = NULL, xz = 0.1) {
  
  tmpBAL <- list()
  
  #Sample W
  xo <- sample(w, numSpecies)
  
  #Retrieving the appropriate persistences
  for(i in unlist(xo, use.names=FALSE)){
    xo[paste("Species", i)] <- community[nrow(community), i] * xz
  }
  
  xo <- c(xo, y)

  tmp <- list()
  for(i in 1:length(xo)){
    if(is.null(tmp[[names(xo[i])]]))
      tmp[[names(xo[i])]] <- xo[[i]]
    else
      tmp[[names(xo[i])]] <- tmp[[names(xo[i])]] + xo[[i]]
  }
  xo <- tmp

  species <- c()
  for(i in names(xo)){
    species <- c(species, as.integer(unlist(strsplit(i, " "))[[2]])-1)
  }

  tmpBAL["Before"] <- list(xo)
  
  A.m <- interactions[[1]][species, species]
  A.ep <- interactions[[2]][species, species]
  A.p <- A.m + A.ep
  A.en <- interactions[[3]][species, species]
  A.c <- interactions[[4]][species, species]
  
  h <- interactions[[5]]
  K <- interactions[[6]][species]
  r <- interactions[[7]][species]
  s <- interactions[[8]][species]
  sK <- s/K
  init.x <- unlist(xo, use.names = FALSE)
  
  delta <- interactions[[10]]
  tl <- interactions[[11]]
  e <- interactions[[12]]
  
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
  steps <- nrow(out)
  
  out <- unlist(out[nrow(out), 2:ncol(out)], use.names=FALSE)
  for(i in 1:length(out))
  {
    xo[i] <- out[i]
  }
  
  tmpBAL["After"] <- list(xo)
  
  count <- length(xo)
  i <- 1
  while(!(i>count)) {
    if(xo[i] < e){
      xo[i] <- NULL
      count <- length(xo)
    }
    else{
      i <- i + 1
    }
  }
  
  tmpBAL["Living"] <- list(xo)
  tmpBAL["Steps"] <- steps
  
  return(list(tmpBAL, xo))
}