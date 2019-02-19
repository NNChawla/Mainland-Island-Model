subsetPath <- function(community, numSpecies, C, replace_sp) {
  
  #print("1")
  
  livingAndDead <- community[nrow(community),2:length(community)] > 10^-5 #getting the status of each species
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
  if(numSpecies > nStar)
    return(print("nI is greater than N*"))
  
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
      
      stepTime <- 100
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
      
      stepTime <- 100
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