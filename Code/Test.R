VectorCvNs <- function(S, C, step, p.m, p.e, s = 1, alpha = 0.5, xo = 10, r = 1, K = 100, h = 20, delta = 0.001, tl = 5000, e = 0.01, replicates = 10) {
  
  N <- 1     ## set the number of replicate webs to make
  matrixSize <- round(C/step)
  
  reps <- 1:replicates
  cSteps <- seq(step, C, step)
  x <- S/(C/step)
  sSteps <- seq(x, sum(rep(x, matrixSize)), x)
  
  for(i in reps) {
    for(C_step in cSteps) {
      for(S_step in sSteps) {
        print(c(i, C_step, S_step))
        tic("Time to Steps")
        z <- fullSim(S_step, C_step, p.m, p.e, s, alpha, xo, r, K, h, delta, tl, e)
        print(c(z[[3]][[9]], z[[2]]))
        debugGraph(z[[1]])
        toc()
      }
    }
  }
  return()
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
  steps <- nrow(out)
  return(list(out, mean(out[nrow(out),2:ncol(out)] > e), list(A.m, A.ep, A.en, A.c, h, K, r, s, steps, delta, tl, e, interaction_index, interaction_type)))
}