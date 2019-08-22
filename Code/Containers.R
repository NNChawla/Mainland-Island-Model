#!/usr/bin/env Rscript

rm(list = ls())

source("QianModel.R")
source("VectorModel.R")
source("VectorMatrixFunctions.R")

S <- 200
C <- 1.0
step <- 0.1
reps <- 5

vcn <- function(x){
  out <- VectorCvNs(S, C, step, x[1], x[2], replicates = reps)
  print(format(object.size(out), units="auto"))
  return(out)
}

probs <- list()
count <- 1
for(p.m in seq(0, 1, 0.1)){
  for(p.e in seq(0, 1, 0.1)) {
    p.c = 1.0 - (p.m + p.e)
    if(p.m+p.e+p.c == 1 && p.c > 0){
      probs[[count]] <- c(p.m, p.e)
      count <- count + 1
    }
  }
}
containers <- lapply(probs, vcn)
print("Containers Done")
print(format(object.size(containers), units="auto"))
saveRDS(containers, "containers.rds")