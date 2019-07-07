source("QianModel.R")
source("VectorModel.R")
source("VectorMatrixFunctions.R")

S <- 200
C <- 1.0
step <- 0.1
reps <- 5
containers <- list()

count <- 1
for(p.m in seq(0, 1, 0.1)){
  for(p.e in seq(0, 1, 0.1)) {
    p.c = 1.0 - (p.m + p.e)
    if(p.m+p.e+p.c == 1 && p.c > 0){
      containers[[count]] <- VectorCvNs(S, C, step, p.m, p.e, replicates = reps)
      count <- count + 1
    }
  }
}

datasets <- list()
for(i in 1:length(containers)){
  datasets[[i]] <- VectorMetaMatrix(containers[[i]], replicates = 5)
}