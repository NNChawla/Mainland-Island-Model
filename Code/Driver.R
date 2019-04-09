require("tidyverse")
require("reshape")
source("FoodWebFunctions.R")
source("Stability_Model.R")
source("MatrixGraphs.R")

S <- c(50, 100, 150, 200, 300)
C <- 0.5
step <- 0.05
containers <- list()

nI <- c(2, 5, 10)
sT <- c(10, 50, 100, 200)

for(i in 1:length(S)){
  step <- 0.05
  containers[[i]] <- tryCatch(CvNs(S[[i]], C, step), error=function(e) NA)
  decrease <- TRUE
  while(is.na(containers[[i]]) && !near(step, 0.0)){
    if(decrease){
      step <- step - 0.01
      containers[[i]] <- tryCatch(CvNs(S[[i]], C, step), error=function(e) NA)
    }
    else{
      containers[[i]] <- 0
    }
  }
  print(c(S[[i]], step))
}

datasets <- multiMatrix(containers, nI, sT)