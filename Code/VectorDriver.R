#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
containerNum <- as.numeric(args[1])

source("QianModel.R")
source("VectorModel.R")
source("VectorMatrixFunctions.R")

reps <- 5

vmm <- function(x){
  out <- VectorMetaMatrix(x, replicates = reps)
  print(format(object.size(out), units="auto"))
  return(out)
}

container <- readRDS("containers.rds")[[containerNum]]

saveRDS(vmm(container), paste("ds", containerNum, sep=""))
print(paste(containerNum, "Complete"))