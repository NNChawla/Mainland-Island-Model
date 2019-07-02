# VectorMetaMatrix <- function(container, nI = 5, graphStep = 1, replicates = 10, e = 0.01, kl = 5000) {
#   
#   mainS <- container$parms[[1]]
#   meanData <- container$mean
#   mainlands <- container$communities
#   interactions <- container$interactions
#   nStarMatrix <- data.frame(matrix(nrow=nrow(meanData), ncol=ncol(meanData)))
#   for(i in 1:ncol(nStarMatrix)){
#     nStarMatrix[i] <- meanData[,i]*as.numeric(colnames(meanData)[i])
#   }
#   colnames(nStarMatrix) <- colnames(meanData)
#   rownames(nStarMatrix) <- rownames(meanData)
#   
#   rownums <- 1:nrow(nStarMatrix)
#   colnums <- 1:ncol(nStarMatrix)
#   
#   cs <- as.numeric(rownames(meanData))
#   ss <- as.numeric(colnames(meanData))
#   
#   lis <- list()
#   stepCounter <- 1
#   
#   for(i in cs) {
#     for(j in ss) {
#       rpNum <- sample(length(mainlands), 1)
#       c <- which(cs == i)
#       s <- which(ss == j)
#       community <- mainlands[[rpNum]][[c]][[s]]
#       interaction <- interactions[[rpNum]][[c]][[s]]
#       if(is.null(community)) {
#         list("NULL", "NULL")
#       }
#       else {
#         
#         mainLiving <- sum(community[nrow(community), 2:ncol(community)] > e)
#         stepSize <- graphStep
#         subset <- list()
#         frames <- list()
#         pldRs <- list()
#         
#         for(k in 1:replicates) {
#           pldRs[[k]] <- VectorPath(community, interaction, nI, i, e, kl)
#           subset[k] <- list(pldRs[[k]])
#           z <- c()
#           for(l in 1:length(subset[[k]])){
#             z <- c(z, lengths(subset[[k]][[l]]["Living"], use.names=FALSE))
#           }
#           z <- z/mainLiving
#           z <- data.frame(z)
#           colnames(z) <- k
#           z["Step"] <- c(1:nrow(z))
#           if(is.na(subset[k]))
#             z[[1]] <- NA
#           frames[[k]] <- z[seq(1, nrow(z), stepSize), ]
#         }
#         
#         frames <- Reduce(function(x, y) merge(x=x, y=y, by="Step", all.y = TRUE, all.x = TRUE), frames)
#         colnames(frames) <- c("Step", paste(i, j, 1:replicates))
#         list(frames, pldRs)
#       }
#     }
#   }
#   
#   meanMat <- list()
#   pathLDRatios <- list()
#   
#   for(i in 1:length(out)){
#     meanMat[[i]] <- list()
#     pathLDRatios[[i]] <- list()
#     for(j in 1:length(out[[i]])){
#       meanMat[[i]][[j]] <- out[[i]][[j]][[1]]
#       pathLDRatios[[i]][[j]] <- out[[i]][[j]][[2]]
#     }
#   }
#   
#   frames <- list()
#   for(i in 1:length(meanMat)){
#     frames[[i]] <- Reduce(function(x, y) merge(x=x, y=y, by="Step", all.y = TRUE), meanMat[[i]])
#   }
#   frames <- Reduce(function(x, y) merge(x=x, y=y, by="Step", all.y = TRUE), frames)
#   
#   for(i in 1:length(pathLDRatios)){
#     for(j in 1:length(pathLDRatios[[i]])){
#       if(sum(is.na(pathLDRatios[[i]][[j]][[1]])) == 1) {
#         pathLDRatios[[i]][[j]] <- NA
#         next
#       }
#       for(k in 1:replicates){
#         z <- c()
#         for(l in 1:length(pathLDRatios[[i]][[j]][[k]])){
#           live <- length(pathLDRatios[[i]][[j]][[k]][[l]][["Living"]])
#           dead <- length(pathLDRatios[[i]][[j]][[k]][[l]][["Before"]])
#           z <- c(z, live/dead)
#         }
#         pathLDRatios[[i]][[j]][[k]] <- z
#       }
#     }
#   }
#   
#   return(list(frames, c(nI, mainS), timeMatrix(frames), pathLDRatios))
# }