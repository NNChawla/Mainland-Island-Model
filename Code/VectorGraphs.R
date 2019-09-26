#original return values
#list(frames, c(modelType, nI, stepTime, mainS), timeMatrix(frames), pathLDRatios)
#new return values
#list(frames, c(nI, mainS), timeMatrix(frames), pathLDRatios)

library(zoo)
library(ggtern)

fixDS <- function(n){
  print(n)
  fileName <- paste("ds", n, sep="")
  dsCopy <- readRDS(fileName)
  dsCopy$TimeMatrix <- tbTimeMatrix(timeMatrix(dsCopy$Frames))
  dsCopy$Data <- NULL
  data <- tbFrames(dsCopy$Frames, 5)
  dsCopy$Steps <- tbSteps(dsCopy$Steps, data)
  return(dsCopy)
}

mergeDatasets <- function(containers){
  dfs <- list()
  for (i in 1:length(containers)) {
    df <- containers[[i]]$Steps
    df <- add_column(df, Mutualism = rep(containers[[i]]$Parms[[4]], nrow(df)),
                     Exploitation = rep(containers[[i]]$Parms[[5]], nrow(df)),
                     Competition = rep(containers[[i]]$Parms[[6]], nrow(df)),
                     .before = 1)
    dfs[[i]] <- df
  }
  dfs <- bind_rows(dfs)
  dfs$eqSteps[is.na(dfs$eqSteps)] <- 1
  return(dfs)
}

probZ <- function(){
  probs <- list()
  count <- 1
  for(p.m in seq(0, 1, 0.1)){
    for(p.e in seq(0, 1, 0.1)) {
      p.c = 1.0 - (p.m + p.e)
      if(p.m+p.e+p.c == 1 && p.c > 0){
        probs[[count]] <- c(p.m, p.e, 1-(p.m+p.e))
        count <- count + 1
      }
    }
  }
  return(probs)
}

forge <- function(dfs) {
  probs <- probZ()
  tables <- vector(mode = "list", length = 5500)
  count <- 1
  
  for(i in seq(0.1, 1, 0.1)) {
    for(j in seq(20, 200, 20)) {
      for(k in 1:length(probs)) {
        tables[[count]] <- meanPath(dfs, connectance = i, richness = j, Mut = probs[[k]][1], Exp = probs[[k]][2], Comp = probs[[k]][3])
        count <- count + 1
        print(c(i, j, k))
      }
    }
  }
  return(tables)
}

meanPath <- function(mats, connectance = 0.1, richness = 20, Mut, Exp, Comp){
  mats <- mats %>% filter(C == connectance & N == richness & Mutualism == Mut & Exploitation == Exp & Competition == Comp)
  mats$eqSteps[is.na(mats$eqSteps)] = 1
  R <- c()
  val <- c()
  step <- c()
  maxVal <- 0
  for(i in 1:nrow(mats)){
    R <- c(R, rep(mats$R[i], mats$eqSteps[i]))
    val <- c(val, rep(mats$value[i], mats$eqSteps[i]))
  }
  
  for(i in 1:length(rle(R)$values)){
    x <- seq(1, rle(R)$lengths[[i]])
    if(length(x) > maxVal) {
      maxVal <- length(x)
    }
    step <- c(step, x)
  }
  
  mats <- as.tibble(data.frame(R = R, value = val, Step = step))
  values <- vector(mode="list", length=maxVal)
  steps <- 1:maxVal
  muts <- rep(Mut, maxVal)
  exps <- rep(Exp, maxVal)
  coms <- rep(Comp, maxVal)
  conns <- rep(connectance, maxVal)
  speciesN <- rep(richness, maxVal)
  
  for(i in steps){
    tmpMat <- mats %>% filter(Step == i)
    values[[i]] <- mean(tmpMat$value)
  }
  
  values <- unlist(values)
  mats <- as.tibble(data.frame(M=muts, E=exps, Comp=coms, C=conns, N=speciesN, Step=steps, Persistence=values))
  return(mats)
  # plt <- ggplot(mats, aes(x=Step, y=N))
  # plt <- plt + geom_line(size = 1.25) + xlab("Step Number") + ylab("Nisle")
  # plt <- plt + ggtitle(paste("Connectance", connectance, "Species #", richness))
  # 
  # return(plt)
}
  
  multiGraph <- function(multiMat, containers){
    massGraphs <- list()
  tGraphs <- list()
  pGraphs <- list()
  count <- 1
  for(i in 1:length(multiMat)){
    for(j in 1:length(multiMat[[i]])){
      for(k in 1:length(multiMat[[i]][[j]])){
        if(is.na(multiMat[[i]][[j]][[k]])){
          next
        }
        massMat <- multiMat[[i]][[j]][[k]]
        meanData <- containers[[i]]$mean
        S <- containers[[i]]$parms[[1]]
        fName <- paste(S, "-", massMat[[2]][[1]], "-", massMat[[2]][[2]], "-", massMat[[2]][[3]], ".png", sep="")
        massGraphs[[count]] <- massGraph(massMat, paths=c("All"))
        tGraphs[[count]] <- t50graph(massMat, rownames(meanData), colnames(meanData))
        pGraphs[[count]] <- pathGraph(massMat)
        ggsave(filename=paste("mass-", fName, sep=""), plot=massGraphs[[count]], dpi=320, width=20, height=10)
        ggsave(filename=paste("t50-", fName, sep=""), plot=tGraphs[[count]], dpi=320, width=20, height=10)
        ggsave(filename=paste("path-", fName, sep=""), plot=pGraphs[[count]], dpi=320, width=20, height=10)
        
        count <- count + 1
      }
    }
  }
  return(list(massGraphs, tGraphs, pGraphs))
}

simplexGraph <- function(containers, connectance, numSpecies) {
  filterFrame <- function(x, cn, ns) {
    x <- x$TimeMatrix %>% filter(C == cn & N == ns)
    
    #return either mean of final persistence or of total steps
    return(mean(x$t100))
  }
  
  pm <- c()
  pe <- c()
  pc <- c()
  
  count <- 1
  for(p.m in seq(0, 1, 0.1)){
    for(p.e in seq(0, 1, 0.1)) {
      p.c = 1.0 - (p.m + p.e)
      if(p.m+p.e+p.c == 1 && p.c >= 0){
        pm[[count]] <- p.m
        pe[[count]] <- p.e
        pc[[count]] <- 1 - (p.m + p.e)
        count <- count + 1
      }
    }
  }
  
  
  values <- unlist(lapply(containers, filterFrame, cn=connectance, ns=numSpecies))
  
  df <- data.frame(M = pm, E = pe, C = pc, N = rep(0, length(pc)))
  
  count <- 1
  for(i in 1:length(pm)) {
    if(df[i,3] != 0 && df[i,3] != 1) {
      df[i,4] = values[count]
      count <- count + 1
    }
  }
  
  palette <- c( "#FF9933", "#002C54", "#3375B2", "#CCDDEC", "#BFBFBF", "#000000")
  simplex <- ggtern(data=df, aes(E, C, M)) + 
    geom_hex_tern(aes(value=log(N)), color="black", alpha=1)
  simplex <- simplex + scale_fill_gradient(name="Scale", low=palette[2], high=palette[4])
  
  #simplex <- simplex + Tarrowlab("Competition") +
  #  Larrowlab("Exploitation") + Rarrowlab("Mutualism") + theme_showarrows()
  
  simplex <- simplex + labs(title=paste(connectance, numSpecies))
  
  return(simplex)
}

simplexArrange <- function(containers, N) {
  
  graphs <- vector(mode = "list", length = 10)
  graphs[[1]] <- simplexGraph(containers, 0.1, N)
  graphs[[2]] <- simplexGraph(containers, 0.2, N)
  graphs[[3]] <- simplexGraph(containers, 0.3, N)
  graphs[[4]] <- simplexGraph(containers, 0.4, N)
  graphs[[5]] <- simplexGraph(containers, 0.5, N)
  graphs[[6]] <- simplexGraph(containers, 0.6, N)
  graphs[[7]] <- simplexGraph(containers, 0.7, N)
  graphs[[8]] <- simplexGraph(containers, 0.8, N)
  graphs[[9]] <- simplexGraph(containers, 0.9, N)
  graphs[[10]] <- simplexGraph(containers, 1.0, N)
  plot <- do.call("grid.arrange", c(graphs, ncol=2, as.table=FALSE))
  return(plot)
}

lineArrange <- function(containers, N) {
  
}

#Line Graph: average replicates of timeStep extended line graphs for each container by padding shorter ones with zeros at the end?
#Mass mean lines for every container for every CvN pairing into singular 55 line graphs?

pdfSimplex <- function(containers) {
  N <- seq(20, 200, 20)
  for(i in 1:length(N)) {
    pdf(paste("t100-log-", N[[i]], ".pdf", sep=""))
    print(simplexArrange(containers, N[[i]]))
    dev.off()
  }
}