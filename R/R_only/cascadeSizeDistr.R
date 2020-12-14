library(igraph)
library(tictoc)
#cascade size distribution on tree:
#input: igraph object gtree with ICM paremters: V(gtree)$pInit are initial activation probabilities,
#network weights: E(gtree)$weight infection probabilities 
treeICM <- function(gtree){
  N <- length(V(gtree))
  V(gtree)$name <- apply(as.array(1:N), 1, function(x) paste("n",x,sep = ""))
  Names <- V(gtree)$name
  pInit <- V(gtree)$pInit
  names(pInit) <- Names
  deg <- degree(gtree,mode = "in")
  kmax <- max(deg)
  root <- V(gtree)$name[which(deg == kmax)[1]]
  tree <- list(c(deg[root]))
  wr <- c(0)
  names(wr) <- root
  treeW <- list(wr)
  treeWback <- list(wr)
  treePinit <- list(c(pInit[root]))
  nodes <- root
  l <- 1
  while((length(nodes) > 0 )){
    l <- l + 1
    nbNodes <- length(nodes)
    newNodes <- c()
    newDegrees <- c()
    newW <- c()
    newWback <- c()
    newPinit <- c()
    for(i in 1:nbNodes){
      nb <- names(neighbors(gtree, nodes[i], mode = 1))
      deg <- degree(gtree,mode = "in") 
      pInit <- V(gtree)$pInit
      names(pInit) <- names(V(gtree))
      newDegrees <- c(newDegrees,deg[nb])
      newPinit <- c(newPinit,pInit[nb])
      wVec <- as.vector(rbind(nb,rep(nodes[i],length(nb))))
      wVecBack <- as.vector(rbind(rep(nodes[i],length(nb)),nb))
      wN <- E(gtree)$weight[get.edge.ids(gtree,wVec)]
      wNback <- E(gtree)$weight[get.edge.ids(gtree,wVecBack)]
      names(wN) <- nb
      names(wNback) <- nb
      newW <- c(newW, wN)
      newWback <- c(newWback, wNback)
      gtree <- delete.vertices(gtree,nodes[i])
      newNodes <- c(newNodes,as.vector(nb[deg[nb]>1]))
    }
    tree[[l]] <- newDegrees
    treeW[[l]] <- newW
    treeWback[[l]] <- newWback
    treePinit[[l]] <- newPinit
    
    nodes <- newNodes
  }
  #treeW: weight of node n is w_pn, treeWback refers to w_np
  levelTree <- length(tree)
  nbrNodesLevel <- unlist(lapply(tree, length))
  maxLevel <- max(nbrNodesLevel)
  Napprox <- N + 1 
  #compute distr:
  zeroPadDim <- nextn(Napprox*2 - 1) 
  transFFT <- exp(-2*pi/zeroPadDim*(0:(zeroPadDim-1))*1i)
  
  for(l in levelTree:2){
    pBnAll <- array(data = 0, dim = c(nbrNodesLevel[l],zeroPadDim))
    pAzeroAll <- array(data = 0, dim = c(nbrNodesLevel[l],zeroPadDim))
    pAsumAll <- array(data = 0, dim = c(nbrNodesLevel[l],zeroPadDim))
    
    childrenCount <- 0
    for(i in 1:nbrNodesLevel[l]){
      #belong to level l  
      k <- tree[[l]][i]
      if(is.null(k) == FALSE){
        if(k > 1){
          children <-  (childrenCount + 1):(childrenCount + k - 1)
          childrenCount <- childrenCount + k - 1
          if(k>2){
            #parent does not fail:
            pwf <- apply(pBnOld[children,],2,prod)*(1-treePinit[[l]][i])
            #later failure of node but no child triggers failure, the parent does
            plfpTmin1 <- apply(pAzeroOld[children,],2,prod)*transFFT*(1-treePinit[[l]][i])
            #later failure of node:
            pfTmin1 <- apply(pAsumOld[children,],2,prod)*transFFT - plfpTmin1
          }else{
            #parent does not fail:
            pwf <- pBnOld[children,]*(1-treePinit[[l]][i])
            #later failure of node but no child triggers failure, the parent does
            plfpTmin1 <- pAzeroOld[children,]*transFFT*(1-treePinit[[l]][i])
            pfTmin1 <- pAsumOld[children,]*transFFT - plfpTmin1
          }
          pBn <- pwf + pfTmin1*(1-treeW[[l]][i])
          pAzero <- pwf*(1-treeWback[[l]][i]) + pfTmin1*(1-treeW[[l]][i]) + plfpTmin1*treeWback[[l]][i]
          pAsum <-   pAzero + pfTmin1*(treeW[[l]][i])
        } else{
          pBn <-(1-treePinit[[l]][i]) + treePinit[[l]][i]*(1-treeW[[l]][i])*transFFT
          pAzero <- (1-treePinit[[l]][i])*(1-treeWback[[l]][i]) + transFFT*((1-treePinit[[l]][i])*treeWback[[l]][i] + treePinit[[l]][i]*(1-treeW[[l]][i]))
          pAsum <-   pAzero + treePinit[[l]][i]*treeW[[l]][i]*transFFT
        }  
        pBnAll[i,] <- pBn
        pAzeroAll[i,] <- pAzero
        pAsumAll[i,] <- pAsum
      }
    }
    pBnOld <- pBnAll
    pAzeroOld <- pAzeroAll
    pAsumOld <- pAsumAll
  }
  #root
  k <- tree[[1]][1]
  children <-  1:k
  if(k>1){
    #parent does not fail:
    pwf <- apply(pBnOld[children,],2,prod)*(1-treePinit[[1]][1])
    #later failure of node but no child triggers failure, the parent does
    plfpTmin1 <- apply(pAzeroOld[children,],2,prod)*transFFT*(1-treePinit[[1]][1])
    #later failure of node:
    pfTmin1 <- apply(pAsumOld[children,],2,prod)*transFFT - plfpTmin1
  }else{
    #parent does not fail:
    pwf <- pBnOld[children,]*(1-treePinit[[1]][1])
    #later failure of node but no child triggers failure, the parent does
    plfpTmin1 <- pAzeroOld[children,]*transFFT*(1-treePinit[[1]][1])
    #later failure of node:
    pfTmin1 <- pAsumOld[children,]*transFFT - plfpTmin1
  }
  distrFFT <- pwf + pfTmin1
  distr <- Re(fft(distrFFT,inverse = TRUE)[1:(N+1)])/zeroPadDim
  
  return(distr) 
}
#helper function to compute cascade size distribution on general network:
#converts relevant variables in tree format for easy message passing from leaves to root
treesTree <- function(gtree){
  N <- length(V(gtree))
  V(gtree)$name <- apply(as.array(1:N), 1, function(x) paste("n",x,sep = ""))#g$tree.nodes
  Names <- V(gtree)$name
  pInit <- V(gtree)$pInit
  names(pInit) <- Names
  deg <- degree(gtree,mode = "in")
  kmax <- max(deg)
  root <- V(gtree)$name[which(deg == kmax)[1]]
  tree <- list(c(deg[root]))
  wr <- c(0)
  names(wr) <- root
  treeW <- list(wr)
  treeWback <- list(wr)
  treePinit <- list(c(pInit[root]))
  nodes <- root
  l <- 1
  while((length(nodes) > 0 )){
    l <- l + 1
    nbNodes <- length(nodes)
    newNodes <- c()
    newDegrees <- c()
    newW <- c()
    newWback <- c()
    newPinit <- c()
    for(i in 1:nbNodes){
      nb <- names(neighbors(gtree, nodes[i], mode = 1))
      deg <- degree(gtree,mode = "in") 
      pInit <- V(gtree)$pInit
      names(pInit) <- names(V(gtree))
      newDegrees <- c(newDegrees,deg[nb])
      newPinit <- c(newPinit,pInit[nb])
      wVec <- as.vector(rbind(nb,rep(nodes[i],length(nb))))
      wVecBack <- as.vector(rbind(rep(nodes[i],length(nb)),nb))
      wN <- E(gtree)$weight[get.edge.ids(gtree,wVec)]
      wNback <- E(gtree)$weight[get.edge.ids(gtree,wVecBack)]
      names(wN) <- nb
      names(wNback) <- nb
      newW <- c(newW, wN)
      newWback <- c(newWback, wNback)
      gtree <- delete.vertices(gtree,nodes[i])
      newNodes <- c(newNodes,as.vector(nb[deg[nb]>1]))
    }
    tree[[l]] <- newDegrees
    treeW[[l]] <- newW
    treeWback[[l]] <- newWback
    treePinit[[l]] <- newPinit
    nodes <- newNodes
  }
  return(list(tree, treeW, treeWback, treePinit))
}
#cascade size distribution on tree (memory reduced) with treeStruct as input instead of igraph object
treeICMred <- function(treeStruct, N){
  tree <- treeStruct[[1]] 
  treeW <- treeStruct[[2]]
  treeWback <- treeStruct[[3]]
  treePinit <- treeStruct[[4]]
  #treeW: weight of node n is w_pn, treeWback refers to w_np
  levelTree <- length(tree)
  nbrNodesLevel <- unlist(lapply(tree, length))
  maxLevel <- max(nbrNodesLevel)
  Napprox <- N + 1
  #compute distr:
  tic("cascade size distr")
  zeroPadDim <- Napprox+1 #nextn(Napprox*2 - 1)
  transFFT <- exp(-2*pi/zeroPadDim*(0:(zeroPadDim-1))*1i)
  
  for(l in levelTree:2){
    pBnAll <- array(data = 0, dim = c(nbrNodesLevel[l],zeroPadDim))
    pAzeroAll <- array(data = 0, dim = c(nbrNodesLevel[l],zeroPadDim))
    pAsumAll <- array(data = 0, dim = c(nbrNodesLevel[l],zeroPadDim))
    
    childrenCount <- 0
    for(i in 1:nbrNodesLevel[l]){
      #belong to level l  
      k <- tree[[l]][i]
      if(is.null(k) == FALSE){
        if(k > 1){
          children <-  (childrenCount + 1):(childrenCount + k - 1)
          childrenCount <- childrenCount + k - 1
          if(k>2){
            #parent does not fail:
            pwf <- apply(pBnOld[children,],2,prod)*(1-treePinit[[l]][i])
            #initial failure of node:
            plfpTmin1 <- apply(pAzeroOld[children,],2,prod)*transFFT*(1-treePinit[[l]][i])
            #later failure of node:
            pfTmin1 <- apply(pAsumOld[children,],2,prod)*transFFT - plfpTmin1
          }else{
            #parent does not fail:
            pwf <- pBnOld[children,]*(1-treePinit[[l]][i])
            #initial failure of node:
            plfpTmin1 <- pAzeroOld[children,]*transFFT*(1-treePinit[[l]][i])
            #later failure of node:
            pfTmin1 <- pAsumOld[children,]*transFFT - plfpTmin1
          }
          pBn <- pwf + pfTmin1*(1-treeW[[l]][i])
          #pAzeroInit <- pwf*(1-treeWback[[l]][i]) + (pifTmin1 + plfTmin1) + plfpTmin1*treeWback[[l]][i]
          pAzero <- pwf*(1-treeWback[[l]][i]) + pfTmin1*(1-treeW[[l]][i]) + plfpTmin1*treeWback[[l]][i]
          pAsum <-   pAzero + pfTmin1*(treeW[[l]][i])
          #pAzeroInit = pAsum!
        } else{
          pBn <-(1-treePinit[[l]][i]) + treePinit[[l]][i]*(1-treeW[[l]][i])*transFFT
          pAzero <- (1-treePinit[[l]][i])*(1-treeWback[[l]][i]) + transFFT*((1-treePinit[[l]][i])*treeWback[[l]][i] + treePinit[[l]][i]*(1-treeW[[l]][i]))
          pAsum <-   pAzero + treePinit[[l]][i]*treeW[[l]][i]*transFFT
        }  
        pBnAll[i,] <- pBn
        pAzeroAll[i,] <- pAzero
        pAsumAll[i,] <- pAsum
      }
    }
    pBnOld <- pBnAll
    pAzeroOld <- pAzeroAll
    pAsumOld <- pAsumAll
  }
  rm(pBnAll, pAzeroAll, pAsumAll)
  gc()
  #root
  k <- tree[[1]][1]
  children <-  1:k
  if(k>1){
    #parent does not fail:
    pwf <- apply(pBnOld[children,],2,prod)*(1-treePinit[[1]][1])
    pfTmin1 <- (apply(pAsumOld[children,],2,prod) - (1-treePinit[[1]][1])*apply(pAzeroOld[children,],2,prod))*transFFT 
  }else{
    #parent does not fail:
    pwf <- pBnOld[children,]*(1-treePinit[[1]][1])
    pfTmin1 <- (pAsumOld[children,] - (1-treePinit[[1]][1])*pAzeroOld[children,])*transFFT
  }
  distrFFT <- pwf + pfTmin1
  rm(pwf, pfTmin1)
  gc()
  distr <- Re(fft(distrFFT,inverse = TRUE)[1:(N+1)])/zeroPadDim
  time <- toc()
  return(list(distr,time)) 
}

### Tree Distribution Approximation ###
#Belief Propagation to compute marginal activation probabilities
BP_ICM <-function(g,deg,N){
  #pij[[i]][j]: P(s_j = 1 ||  s_i = 0) * w_{ji}
  pInit <- V(g)$pInit
  nb <- neighbors(g, 1, mode = "in")
  pij <- V(g)$pInit[nb]
  names(pij) <- V(g)$name[nb]
  #wij[[i]][j]: w_{ji} 
  wvec <- as.vector(rbind(names(nb),rep("n1",length(nb)))) 
  wij <- E(g)$weight[get.edge.ids(g,wvec)]
  names(wij) <- names(nb)
  indVec <- 1:N
  names(indVec) <- V(g)$name
  pijList <- list(pij*wij)
  wijList <- list(wij)
  for(i in 2:N){
    nb <- neighbors(g, i, mode = "in")
    pij <- V(g)$pInit[nb]
    names(pij) <- V(g)$name[nb]
    
    wvec <- as.vector(rbind(names(nb),rep(V(g)$name[i],length(nb)))) 
    wij <- E(g)$weight[get.edge.ids(g,wvec)]
    names(wij) <- names(nb)
    wijList[[i]] <- list(wij)
    pijList[[i]] <- list(pij*wij)
  }
  
  pijListOld <- pijList
  for(t in 1:10){
    prodVec <- unlist(sapply(1:N, function(x, y, p) (1-p[x])*prod(1-unlist(y[[x]])), y=pijList, p = pInit))
    for(i in 1:N){
      pij <- unlist(pijList[[i]])
      wij <- unlist(wijList[[i]])
      nb <- neighbors(g, i, mode = "in")
      nbInd <- indVec[nb]
      for(j in 1:deg[i]){
        pij[j] <- (1- prodVec[nb[j]]/(1-(unlist(pijListOld[[nbInd[j]]])[paste("n",i,sep="")])))*wij[j]
      }
      pijList[[i]] <- list(pij)
    }
    pijListOld <- pijList
  }
  #avg failure prob:
  #failProb <- 1 - unlist(sapply(1:N, function(x, y, p) (1-p[x])*prod(1-unlist(y[[x]])), y=pijList, p = pInit))
  return(pijList)
}
##
#Compute Maximum Spanning Tree and updated ICM for TDA:
trees <- function(g,FailProbBeforeNode,deg,N){
  set.seed(20)
  pInit <- V(g)$pInit
  names(pInit) <- V(g)$name
  gud <- as.undirected(g)
  gtree <- mst(gud, weights = -E(gud)$weight, algorithm = "prim")
  nbTree <- neighborhood(gtree,1,V(gtree),"in") 
  nbOrig <- neighborhood(g,1,V(g),"in") 
  degOrig <- deg
  delNeighbors <- list()
  for(i in 1:N){
    delNeighbors <- append(delNeighbors,list(setdiff(nbOrig[[i]],nbTree[[i]])))
  }
  deg <- degree(gtree)
  vec <- degOrig-deg
  initFailProb <- pInit 
  names(initFailProb) <- V(g)$name
  for(i in 1:N){
    if(vec[i]>0){
      initFailProb[i] <- 1- (1-pInit[i])*prod(1-unlist(FailProbBeforeNode[[i]])[V(g)$name[unlist(delNeighbors[i])]])
    }
  }
  pInit <- initFailProb
  V(gtree)$pInit <- pInit
  rm(initFailProb)
  
  ### Build data structure of MST
  deg <- degree(gtree,mode = "in")
  kmax <- max(deg)
  #root <- V(gtree)$name[which(deg == kmax)[1]]
  root <-  V(gtree)$name[1]
  tree <- list(c(deg[root]))
  wr <- c(0)
  names(wr) <- root
  treeW <- list(wr)
  treeWback <- list(wr)
  treePinit <- list(c(pInit[root]))
  nodes <- root
  l <- 1
  while((length(nodes) > 0 )){
    l <- l + 1
    nbNodes <- length(nodes)
    newNodes <- c()
    newDegrees <- c()
    newW <- c()
    newWback <- c()
    newPinit <- c()
    for(i in 1:nbNodes){
      nb <- names(neighbors(gtree, nodes[i], mode = "in"))
      deg <- degree(gtree,mode = "in") 
      pInit <- V(gtree)$pInit
      names(pInit) <- names(V(gtree))
      newDegrees <- c(newDegrees,deg[nb])
      newPinit <- c(newPinit,pInit[nb])
      
      wVec <- as.vector(rbind(nb,rep(nodes[i],length(nb))))
      wVecBack <- as.vector(rbind(rep(nodes[i],length(nb)),nb))
      wN <- E(g)$weight[get.edge.ids(g,wVec)]
      wNback <- E(g)$weight[get.edge.ids(g,wVecBack)]
      names(wN) <- nb
      names(wNback) <- nb
      newW <- c(newW, wN)
      newWback <- c(newWback, wNback)
      
      gtree <- delete.vertices(gtree,nodes[i])
      newNodes <- c(newNodes,as.vector(nb[deg[nb]>1]))
    }
    tree[[l]] <- newDegrees
    treeW[[l]] <- newW
    treeWback[[l]] <- newWback
    treePinit[[l]] <- newPinit
    nodes <- newNodes
  }
  return(list(tree,treePinit,treeW,treeWback))
}
#helper functions for convolutions with FFT:
embed <- function(x,N){
  #increase size
  n = length(x)
  return(fft(c(fft(x, inverse = TRUE)/n, rep(0, N-n))))
}
realemb <- function(x){
  return(fft(x, inverse = TRUE)/length(x))
}
#cascade size distribution for general network g (igraph object)
#input: igraph object gtree with ICM paremters: V(gtree)$pInit are initial activation probabilities,
#network weights: E(gtree)$weight infection probabilities 
ICMmem <- function(g){
  N <- length(V(g))
  V(g)$name <- apply(as.array(1:N), 1, function(x) paste("n",x,sep = ""))
  pInit <- V(g)$pInit
  deg <- degree(g,mode = "in")
  tic("BP")
  FailProbBeforeParent <- BP_ICM(g,deg,N) 
  tBP <- toc()
  treeStruct <- trees(g,FailProbBeforeParent,deg,N) #includes MST calculation etc.
  tree <- treeStruct[[1]]
  treePinit <- treeStruct[[2]]
  treeW <- treeStruct[[3]]
  treeWback <- treeStruct[[4]]
  rm(g, deg)
  gc()
  
  levelTree <- length(tree)
  nbrNodesLevel <- unlist(lapply(tree, length))
  maxLevel <- max(nbrNodesLevel)
  Napprox <- N + 1 #sum(nbrNodesLevel) + 1 #N
  tic("distr ICM mem")
  #compute distr:
  zeroPadDim <- Napprox + 1 #nextn(Napprox*2 - 1)
  transFFT <- exp(-2*pi/zeroPadDim*(0:(zeroPadDim-1))*1i)
  for(l in levelTree:2){
    pBnAll <- array(data = 0, dim = c(nbrNodesLevel[l],zeroPadDim))
    pAzeroAll <- array(data = 0, dim = c(nbrNodesLevel[l],zeroPadDim))
    pAsumAll <- array(data = 0, dim = c(nbrNodesLevel[l],zeroPadDim))
    
    childrenCount <- 0
    for(i in 1:nbrNodesLevel[l]){
      #belong to level l  
      k <- tree[[l]][i]
      if(is.null(k) == FALSE){
        if(k > 1){
          children <-  (childrenCount + 1):(childrenCount + k - 1)
          childrenCount <- childrenCount + k - 1
          if(k>2){
            #parent does not fail:
            pwf <- apply(pBnOld[children,],2,prod)*(1-treePinit[[l]][i])
            #initial failure of node:
            #later failure of node but no child triggers failure, the parent does
            plfpTmin1 <- apply(pAzeroOld[children,],2,prod)*transFFT*(1-treePinit[[l]][i])
            #later failure of node:
            pfTmin1 <- apply(pAsumOld[children,],2,prod)*transFFT - plfpTmin1
            #
          }else{
            #parent does not fail:
            pwf <- pBnOld[children,]*(1-treePinit[[l]][i])
            #later failure of node but no child triggers failure, the parent does
            plfpTmin1 <- pAzeroOld[children,]*transFFT*(1-treePinit[[l]][i])
            #later failure of node:
            pfTmin1 <- pAsumOld[children,]*transFFT - plfpTmin1
          }
          pBn <- pwf + pfTmin1*(1-treeW[[l]][i])
          pAzero <- pwf*(1-treeWback[[l]][i]) + pfTmin1*(1-treeW[[l]][i]) + plfpTmin1*treeWback[[l]][i]
          pAsum <-   pAzero + pfTmin1*(treeW[[l]][i])
        } else{
          pBn <-(1-treePinit[[l]][i]) + treePinit[[l]][i]*(1-treeW[[l]][i])*transFFT
          pAzero <- (1-treePinit[[l]][i])*(1-treeWback[[l]][i]) + transFFT*((1-treePinit[[l]][i])*treeWback[[l]][i] + treePinit[[l]][i]*(1-treeW[[l]][i]))
          pAsum <-   pAzero + treePinit[[l]][i]*treeW[[l]][i]*transFFT
        } 
        pBnAll[i,] <- pBn
        pAzeroAll[i,] <- pAzero
        pAsumAll[i,] <- pAsum
      }
    }
    pBnOld <- pBnAll
    pAzeroOld <- pAzeroAll
    pAsumOld <- pAsumAll
  }
  rm(pBnAll, pAzeroAll, pAsumAll)
  gc()
  #root
  k <- tree[[1]][1]
  children <-  1:k
  if(k>1){
    #parent does not fail:
    pwf <- apply(pBnOld[children,],2,prod)*(1-treePinit[[1]][1])
    #parent does fail:
    pfTmin1 <- (apply(pAsumOld[children,],2,prod) - (1-treePinit[[1]][1])*apply(pAzeroOld[children,],2,prod))*transFFT 
  }else{
    pwf <- pBnOld[children,]*(1-treePinit[[1]][1])
    pfTmin1 <- (pAsumOld[children,] - (1-treePinit[[1]][1])*pAzeroOld[children,])*transFFT
  }
  distrFFT <- pwf + pfTmin1
  distr <- Re(fft(distrFFT,inverse = TRUE)[1:(N+1)])/zeroPadDim
  distr[1] <- prod(1-pInit)
  distr[-1] <- distr[-1]/sum(distr[-1])*(1-distr[1])
  tD <- toc()
  return(list(distr, tBP, tD))   
}

### conditional failure probabilities ###
condFailProb_ICM <- function(g){
  N <- length(V(g))
  V(g)$name <- apply(as.array(1:N), 1, function(x) paste("n",x,sep = ""))#g$tree.nodes
  Names <- V(g)$name
  pInit <- V(g)$pInit
  names(pInit) <- Names
  deg <- degree(g,mode = "in")
  FailProbBeforeParent <- BP_ICM(g,deg,N) 
  treeStruct <- trees(g,FailProbBeforeParent,deg,N) 
  tree <- treeStruct[[1]]
  treePinit <- treeStruct[[2]]
  treeW <- treeStruct[[3]]
  treeWback <- treeStruct[[4]]
  parentTree <- tree
  levelTree <- length(tree)
  nbrNodesLevel <- unlist(lapply(tree, length))
  maxLevel <- max(nbrNodesLevel)
  Napprox <- N + 1
  #compute distr:
  zeroPadDim <- nextn(Napprox*2 - 1) #nextn(Napprox) #nextn(Napprox*2 - 1)
  transFFT <- exp(-2*pi/zeroPadDim*(0:(zeroPadDim-1))*1i)
  #memorize pBn for every node:
  matwf <- array(data = 0, dim = c(N, zeroPadDim))
  matlfTmin1 <- array(data = 0, dim = c(N, zeroPadDim))
  matfTmin1 <- array(data = 0, dim = c(N, zeroPadDim))
  rownames(matwf) <- Names
  rownames(matlfTmin1) <- Names
  rownames(matfTmin1) <- Names
  
  for(l in levelTree:2){
    pBnAll <- array(data = 0, dim = c(nbrNodesLevel[l],zeroPadDim))
    pAzeroAll <- array(data = 0, dim = c(nbrNodesLevel[l],zeroPadDim))
    pAsumAll <- array(data = 0, dim = c(nbrNodesLevel[l],zeroPadDim))
    
    childrenCount <- 0
    for(i in 1:nbrNodesLevel[l]){
      #belong to level l  
      k <- tree[[l]][i]
      node <- names(tree[[l]][i])
      if(is.null(k) == FALSE){
        if(k > 1){
          children <-  (childrenCount + 1):(childrenCount + k - 1)
          childrenCount <- childrenCount + k - 1
          parentVec <- rep(node, k-1)
          names(parentVec) <- names(tree[[l+1]][children])
          parentTree[[l+1]][children] <- parentVec
          if(k>2){
            #parent does not fail:
            pwf <- apply(pBnOld[children,],2,prod)*(1-treePinit[[l]][i])
            #later failure of node but no child triggers failure, the parent does
            plfpTmin1 <- apply(pAzeroOld[children,],2,prod)*transFFT*(1-treePinit[[l]][i])
            #later failure of node:
            pfTmin1 <- apply(pAsumOld[children,],2,prod)*transFFT - plfpTmin1
          }else{
            #parent does not fail:
            pwf <- pBnOld[children,]*(1-treePinit[[l]][i])
            #later failure of node but no child triggers failure, the parent does
            plfpTmin1 <- pAzeroOld[children,]*transFFT*(1-treePinit[[l]][i])
            #later failure of node:
            pfTmin1 <- pAsumOld[children,]*transFFT - plfpTmin1
          }
          matwf[node, ] <- pwf
          matlfTmin1[node, ] <-  plfpTmin1
          matfTmin1[node, ] <-  pfTmin1 
          pBn <- pwf + pfTmin1*(1-treeW[[l]][i])
          pAzero <- pwf*(1-treeWback[[l]][i]) + pfTmin1*(1-treeW[[l]][i]) + plfpTmin1*treeWback[[l]][i]
          pAsum <-   pAzero + pfTmin1*(treeW[[l]][i])
        } else{
          pBn <- (1-treePinit[[l]][i]) + treePinit[[l]][i]*(1-treeW[[l]][i])*transFFT
          pAzero <- (1-treePinit[[l]][i])*(1-treeWback[[l]][i]) + transFFT*((1-treePinit[[l]][i])*treeWback[[l]][i] + treePinit[[l]][i]*(1-treeW[[l]][i]))
          pAsum <-   pAzero + treePinit[[l]][i]*treeW[[l]][i]*transFFT
        } 
        pBnAll[i,] <- pBn
        pAzeroAll[i,] <- pAzero
        pAsumAll[i,] <- pAsum
      }
    }
    pBnOld <- pBnAll
    pAzeroOld <- pAzeroAll
    pAsumOld <- pAsumAll
  }
  #root
  k <- tree[[1]][1]
  node <- names(tree[[1]][1])
  children <-  1:k
  parentVec <- rep(node, k)
  names(parentVec) <- names(parentTree[[2]])
  parentTree[[2]] <- parentVec 
  if(k>1){
    #parent does not fail:
    pwf <- apply(pBnOld[children,],2,prod)*(1-treePinit[[1]][1])
    #later failure of node but no child triggers failure, the parent does
    plfpTmin1 <- apply(pAzeroOld[children,],2,prod)*transFFT*(1-treePinit[[1]][1])
    #later failure of node:
    pfTmin1 <- apply(pAsumOld[children,],2,prod)*transFFT - plfpTmin1
  }else{
    #parent does not fail:
    pwf <- pBnOld[children,]*(1-treePinit[[1]][1])
    #later failure of node but no child triggers failure, the parent does
    plfpTmin1 <- pAzeroOld[children,]*transFFT*(1-treePinit[[1]][1])
    #later failure of node:
    pfTmin1 <- pAsumOld[children,]*transFFT - plfpTmin1
  }
  matwf[node, ] <- pwf
  matfTmin1[node, ] <-  pfTmin1
  matlfTmin1[node, ] <-  plfpTmin1
  distrFFT <- pwf + pfTmin1
  distr <- Re(fft(distrFFT,inverse = TRUE)[1:(N+1)])/zeroPadDim
  distr <- ifelse(distr < 10^(-15), 0, distr)
  #analysis only reasonable for support
  support <- ifelse(distr > 0, 1, 0)
  #go back in tree from root to leaves:
  divDistr <- ifelse(distr > 0, 1/distr, 1)
  condFailProb <- array(data = 0, dim = c(N, N+1))
  rownames(condFailProb) <- Names
  d <-  Re(fft(pfTmin1,inverse = TRUE)[1:(N+1)])/zeroPadDim
  condFailProb[node,] <- ifelse(d > 10^(-15), d, 0)*divDistr
  #parent -> node
  for(l in 2:levelTree){
    for(i in 1:nbrNodesLevel[l]){
      #belong to level l  
      k <- tree[[l]][i]
      parent <- parentTree[[l]][i]
      node <- names(parent)
      if(k>1){
        #node data 
        pAzeroNode <- matwf[node,]*(1-treeWback[[l]][i]) + matfTmin1[node,]*(1-treeW[[l]][i]) + matlfTmin1[node,]*treeWback[[l]][i]
        pAsumNode <- pAzeroNode + matfTmin1[node,]*(treeW[[l]][i])
        pBnNode <- matwf[node,] + matfTmin1[node,]*(1-treeW[[l]][i])
        #parent given node i
        wfParent <- matwf[parent,]/pBnNode
        lfTmin1Parent <- matlfTmin1[parent,]/pAzeroNode
        fTmin1Parent <- (matfTmin1[parent,] +  matlfTmin1[parent,])/pAsumNode - lfTmin1Parent
        #parent:
        pBn <- wfParent + fTmin1Parent*(1-treeWback[[l]][i])
        pAzero <- wfParent*(1-treeW[[l]][i]) + fTmin1Parent*(1-treeWback[[l]][i]) + lfTmin1Parent*treeW[[l]][i]
        pAsum <-   pAzero + fTmin1Parent*(treeWback[[l]][i]) 
        # make node a root
        matwf[node,] <- matwf[node,]*pBn
        matfTmin1[node,] <-  (matfTmin1[node,] + matlfTmin1[node,])*pAsum 
        matlfTmin1[node,] <- matlfTmin1[node,]*pAzero 
        matfTmin1[node,] <-  matfTmin1[node,] - matlfTmin1[node,] 
        d <- Re(fft(matwf[node,],inverse = TRUE)[1:(N+1)])/zeroPadDim
        condFailProb[node,] <- 1-ifelse(d > 10^(-16), d, 0)*divDistr
      }else{
        #Node
        pBnNode <- (1-treePinit[[l]][i]) + treePinit[[l]][i]*(1-treeW[[l]][i])*transFFT
        pAzeroNode <- (1-treePinit[[l]][i])*(1-treeWback[[l]][i]) + transFFT*((1-treePinit[[l]][i])*treeWback[[l]][i] + treePinit[[l]][i]*(1-treeW[[l]][i]))
        pAsumNode <-   pAzeroNode + treePinit[[l]][i]*treeW[[l]][i]*transFFT
        #parent given node i
        lfTmin1Parent <- matlfTmin1[parent,]/pAzeroNode
        fTmin1Parent <- (matfTmin1[parent,] +  matlfTmin1[parent,])/pAsumNode - lfTmin1Parent
        wfParent <- matwf[parent,]/pBnNode 
        #parent
        pBn <- wfParent + fTmin1Parent*(1-treeWback[[l]][i])
        # make node a root
        pwf <- pBn*(1-treePinit[[l]][i])
        d <- Re(fft(pwf,inverse = TRUE)[1:(N+1)])/zeroPadDim
        condFailProb[node,] <- 1-ifelse(d > 10^(-16), d, 0)*divDistr
      }
    }
  }
  #rounding because of numerical inprecision
  condFailProb <- apply(condFailProb,c(1,2),function(x){min(max(x,0),1)})
  return(list(distr,condFailProb, support))   
}
#histogram of simulation data:
collectDistr <- function(data,N){
  distr <- hist(data,right = TRUE, breaks = seq(-0.5/N,1+0.5/N,1/N))$counts/length(data)
  return(distr)
}
#sample runs times from the cascade size distribution for the ICM and network g
simulateICM <- function(g, runs){
  pInit <- V(g)$pInit
  wMat <- get.adjacency(g,sparse = FALSE,attr = "weight")
  N <- length(pInit)
  nodes <- 1:N
  cs <- rep(0,runs)
  for(r in 1:runs){
    stat <- runif(N, min = 0, max = 1)
    failed <- nodes[stat <= pInit]
    count <- 0
    nbrFailed <- length(failed)
    while(count < nbrFailed){
      count <- count+1
      i <- failed[count]
      nb <- setdiff(nodes[wMat[i,]>0], failed)
      nn <- length(nb)
      if(nn>0){
        stat <- runif(nn, min = 0, max = 1)
        new <- nb[stat <= wMat[i,nb]]
        nbrFailed <- nbrFailed + length(new)
        failed <- c(failed, new)
      }
    }
    cs[r] <- nbrFailed
  }
  distr <- hist(cs,0:(N+1)-0.5,plot = FALSE)$counts
  distr <- distr/runs
  return(distr)
}

