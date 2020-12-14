library(igraph)
library(tictoc)

randTree <- function(N){
  g <- graph.empty(n=N, directed=FALSE)
  g <- add_edges(g, c(1, 2))
  for(i in 3:N){
    source <- sample(i-1, 1)
    g <- add_edges(g, c(source, i))
  }
  return(g)
}
#configuration model random graph with scale free degree distribution and target average degree z
helpDegSeq <- function(l, z, gamma, dmax){
  C = sum(((l+1):dmax)^(-gamma+1))/sum(((l+1):dmax)^(-gamma))
  out = (z - C)/(0.5*l*(l+1)-l*C)
  return(out)
}
randConfig <- function(N, z){
  #generate configuration model networks of different sizes with prescribed degree sequence and small average degree
  #exponent
  gamma <- 3#3
  dmax <- N
  p <- 1:dmax
  p <- p^(-gamma)
  y = sapply(1:(N-1), function(x){helpDegSeq(x, z, gamma, dmax)})
  l <- which(y*N >= 1, arr.ind = TRUE)[1]
  p[1:l] <- y[l]
  C = (1- l*y[l])/sum(((l+1):dmax)^(-gamma))
  p[(l+1):dmax] <- C*((l+1):dmax)^(-gamma)
  q <-round(p*N*(p*N >= 1 ))
  h <- N-sum(q) 
  dmaxEffOrig <- which.min(q)-1 
  if(h>0){
    if(ceiling(h/2)-1 <= dmax){
      q[(dmaxEffOrig-floor(h/2)):(dmaxEffOrig+ceiling(h/2)-1)] <- 1 + q[(dmaxEffOrig-floor(h/2)):(dmaxEffOrig+ceiling(h/2)-1)]
    }else{
      q[(dmaxEffOrig-floor(h/2)):dmax] <- 1 + q[(dmaxEffOrig-floor(h/2)):dmax]
      h <- N-sum(q) 
      q[1:h] <- q[1:h] + 1
    }
  }
  dmaxEff <- which.min(q) -1  
  dSeq <- array(data = 0, dim = c(1,N))
  n <- 0
  for(i in 1:dmaxEff){
    dSeq[(n+1):(n+q[i])] <- i
    n <- n + q[i]
  }
  t <-  try(gg<-degree.sequence.game(dSeq,method = "simple.no.multiple"))
  if("try-error" %in% class(t)){  
    dSeq[1] <- dSeq[1]+1 
    t2 <- try(gg <- degree.sequence.game(dSeq,method = "simple.no.multiple"))
    if("try-error" %in% class(t2)){
      dSeq[1] <- dSeq[1]+1 
      t3 <- try(gg <- degree.sequence.game(dSeq,method = "simple.no.multiple"))
      if("try-error" %in% class(t3)){
        gg <- degree.sequence.game(dSeq,method = "simple.no.multiple")
      }
    }
  }
  c <- components(gg)
  g <- induced_subgraph(gg,V(gg)[c$membership==which(c$csize == max(c$csize), arr.ind = TRUE)])
  return(g)
}
#random network with clustering as referenced in paper
randClust <- function(N, d, o){
  g <- graph.empty(n=N, directed=FALSE)
  g <- add_edges(g, c(1, 2))
  for(i in 3:N){
    v <- i
    deg <- degree(g)[1:(i-1)]
    potNeighb <- 1:(i-1)
    for(j in 1:min(i-1,d)){
      degLoc <- deg[potNeighb]
      p <- degLoc/sum(degLoc)
      if(length(potNeighb) > 1){
        u <- sample(potNeighb,1,prob = p)
      }else{
        u <- potNeighb
      }
      g <- add_edges(g, c(v, u))
      potNeighb <- setdiff(potNeighb, u)
    }
    for(j in 1:o){
      nb <- as_ids(neighbors(g, v))
      if(length(nb) > 1){
        u <- sample(nb,1)
      }else{
        u <- nb[1]
      }
      nbRest <- setdiff(nb,u)
      nbrNb <- length(nbRest) 
      if(nbrNb > 0){
        if(nbrNb > 1){
          w <- sample(nbRest, 1)  
        }else if(nbrNb == 1){
          w <- nbRest[1]
        }
        g <- add_edges(g, c(u, w))
      }
    }
  }
  g <- simplify(g)
  c <- components(g)
  g <- induced_subgraph(g,V(g)[c$membership==which(c$csize == max(c$csize), arr.ind = TRUE)])
  return(g)
}
