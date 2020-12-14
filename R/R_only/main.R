library(igraph) 
#library(plotrix)
#library(optrees)
#library(data.table)
source("cascadeSizeDistr.R")
source("randomGraphs.R")
set.seed(42)
########################################################
###### Comparison of SDP and TDA with simulations  #####
########################################################
#Tree:
#create networks with different cascade model parameters and compute the corresponding cascade size distributions
#### Trees ####
#random tree consisting of 10 nodes:
N <- 10
g <- randTree(N)
g <- as.directed(g)
plot(g)
#define exemplary cascade model
#initial activation probabilities of nodes:
V(g)$pInit <- rep(0.1,N)
E(g)$weight <- rep(0.2, length(E(g)))

#compute cascade size distribution:
#Faster implementation for trees:
d <- treeICM(g) 
#Function for general networks works as well:
d <- ICMmem(g)[[1]] 
plot((0:N)/N, d, xlab=expression(rho), ylab="P")

#Comparison with sampling:
#number of samples:
m <- 100000
dapprox <- simulateICM(g,m)
plot((0:N)/N, d, xlab=expression(rho), ylab="P", type = "o")
points((0:N)/N, dapprox, col="red")

#### General networks #####
#network drawn from configuration model with scale free degree distribution
N <- 10 #number of nodes
z <- 2 #target average degree
g <- randConfig(10, z)
N <- length(V(g))
g <- as.directed(g)
V(g)$pInit <- rep(0.1,N)
E(g)$weight <- rep(0.2, length(E(g)))
#Plot network:
plot(g)
#Cascade size distribution
#approximation by TDA:
d <- ICMmem(g)[[1]]
m <- 100000
#sampling:
dapprox <- simulateICM(g,m)
#plot:
plot((0:N)/N, d, xlab=expression(rho), ylab="P", type = "o")
points((0:N)/N, dapprox, col="red")

############################################################
########### Conditional activation probabilities ###########
############################################################
out <- condFailProb_ICM(g)
support <- (0:N)/N
#cascade size distribution:
out[[1]]
#conditional activation probabilities:
out[[2]]
#activation probability of node 1 conditional on the cascade size rho = support[2] = 1/N:
out[[2]][1,2] 
#Sometimes, the conditional probabilities are so small that they cannot be computed with numerical accuracy
#out[[3]] indicated with a 1, which parts of the support are provided as output
support <- support[out[[3]]==1]

#Let us compare the conditional activation probabilities of node 1 and 2:
plot(support, out[[2]][1,] , xlab=expression(rho), ylab=expression(paste(plain(P),"( s=1 | ", rho, " )")), type = "o")
points(support, out[[2]][2,] , xlab=expression(rho), ylab=expression(paste(plain(P),"( s=1 | ", rho, " )")), type = "o", col="red")

