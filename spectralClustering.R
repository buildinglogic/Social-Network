rm(list = ls())

library(igraph)
library(RSpectra)
library(expm) # sqrtm of matrix
library(matrixcalc) # matrix operation
library(stats) # for k-means

# first
# generate graph data with buildin stochastic blockmodel in igraph
# the number of groups k=  5,  
# number  of  nodes  within  each  group s=  10,  
# probability  of  an  edge  betweengroups r= 0.1, 
# probability of an edge within a groupp+r= 0.4.


s_pool <- c(10, 20, 30, 40, 50)

r <- 0.1
p <- 0.3
k <- 5
s <- 10

misclass(k, s, r, p)

result <- c()
for(s in s_pool) {
  aver_misclass <- 0
  num_repeat <- 1
  for(i in 1:num_repeat) {
    aver_misclass = aver_misclass + (1/num_repeat) * misclass(k, s, r, p)
  }
  result <- c(result, aver_misclass)
}
n_pool <- s_pool * 5
plot(result ~ n_pool, type = "o", ylab = "number of misclass vertices", xlab = "n", main = "number of misclassified as n grow")

misclass <- function(k, s, r, p) {

  n <- k * s
  pref.matrix <- matrix(data = rep(r, k*k), nrow = k)
  diag(pref.matrix) <- p
  g <- sample_sbm(n, pref.matrix, block.sizes = rep(s,k), directed = FALSE, loops = FALSE)
  
  # num_groups <- k
  # pal <- rainbow(num_groups, alpha=.5)
  # par(mar=c(3,3,1,1),mgp=c(1.75,.75,0),mfrow=c(1,1))
  # plot(g, main = "stochastic blockmodel graph", edge.arrow.size=.5, vertex.label.cex= 0.9, vertex.label.color="black", vertex.label.dist=1.5,
  #      vertex.color=pal)
  
  
  # implement the spectral clustering algorithm
  A <- as_adjacency_matrix(g, type = "both", sparse = FALSE)
  D <- diag(apply(A, 1, sum))
  L <- matrix.inverse(sqrtm(D)) %*% A %*% matrix.inverse(sqrtm(D)) # return 0 matrix???
  
  # # computes M^power
  # "%^%" <- function(M, power)
  #   with(eigen(M), vectors %*% (values^power * solve(vectors)))
  # 
  # L <- (D %^% (-1/2)) %*% A %*% (D %^% (-1/2))  # normalized Laplacian
  
  # use RSpectra library
  # return the k eigenvalues/vectors with largest absolute values
  res = eigs(L, k, which = "LM")  # "LM" is the default
  X <- res$vectors
  
  # evl <- eigen(L, symmetric =TRUE)
  # head(evl$values, k)
  # # evl$vectors # columns are eigenvectors
  # # take the k vectors to the k smallest eigenvalues
  # X  <- evl$vectors[,(ncol(evl$vectors)-k+1):ncol(evl$vectors)]
  
  # k-mean to clustering
  km <- kmeans(X, centers=k)
  sigma <- km$cluster # which group belongs to for each vertex
  # plot(g, main = "spectral clustering modeling", edge.arrow.size=.5, vertex.label.cex= 0.9, vertex.label.color="black", vertex.label.dist=1.5,
  # vertex.color=pal[sigma])


  # clustering from population
  Z <- matrix(0, n, k)
  Z[1:s, 1] = 1
  Z[(s+1):(2*s), 2] = 1
  Z[(2*s+1):(3*s), 3] = 1
  Z[(3*s+1):(4*s), 4] = 1
  Z[(4*s+1):(5*s), 5] = 1
  
  E_AZ <- Z %*% pref.matrix %*% t(Z)
  D_wave <- diag(apply(E_AZ, 1, sum))
  L_wave <- matrix.inverse(sqrtm(D_wave)) %*% E_AZ %*% matrix.inverse(sqrtm(D_wave)) # return 0 matrix???
  
  res_wave = eigs(L_wave, k, which = "LM")  # "LM" is the default
  X_wave <- res_wave$vectors
  km_wave <- kmeans(X_wave, centers=k)
  
  
  # centroids
  Cbig <- km$centers[km$cluster,]
  Zmu <- km_wave$centers[km_wave$cluster,]
  s_vd <- svd(t(X_wave) %*% X)
  Orotate <- s_vd$u %*% t(s_vd$v)
  
  check_good <- TRUE
  while(check_good){
    km <- kmeans(X,k)
    Cbig <- km$centers[km$cluster,]
    check_good <- (sum((X-Cbig)^2)> sum((X-Zmu%*%Orotate)^2))
  }
  
  # calculate misclustering rate
  P <- max(diag(t(Z) %*% Z))
  M <- c()
  for(i in 1:n) {
    if( sum((Cbig[i,] - Zmu[i,] %*% Orotate)^2) >= 1/(2 * P) ) {
      M <- c(M, i)
    }
  }
  
  numbermisclasified <- length(M)
  misclusteringRate <- length(M) / n
  return(numbermisclasified)
}



