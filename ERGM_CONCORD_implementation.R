# This is the implementation for course assignments of STA 650


library(igraph)
library(statnet)

# clean the workspace
rm(list=ls())


## Part I: Analyze the given Network File data

# load data, named as g.sim
load("HW3_network.RData")
ls() # show datas
g.sim
plot(g.sim, main="g.sim network")

# a) number of edges, 2-star, 3-star, and triangles
g.sim.sum <- summary(g.sim ~ edges + kstar(2:3) + triangles) ####

# c) model with ergm
g.sim.01 <- ergm(g.sim ~ edges) # fit with edges only model
summary(g.sim.01)
g.sim.01.sim <- simulate(g.sim.01, nsim=1000) # simulate 1000 graphs
class(g.sim.01.sim)
sum.list <- summary(g.sim.01.sim~edges + kstar(2:3) + triangles) # view summary
# plot the simulation results
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0),mfrow=c(2,2))
hist(sum.list[,2], xlab = "2-star number", main = "simulated 2-star histogram")
abline(v = g.sim.sum[2],col="red") # add reference lien
hist(sum.list[,3], xlab = "3-star number", main = "simulated 3-star histogram")
abline(v = g.sim.sum[3],col="red")
hist(sum.list[,4], xlab = "triangle number", main = "simulated triangles histogram")
abline(v = g.sim.sum[4],col="red")
# recorde AIC and BIC from summary of fitting
AIC(g.sim.01)
BIC(g.sim.01)

# e) fit with extra statistics
g.sim.02 <- ergm(g.sim ~ edges + kstar(2)) 
summary(g.sim.02) # best with smaller AIC and BIC

# f) check the best fit from e)
g.sim.02 <- ergm(g.sim ~ edges + kstar(2)) 
summary(g.sim.02) # best with smaller AIC and BIC
g.sim.02.sim <- simulate(g.sim.02, nsim=1000) # simulate 1000 graphs
sum.list.02 <- summary(g.sim.02.sim~edges + kstar(2:3) + triangles) # view summary
# plot the simulation results
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0),mfrow=c(2,2))
hist(sum.list.02[,2], xlab = "2-star number", main = "simulated 2-star histogram")
abline(v = g.sim.sum[2],col="red") # add reference lien
hist(sum.list.02[,3], xlab = "3-star number", main = "simulated 3-star histogram")
abline(v = g.sim.sum[3],col="red")
hist(sum.list.02[,4], xlab = "triangle number", main = "simulated triangles histogram")
abline(v = g.sim.sum[4],col="red")
# recorde AIC and BIC from summary of fitting
AIC(g.sim.02)
BIC(g.sim.02)

# g). fit with other statistics, expanded for e)
g.sim.03 <- ergm(g.sim ~ edges + kstar(2:3))
AIC(g.sim.03)
BIC(g.sim.03)
summary(g.sim.03)

g.sim.04 <- ergm(g.sim ~ edges + kstar(2:3) + triangles)
summary(g.sim.04)
AIC(g.sim.04)
BIC(g.sim.04)





### Part 2: ERGMs in action. Simulate 1000 networks of 100 or 1000 nodes


g.use2 <- network(100,density=0.1,directed=FALSE,seed=1)
g.sim2 <- simulate(~edges+kstar(2),
                   nsim=1000, coef=c(-1,1),
                   basis=g.use2,
                   control=control.simulate(
                     MCMC.burnin=1000,
                     MCMC.interval=100),seed=2)
class(g.sim2)
#length(g.sim2)
g.sum2.sum <- summary(g.sim2 ~ edges + kstar(2))
#g.sum2.sum # 1000 x 2 list
hist(g.sum2.sum[,1], xlab = "edges", main = "simulated edges histogram")
hist(g.sum2.sum[,2], xlab = "2-star number", main = "simulated 2-star histogram")





### Part 3: Structural equivalence: implement the CONCOR algorithm and apply it to the Karate club network.



library(igraphdata) # for the karte data

# CONCOR implementation

# return two subsets from a socia_matrix
graphDivide <- function(socia_matrix)
{
  # need to check base case here
  tol <- 0.0001
  if( (nrow(socia_matrix) == 1 && ncol(socia_matrix) == 1) || abs(prod(socia_matrix) - 1) < tol ) {
    return(list(group1 = rownames(socia_matrix), group2 = NULL))
  }
  
  cor_matrix <- round(cor(socia_matrix, use = "pairwise"), 3) # return correlation matrix of sociamatrix
  while(abs(prod(cor_matrix) - 1) > tol) {
    cor_matrix <- round(cor(cor_matrix, use = "pairwise"), 3)
  }
  
  #permute to find group1 and group2
  g1 <- names(which(cor_matrix[, 1] == 1))
  g2 <- names(which(cor_matrix[, 1] == -1))
  return(list(group1 = g1, group2 = g2))
}

# divide groups with depth times
ConcorIntoSubGraphs <- function(socia_matrix, curtGroups, Depth) # curtGroups is a list of list
{
  curtLayerGroups <- curtGroups
  for(i in 1:Depth) {
    nextLayerGroups <- list()
    for( j in 1:length(curtLayerGroups)) {
      curtGroup <- curtLayerGroups[[j]]
      index_row <- rownames(socia_matrix) %in% curtGroup
      index_col <- colnames(socia_matrix) %in% curtGroup
      curt_sociaMatrix <- subset(socia_matrix, subset = index_row, select = index_col)
      gd <- graphDivide(curt_sociaMatrix)
      if(!is.null(gd$group1)) {
        nextLayerGroups[[length(nextLayerGroups) + 1]] <- gd$group1
      }
      if(!is.null(gd$group2)) {
        nextLayerGroups[[length(nextLayerGroups) + 1]] <- gd$group2
      }
    }
    curtLayerGroups <- nextLayerGroups
    #print(nextLayerGroups)
  }
  return(curtLayerGroups)
}

# input: adj.matrix, depth of CONCORD algo
# output: a colored graph with diff group with diff color
concordWithAdjmatrix <- function(adjMatrix, depth) {
  vertices <- rownames(adj.matrix)
  groupsList <- list(vertices) # list of vertices list
  groupsList <- ConcorIntoSubGraphs(adj.matrix, groupsList, depth)
  
  #print(groupsList)
  v_groupnum <- findIndex(groupsList, names(V(karate)))
  num_groups <- length(groupsList)
  pal <- rainbow(num_groups, alpha=.5)
  par(mar=c(3,3,1,1),mgp=c(1.75,.75,0),mfrow=c(1,1))
  plot(karate, main = "karate original group", edge.arrow.size=.5, vertex.label.cex= 0.9, vertex.label.color="black", vertex.label.dist=1.5)
  plot(karate, main = paste("karate after", depth,"divide", sep = " "), edge.arrow.size=.5, vertex.label.cex= 0.9, vertex.label.color="black", vertex.label.dist=1.5,
       vertex.color=pal[v_groupnum] )

  return(groupsList)
  # return a list of groups, each group containing vertices' names
}

# return the group # the vertex belong to
findIndex <- function(groups, vertices) {
  loc <- rep(0, length(vertices))
  len <- length(groups)
  for(i in 1:length(loc)) {
    for(j in 1:len) {
      if(vertices[i] %in% groups[[j]]) {
        loc[i] = j
      }
    }
  }
  return(loc)
}

# CONCOR program test #
data("karate")
karate

adj.matrix <- as_adjacency_matrix(karate, type = "both", sparse = FALSE)
depth <- 2
gpList <- concordWithAdjmatrix(adj.matrix, depth)
gpList

