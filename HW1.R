library(igraph)
library(igraphdata)

#HW1: find adjacency matrix, edge list, degree distribution
g1 <- graph( edges=c(1,2, 2,3, 3,1, 4,1, 1,5, 4,5, 5,4, 5,3), n=5, directed=T ) # an undirected graph with edges
plot(g1)

#adjacency matrix
g1[] # sparse adjacency matrix
get.adjacency(g1,sparse=FALSE) #adjacency matrix, not sparse

#edge list
V(g1) # vertices
E(g1) # edge list
as_edgelist(g1) # edge list in table

#degree distribution
# Node degrees
# 'degree' has a mode of 'in' for in-degree, 'out' for out-degree,
# and 'all' or 'total' for total degree. 
deg <- degree(g1, mode="all")
hist(deg, breaks=1:vcount(g1)-1, main="Histogram of node degree")


#### p3
# plot from literal
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0),mfrow=c(2,2))
plot(a <- graph_from_literal(a--b, b--c, c--d, d--e), main = "a") # graph a
plot(b <- graph_from_literal(a--b, b--c, c--d, c--e), main = "b") # graph b
plot(c <- graph_from_literal(a--b, a--c, a--d, a--e), main = "c") # graph c
plot(d <- graph_from_literal(a--b, b--c, c--d, d--e, e--a), main = "d") # graph d
dev.off
# in graph "a", vertex c is the most central one, since it is in the middle
# in graph "b", vertex c, it not only has most degree, but also in the middle
# in graph "c", vertext a
# in grah "d", equal central

#Report in one table per graph the degree, closeness and betweenness centralities
#### centrality and centralization

# for loop

for(g in c("a","b", "c", "d")) {
  # Degree (number of ties)
  print(paste("graph : ", g))
  
  print("degree centrality")
  #print(degree(get(g), mode="all"))
  print(centr_degree(get(g), mode="all", normalized=T))
  # when norm is true, the centr will devide the theoreteical_max (n)(n-1), which
  # should be (n-1)(n-2)
  
  print("closeness centrality")
  # Closeness (centrality based on distance to others in the graph)
  # Inverse of the node's average geodesic distance to others in the network
  #print(closeness(get(g), mode="all", weights=NA) )
  print(centr_clo(get(g), mode="all", normalized=T) )
  
  print("betweenness centrality")
  # Betweenness (centrality based on a broker position connecting others)
  # (Number of geodesics that pass through the node or the edge)
  #print(betweenness(get(g), directed=F, weights=NA))
  #print(edge_betweenness(get(g), directed=F, weights=NA))
  print(centr_betw(get(g), directed=F, normalized=T))
}



