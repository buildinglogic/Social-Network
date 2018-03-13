rm(list=ls())
library(igraph)

## return the Durfee number of a non-increasing degree seq
# algo: use a binary search to find the biggest j, with d_j >= j -1
Durfee <- function(d) { # d is a sorted seq, non-increasing
  start <- 1
  end <- length(d)
  if(end == 1) {return(1)}
  while(start < end - 1) {
    mid <- floor( (start + end) / 2 )
    if(d[mid] >= mid - 1) { start <- mid } 
    else {end <- mid - 1}
  }
  if(d[start + 1] >= start) {start <- start + 1}
  return(start)
}

#### test #####
# d = c(3, 1, 1, 1)
# df <- Durfee(d)

## isgraphical(d): implement Erodos-Gallai algo
# optimize with Durfee number, only validate the first m equations
isgraphical <- function(d) {
  if(sum(d) %% 2 == 1) return(FALSE) # if sum is odd, not graphical
  d_sorted = sort(d, decreasing = TRUE)
  n = length(d)
  if(n == 0) {return(TRUE)}
  if(d_sorted[n] < 0) {return(FALSE)} 
  dfnum = Durfee(d_sorted)
  presum <- rep(0, n)
  minflip <- rep(0, n) # pos of flip of min(k, d_i)
  
  presum[1] = d_sorted[1]
  for(i in 2:n) {
    presum[i] = presum[i-1] + d_sorted[i]
  }
  j <- n
  for(k in 1:n) {
    while(j > 0 && d_sorted[j] < k) {j = j - 1}
    minflip[k] = j # the last d_j >= k
  }
  
  for(k in 1:dfnum) {
    min_k_num = max(c(0, minflip[k] - k)) # terms with k smaller in \sum_k+1 ^n min(k, d_i) in Erodos-Gallai
    postsum = k*min_k_num + presum[n] - presum[k + min_k_num]
    if(presum[k] > k*(k-1) + postsum) {return(FALSE)}
  }
  return(TRUE)
}

#### test ####
# d = c(1, 1, 1, 3)
# g <- isgraphical(d)
# View(d) # sort is not in place sorting

## find the d_i with min positive, otherwise return 0
findMinPositive <- function(d) {
  m = 1;
  for(i in 1:length(d)) {
    if(d[m] == 0 || (d[i] > 0 && d[i] < d[m])) {
      m = i;
    }
  }
  if(d[m] == 0) return(0) # no positive numer
  else return(m)
}

## compute wether the edge is already been assigned
isVisitedEdge <- function(u, v, edges) {
  m <- nrow(edges)
  if(m == 0) return(FALSE) # in case out of bounds
  #if(u > m || v > m) return("error: the vertices are not in the graph")
  for(i in 1:m) {
    if( (edges[i,1] == u && edges[i,2] == v) || (edges[i,2] == u && edges[i,1] == v) ) {
      return(TRUE)
    }
  }
  return(FALSE)
}

## after remove an edge {i, j}, d is still graphical or not
isgraphicalLessEdge <- function(i, j, d) {
  tempd <- d;
  tempd[i] <- tempd[i] - 1;
  tempd[j] <- tempd[j] - 1;
  # View(d) # tempd is another copy, not in place
  if(isgraphical(tempd)) return(TRUE)
  else return(FALSE)
}

## compute candidates list J = (j != i, {i, j} \in E and Oij d is graphical)
findCandidateList <- function(v, d, edges) {
  candidatelist <- numeric();
  numCandidate <- 0
  n <- length(d)
  for(i in 1:n) {
    if(i != v && d[i] > 0 && !isVisitedEdge(v, i, edges) && isgraphicalLessEdge(v, i, d)) {
      numCandidate = numCandidate + 1;
      candidatelist[numCandidate] = i;
    }
  }
  return(candidatelist)
}

# #### test ####
# d = c(2, 2, 0,2)
# edges <- matrix(nrow=0, ncol=2)
# can <- findCandidateList(1,d,edges)
# g <- isgraphical(c(2, 2, 0, 2))

## find one neighboring edge
findOneNeighbor <- function(v, d, edges) {
  candlist = findCandidateList(v, d, edges)
  # print(candlist)
  len = length(candlist)
  if(len == 0) {stop("candlist is empty")}
  candDegreeList <- rep(0, len) # degrees of candidates
  for(i in 1:len) {
    candDegreeList[i] = d[candlist[i]]
  }
  
  cd_len = length(candDegreeList) 
  j <- sample(1:len,1,replace=FALSE,prob=candDegreeList) # pick a cand proportional to its degree
  neighbor <- candlist[j]
  
  newd <- d;
  newd[v] = newd[v] - 1
  newd[neighbor] = newd[neighbor] - 1
  newedges <- edges
  newedges <- rbind(edges,c(v,neighbor))
  
  return(list(new.degrees = newd, new.edges = newedges, prob = candDegreeList[j]/sum(candDegreeList)))
}

# #### test ####
# d = c(1, 1, 0, 2)
# edges <- matrix(nrow=0, ncol=2)
# can <- findCandidateList(1,d,edges)
# result <- findOneNeighbor(1, d, edges)
# View(result$new.edges)


## findAllNeighbors() : find all neighbors for a given vertex
## then assign connections to neighbors to the output edges
findAllNeighbors <- function(v, d, edges) {
  tempD <- d
  tempE <- edges
  p <- 1
   while(tempD[v] > 0) {
     removeOneEdge <- findOneNeighbor(v, tempD, tempE)
     tempD <- removeOneEdge$new.degrees
     tempE <- removeOneEdge$new.edges
     p <- p * removeOneEdge$prob
   }
  return(list(new.degrees = tempD, new.edges = tempE, prob = p))
}

# #### test ####
# d = c(2, 2, 0, 2)
# edges <- matrix(nrow=0, ncol=2)
# can <- findCandidateList(1,d,edges)
# result <- findAllNeighbors(1, d, edges)
# View(result$new.edges)
# View(result$new.degrees)
# View(result$prob)

generateSeqRandomGraph <- function(d)
{
  if (isgraphical(d) == FALSE) return("error: the input is not graphical")
  n <- length(d)
  edges <- matrix(nrow=0, ncol=2)
  degrees <- d
  edgeSeqProb <- 1 # the prob of generate the edge sequence, sigma(Y)
  equivClassSize <- 1 # c(Y)
  while ((v <- findMinPositive(degrees)) > 0) {
    removeOneVertex <- findAllNeighbors(v, degrees, edges)
    degrees <- removeOneVertex$new.degrees
    edges <- removeOneVertex$new.edges
    edgeSeqProb <- edgeSeqProb * removeOneVertex$prob
    equivClassSize <- equivClassSize * factorial(degrees[v])
    # print(degrees)
  }
  return(list(n=n,degrees=d,edges=edges,c_Y=equivClassSize,sigma_Y=edgeSeqProb,weight=1/(equivClassSize*edgeSeqProb)))
}

## return a estimate graph number vector
estimateGraphNum <- function(d, runs) {
  graphNum <- numeric()
  count <- 0
  for(i in 1:runs) {
    result <- generateSeqRandomGraph(d)
    count = count + 1
    graphNum[count] <- result$weight
  }
  return(graphNum)
}

# #### test ####
# d = c(2, 2, 1, 3)
# result <- generateSeqRandomGraph(d)
# edges <- result$edges
# g <- graph_from_edgelist(edges, directed = FALSE)
# par(mar=c(3,3,1,1),mgp=c(1.75,.75,0),mfrow=c(1,2))
# plot(g)

### test generate one random graph with karate data ###
par(mar=c(3,3,1,1),mgp=c(1.75,0.75,0),mfrow=c(1,2))

data(karate)
karate
plot(karate, main = "karate club")

d <- degree(karate)
result <- generateSeqRandomGraph(d)
edges <- result$edges
g <- graph_from_edgelist(edges, directed = FALSE)
plot(g, main = "a random graph from karate club")

### estimate the number of graphs and MC error ####
# let pi be uniform, then E(weight) = number of graphs
runs <- 20
numGraph <- estimateGraphNum(d, runs)
hat_mu <- mean(numGraph) # estimate
mcerror <- var(numGraph) # Mento carlo error
print(paste("estimated num of graphs:", hat_mu, "and Mento Carlo error", mcerror, "with run = ", runs,seq = " "))





