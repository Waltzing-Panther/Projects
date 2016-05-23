library(ggplot2)
library(pracma)

wine = read.table("wine.data")

grpSelector = function(repSet, x){
  x = as.numeric(x)
  tmp = t(apply(repSet, 1, function(y){y - x}))
  d = apply(tmp, 1, function(y) Norm(y))
  grp = as.numeric(which(d == min(d)))
  return(grp)
}

chg = function(repSet1,repSet2){
  distVec = apply(repSet1 - repSet2, 1, function(x) Norm(x))
  return(sum(distVec))
}

calcCenter = function(data,k,items){
  cent = matrix(nrow = k, ncol = length(items))
  for(j in 1:k){
    index = which(data$clusNum == j)
    cent[j,] = colMeans(data[index,items])
  }
  return(cent)
}

# data = data set of interset
# items = vector of columns in data set
# k = number of clusters
kmeans2 = function(data,items,k){
  #Initialize Representatives 
  n = dim(data)[1]
  startIndex = sample(1:n,k, replace = FALSE)
  repSet = as.matrix(data[startIndex,items])
  clusNum = rep(0,n)
  data = cbind(data,clusNum)
  
  for(i in 1:10){
    #Create clusters
    grpVec = apply(data[,items], 1, function(x){grpSelector(repSet,x)})  
    data$clusNum = grpVec
    
    #Calculate improved representatives
    repSet2 = calcCenter(data,k,items)
    delta = chg(repSet,repSet2)
    if(delta < .001)
      break
    else
      repSet = repSet2
  }
  return(data)
}

k = 3
items = c(2,3)
data = kmeans2(wine, items, k)
data$clusNum = as.factor(data$clusNum)
graph = ggplot(data, aes(V2,V3, colour = clusNum)) + geom_point()
graph

data2 = wine
data2$clusNum = kmeans(wine[,items], centers = 3)$cluster
data2$clusNum = as.factor(data2$clusNum)
graph = ggplot(data2, aes(V2,V3, colour = clusNum)) + geom_point()
graph

system.time(kmeans2(wine, items, k))
system.time(kmeans(wine[,items], centers = 3))
