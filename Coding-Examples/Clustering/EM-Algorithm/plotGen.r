library(grid)
library(gridExtra)
library(animation)
library(foreach)
library(parallel)

# drawplot creates a cluster plot of 2-dimensional data. The EM-Algorithm
#   for Gaussian mixture models is used. 
# input:
#   data: a 2 column matrix
#   xres: courseness parameter for x-axis
#   yres: courseness parameter for y-axis
#   clusterCount: number of clusters to be used
#   seed: seed for random generator
#   iter: number of iterations for the EM-Algorithm
drawplot = function(data, clusterCount, xres, yres, seed, iter){
  set.seed(seed)
  rlenX = abs(min(data[,1]) - max(data[,1]))
  rlenY = abs(min(data[,2]) - max(data[,2]))
  marginsX = .05 * rlenX
  marginsY = .05 * rlenY
  
  X=seq(min(data[,1]) - marginsX, max(data[,1]) + marginsX, rlenX/xres)
  Y=seq(min(data[,2]) - marginsY, max(data[,2]) + marginsY, rlenY/yres)
  
  X = rep(X,length(Y))
  X = sort(X)
  
  plotPoints = matrix(nrow = length(X), ncol = clusterCount + 2)
  plotPoints = as.data.frame(plotPoints)

  tmp = matrix(nrow = length(X), ncol = 2)
  tmp[,1] = X
  tmp[,2] = Y
  
  res = EMest(data, clusterCount, maxiter = iter, type = 2)
  mu = res$mu
  sigma = res$sigma
  data$cluster = cluster(res$mem)
  data$cluster = as.factor(data$cluster)
  plotPoints[,1] = X
  plotPoints[,2] = Y
  for(i in 1:clusterCount){
    plotPoints[,i+2] = apply(tmp,1,function(x){ dmvnorm(x, mean = mu[,i], sigma = sigma[,,i])})
  }
  
  colnames(plotPoints)[1:2] = c("x","y")
  colnames(data)[1:2] = c("x","y")
  
  p1 = ggplot() + geom_point(data = data, aes_string("x","y", colour = "cluster"))
  
  for(i in 1:clusterCount){
    p1 = p1 + stat_contour(data = plotPoints,
                           aes_string("x","y", z = colnames(plotPoints)[i+2]),
                           bins = 8) 
  }
  
  p1 = p1 + coord_cartesian(xlim=c(min(data[,1]) - marginsX, max(data[,1]) + marginsX)
                            , ylim=c(min(data[,2]) - marginsY, max(data[,2]) + marginsY))
  p1 = p1 + ggtitle(paste(iter, "Iterations"))
  p1 = p1 + theme(plot.title = element_text(size = 22))
  
  return(p1)
}


# dWrapper creates a gif of drawplot graphs evaluate at the values indicated
#   in iter.
dWrapper = function(data, clusterCount, seed, xres, yres, iter){
  no_cores = detectCores() - 1
  cl = makeCluster(no_cores, type = "FORK")
  
  clusterEvalQ(cl, library(mvtnorm))
  clusterEvalQ(cl, library(ggplot2))
  clusterExport(cl, c("data","clusterCount", "seed", 
                      "xres", "yres", "drawplot",
                      "EMest", "EMest1", "EMest2",
                      "calcMemProb", "cluster"))
  retList = parLapply(cl, iter, function(x) drawplot(data, clusterCount, xres, yres,seed, x))
  stopCluster(cl)
  
  return(retList)
}



