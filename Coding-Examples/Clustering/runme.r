source("EM-Algorithm/EM.r")
source("EM-Algorithm/plotGen.r")
source("K-Means/kmeans.r")

# clustering of old faithful data using EM-Algorithm
data = faithful
clusterCount = 2
xres = 100
yres = 100
seed = 100
iter = 20
drawplot(data, clusterCount, xres, yres, seed, iter)

# clustering of wine data using EM-Algorithm. The variables malic acid and ash are used.
# for a short description of the wine data set see:
# https://archive.ics.uci.edu/ml/datasets/Wine
data = read.csv("wine.data")
data = data[,c(2,3)]
clusterCount = 3
xres = 100
yres = 100
seed = 100
iter = 100
drawplot(data, clusterCount, xres, yres, seed, iter)
  
