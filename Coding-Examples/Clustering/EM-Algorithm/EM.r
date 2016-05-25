library(ggplot2)
library(mvtnorm)

# calcMemProb calculates the membership probabilities
# input:
#   mu: matrix containing cluster means
#   sigma: stacked covariance matrices
#   mix: mixing probabilities
#   y: matrix of observed values (observations given by column)
#   n: number of observations 
#   m: number of clusters
# output:
#   M: matrix of membership probabilities for each point
calcMemProb = function(mu,sigma,mix,y,n,m){
  M = matrix(nrow = m, ncol = n)
  for(i in 1:n){
    for(j in 1:m){
      num = dmvnorm(y[,i], mean = mu[,j], sigma = sigma[,,j]) * mix[j]
      den = sapply(1:m, function(x){dmvnorm(y[,i], mean = mu[,x],
                                          sigma = sigma[,,x]) * mix[x]})
      M[j,i] = num/sum(den)
    }
  }
  return(M)
}

# cluster determines which cluster each point belongs to
# inputs:
#   M: matrix of membership probabilities for each point
# output:
#   vector of values from 1 to (number of clusters) that 
#     indicates the cluster that the point belongs to
cluster = function(M){
  return(sapply(1:dim(M)[2],function(x) which.max(M[,x])))
}


# EMest1 uses EM-algorithm to estimate parameters in gaussian
#   mixture model
# input:
#   data: matrix of observations (observations given by row)
#   clusterNum: number of clusters to use
#   delta: precision parameter (e.g delta = .01) ensures 
#      algorithm only exits when difference for each
#      parameter estimate between sucessive iterations 
#      is less that .01
#   maxiter: maximum number of iterations before termination
EMest1 = function(data,clusterNum, delta, maxiter){
  data = unname(data)
  n = dim(data)[1]
  m = clusterNum
  d = dim(data)[2]
  obs = t(data)
  
  #Initialize guess
  mu = matrix(0,nrow = d, ncol = m)
  mu = unname(obs[,sample(1:n,m)])
  
  sigma = array(1,dim=c(d,d,m))
  for(i in 1:m){
    sigma[,,i] = diag(apply(obs,1,var),d)
  }
  mix = rep(1/2, m)
  
  if (maxiter == 0){
    M = calcMemProb(mu,sigma,mix,obs,n,m)
    return(list(mix = mix, mu = mu, sigma = sigma, mem = M))
  }
  
  for(k in 1:maxiter){
    mix_old = mix
    mu_old = mu
    sigma_old = sigma
    
    M = calcMemProb(mu,sigma,mix,obs,n,m)
    
    for(j in 1:m){
      mix[j] = sum(M[j,])/n
    }
    
    for(j in 1:m){
      for(r in 1:d){
        mu[r,j] = sum(M[j,] * obs[r,])
        mu[r,j] = mu[r,j]/sum(M[j,])
      }
    }
    
    for(j in 1:m){
      tmp = matrix(0,d,d)
      for(i in 1:n){
        tmp = tmp + M[j,i]*((obs[,i] - mu[,j]) %*% t((obs[,i] - mu[,j])))
      }
      tmp = tmp/sum(M[j,])
      sigma[,,j] = tmp
    }
    
    m1 = max(abs(mu_old - mu))
    m2 = max(abs(sigma_old - sigma))
    m3 = max(abs(mix_old - mix))
    
    if(m1 < delta && m2 < delta &  m3 < delta){
      print(paste("Number of iterations: ", k))
      return(list(mix = mix, mu = mu, sigma = sigma, mem = M))
    }
  }
  M = calcMemProb(mu,sigma,mix,obs,n,m)
  
  warning("desired convergence not achieved within maximum number of iterations")
  return(list(mix = mix, mu = mu, sigma = sigma, mem = M))
}


# EMest2 uses EM-algorithm to estimate parameters in gaussian
#   mixture model
# input:
#   data: matrix of observations (observations given by row)
#   clusterNum: number of clusters to use
#   maxiter: number of iterations to run before termination
EMest2 = function(data,clusterNum, maxiter){
  data = unname(data)
  n = dim(data)[1]
  m = clusterNum
  d = dim(data)[2]
  obs = t(data)
  
  #Initialize guess
  mu = matrix(0,nrow = d, ncol = m)
  mu = unname(obs[,sample(1:n,m)])
  
  sigma = array(1,dim=c(d,d,m))
  for(i in 1:m){
    sigma[,,i] = diag(apply(obs,1,var),d)
  }
  mix = rep(1/2, m)
  
  if (maxiter == 0){
    M = calcMemProb(mu,sigma,mix,obs,n,m)
    return(list(mix = mix, mu = mu, sigma = sigma, mem = M))
  }
  
  for(k in 1:maxiter){
    M = calcMemProb(mu,sigma,mix,obs,n,m)
    
    for(j in 1:m){
      mix[j] = sum(M[j,])/n
    }
    
    for(j in 1:m){
      for(r in 1:d){
        mu[r,j] = sum(M[j,] * obs[r,])
        mu[r,j] = mu[r,j]/sum(M[j,])
      }
    }
    
    for(j in 1:m){
      tmp = matrix(0,d,d)
      for(i in 1:n){
        tmp = tmp + M[j,i]*((obs[,i] - mu[,j]) %*% t((obs[,i] - mu[,j])))
      }
      tmp = tmp/sum(M[j,])
      sigma[,,j] = tmp
    }
  }
  
  M = calcMemProb(mu,sigma,mix,obs,n,m)
  return(list(mix = mix, mu = mu, sigma = sigma, mem = M))
}

# EMest is a wrapper function for EMest1 and Emest2 which
#   allows us to specify either the delta term or the maxiter term.
EMest = function(data,clusterNum, delta = .01, maxiter = 10, type = 1){
  if(type == 1)
    return(EMest1(data,clusterNum, delta, maxiter))
  else if(type == 2)
    return(EMest2(data,clusterNum, maxiter))
}

