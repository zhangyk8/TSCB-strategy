## Functional programming of the two-step clustering-based Bayesian network structrue learning strategy. 
# (Standard Version)
cluster_BNL = function(data, stand = TRUE, method = "average", N = floor(sqrt(2)*ncol(data)/2), algo = gs, ...){
  
  # The function takes a data frame as an input. 
  # "stand" is a logical variable indicating whether a standization should be performed. 
  # "method" is used for the clustering. 
  # "N" tuned the number of groups into which the variables are divided.
  
  require(bnlearn)
  require(cluster)
  require(infotheo)
  
  # Identify the score-based method
  if(names(as.list.function(algo))[2] == "start"){
    data1 = as.data.frame(lapply(data, as.numeric))
    if(stand){
      data1 = apply(data1, MARGIN = 2, function(column){
        column - mean(unique(column))
      })
      dissim = as.dist(matrix(1, nrow = ncol(data1), ncol = ncol(data1)) - cor(data1))
    }else{
      dissim = as.dist(matrix(1, nrow = ncol(data1), ncol = ncol(data1)) - cor(data1))
    }
  }else{
    
    if(all(sapply(data, class) == "factor")){
      dissim = as.dist(-mutinformation(data))
    }else if(all(sapply(data, class) == "numeric")){
      dissim = as.dist(matrix(1, nrow = ncol(data), ncol = ncol(data)) - cor(data))
    }else{
      data1 = as.data.frame(lapply(data, as.numeric))
      data1 = apply(data1, MARGIN = 2, function(column){
        column - mean(unique(column))
      })
      dissim = as.dist(matrix(1, nrow = ncol(data1), ncol = ncol(data1)) - cor(data1))
    }
  }
  
  dissim[is.na(dissim)] = 2
  cluster1 = agnes(dissim, diss = inherits(dissim, "dist"), stand = FALSE, method = method)
  Group = cutree(cluster1, k = N)
  
  whiteArcs = data.frame(from = c(), to = c())
  for(i in 1:N){
    index = (Group == i)
    
    if(sum(index) > 1){
      bn.cluster = algo(data[ ,index])
      whiteArcs = rbind(whiteArcs, bn.cluster$arcs)
    }
  }
  
  if(nrow(whiteArcs) == 0){
    whiteArcs = NULL
    bn.cl1 = algo(data, whitelist = whiteArcs)
    return(bn.cl1)
  }
  
  if(names(as.list.function(algo))[2] == "start"){
    init_net = empty.graph(names(data))
    arcs(init_net, check.cycles = FALSE) = as.matrix(whiteArcs)
    bn.cl1 = algo(data, start = init_net)
  }else{
    bn.cl1 = algo(data, whitelist = whiteArcs)
  }
  
  return(bn.cl1)
}



## Fast Version

cluster_BNL = function(data, stand = TRUE, method = "average", N = floor(sqrt(2)*ncol(data)/2), algo = gs, ...){
  
  require(bnlearn)
  require(cluster)
  
  if(sum(sapply(data, class) != "factor") > 0) stop("The class of the data frame should be factor!")
  
  data1 = as.data.frame(lapply(data, as.numeric))
  if(stand & any(sapply(data, class) == "factor")){
    data1 = apply(data1, MARGIN = 2, function(column){
      column - mean(unique(column))
    })
    dissim = as.dist(matrix(1, nrow = ncol(data1), ncol = ncol(data1)) - cor(data1))
  }else{
    dissim = as.dist(matrix(1, nrow = ncol(data1), ncol = ncol(data1)) - cor(data1))
  }
  
  dissim[is.na(dissim)] = 2
  cluster1 = agnes(dissim, diss = inherits(dissim, "dist"), stand = FALSE, method = method)
  Group = cutree(cluster1, k = N)
  
  whiteArcs = data.frame(from = c(), to = c())
  for(i in 1:N){
    index = (Group == i)
    
    if(sum(index) > 1){
      bn.cluster = algo(data[ ,index])
      whiteArcs = rbind(whiteArcs, bn.cluster$arcs)
    }
  }
  
  if(nrow(whiteArcs) == 0){
    whiteArcs = NULL
    bn.cl1 = algo(data, whitelist = whiteArcs)
    return(bn.cl1)
  }
  
  if(names(as.list.function(algo))[2] == "start"){
    init_net = empty.graph(names(data))
    arcs(init_net, check.cycles = FALSE) = as.matrix(whiteArcs)
    bn.cl1 = algo(data, start = init_net, whitelist = whiteArcs)
  }else{
    bn.cl1 = algo(data, whitelist = whiteArcs)
  }
  
  return(bn.cl1)
}