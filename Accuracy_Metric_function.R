# Helper function: Compute the accuracy of the resulted Bayesian network

ACC = function(actual, learning, n){
  # "actual" should be a data frame whose two column are "from" and "to". So does "learning".
  
  actual_arcs = paste(actual[,1], actual[,2])
  learning_arcs = paste(learning[,1], learning[,2])
  
  TP = length(intersect(actual_arcs, learning_arcs))
  TN = n*(n-1) - length(union(actual_arcs, learning_arcs))
  return((TP+TN)/(n*(n-1)))
}