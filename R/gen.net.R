# Function used to generate network Gx.
gen.net = function(p, nc = 10, pc = 0.4, pcc = 0.13, psn = 0.04, normalize = FALSE, plot.graph = FALSE){
  ## nc is the number of predictors in each cluster
  ## pc is the probability of connection between variables belonging to the same cluster
  ## pcc is the probability of connection between signal-carriers or non-signal-carriers belonging to different clusters
  ## psn is the probability of connection between a signal-carrier and a non-sginal-carrier.
  
  gx.mat = matrix(0,p,p)    ## Adjacency matrix
  idx.cluster = rep(1:(p/nc), each = nc)
  
  idx.s = c(1,2) # belong to signal carriers or non-signal carriers
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      if(idx.cluster[i]==idx.cluster[j]){
        prob = pc
      }else{
        if(sum(idx.cluster[c(i,j)]%in%idx.s)==1 ){
          ## connection between a signal carrier and a non-signal carrier
          prob = psn
        }else{
          prob = pcc
        }
      }
      gx.mat[i,j] = gx.mat[j,i] = rbinom(1,1,prob)
    }
  }
  colnames(gx.mat) = rownames(gx.mat) = 1:p
  
  if(plot.graph == TRUE){
    require(igraph)
    net = graph_from_adjacency_matrix(gx.mat)
    colrs = c("gray50", "tomato", "gold")
    V(net)$color = colrs[V(net)$media.type]
    plot.igraph(net, edge.arrow.size=0.05, vertex.size = 5, vertex.label.cex = 0.5,vertex.color	=idx.cluster)
    
    # palf = colorRampPalette(c("gold", "dark orange"))
    # heatmap(gx.mat, Rowv = NA, Colv = NA, col = palf(100), scale="none", margins=c(10,10) )
  }
  
  ## 
  A = gx.mat  # adjacency matrix
  A_deg = apply(abs(A),1,sum)
  if(normalize == TRUE){
    L = diag(1, p ) - diag(1/sqrt(A_deg)) %*% A %*% diag(1/sqrt(A_deg)) ## Symmetric normalized Laplacian matrix
  }else{
    L = diag(A_deg) - A  ## Laplacian matrix
  }
  
  return(list( L = L, A = A))
  
}