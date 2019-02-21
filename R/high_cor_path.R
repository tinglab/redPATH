require(mclust)
require(combinat)
require(MASS)
require(dplyr)
#library(fpc)
#library(ICSNP)
#library(flexmix)
#library(mvtnorm)



#source("ai_trial_new.R")

# FORCE HMT
permutations <- function(n){
  if(n==1){
    return(matrix(1))
  } else {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}

mat2lst <- function(m) {
  m <- t(m)
  dm <- as.data.frame(m)
  as.list(dm)
}


path_length <-function(matrix, order){
  l1 <- sum(matrix[cbind(order[-nrow(matrix)], order[-1])])
  
  return(l1)
}

get_permn <- function(n){
  p1 = permutations(n)
  output = mat2lst(p1)
  return(output)
}
# For Calculation of base path:
#' Hamiltonian path - by brute force
#' 
#' This function calculates the hamiltonian path by brute force. It's only scalable up to N = 10.
#'@param distance_matrix A N by N distance matrix (KL distance in this case). Works on any other asymmetrical / symmetrical distance matrices
#'@keywords hamiltonian_path_force
#'@export
#'@section redPATH pseudotime:
#'@examples
#'hamiltonian_path_as(kl_cpp(data))
#'
hamiltonian_path_force <- function(distance_matrix){
  tmp <- c()
  matrix <- as.matrix(distance_matrix)
  # enumerate all possible paths
  paths <- get_permn(nrow(matrix))
  path_total_lengths <- sapply(paths, function(p) path_length(matrix, p))

  hamilton_path <- paths[[which.min(path_total_lengths)]]
  
  min_dist <- min(path_total_lengths)
  
  return(c(min_dist, hamilton_path))
}


corr_dist <- function(x, y = NULL, method = NULL, use = "everything")
{
  x <- t(as.matrix(x))
  if(!is.null(y))
    y <- t(as.matrix(y))
  1 - (stats::cor(x, y) + 1) / 2
}


normalizeResult2 <- function(clust_num_EM, result, tourlength, ordIndex)
{
  new_result <- result
  new_result <- (new_result) / clust_num_EM  # *tourlength

  return(new_result)
  
}

unnormalize_result <- function(clust_num_EM, result, tourlength, ordIndex)
{
  new_result <- result
  new_result <- (new_result * clust_num_EM) - 0.5   # *tourlength

  return(new_result)
  
}

high_cor_path <- function(clust_num_EM, cell_gene_exprs, stageIdxSmp, HMT_iterations, TSPMethod = "arbitrary_insertion", clustMethod = "GMM", rankType = "class", dist_type = "MHL", debug = FALSE, validate = F)
{
  
  if (clust_num_EM >= dim(cell_gene_exprs)[1])
  {
    cluster_result <- c(1: dim(cell_gene_exprs)[1])
  }
  else if (clustMethod == "Corr")
  {
    dist_mat <- corr_dist(cell_gene_exprs)
    scaled_distance <- .Call(stats:::C_DoubleCentre, dist_mat^2)
    fit_hclust <- hclust(as.dist(ED2(scaled_distance)))
    cluster_result <- as.integer(cutree(fit_hclust, clust_num_EM))
    
  }
  else if (clustMethod == "GMM")
  {
    require(mclust)
    fit_EM = Mclust(cell_gene_exprs, G = clust_num_EM, verbose = F)
    cluster_result <- fit_EM$classification
  }
  
  
  #message("Finished Cluster")
  
  obj_mean_result <- cluster_result
  obj_mean_data <- cell_gene_exprs
  cls_means <- c()
  

  for(cluster_id in c(1:clust_num_EM))
  {
    cluster_i = which(obj_mean_result == cluster_id)
    
    if(length(cluster_i) == 1)
    {
      current_cls_mean <- obj_mean_data[cluster_i, ]

    }
    else
    {
      current_cls_mean <- apply(obj_mean_data[cluster_i, ], 2, mean)
    }
    cls_means <- cbind(cls_means, current_cls_mean)
    
  }

  #print(dim(cls_means))
  ###################################################################################
  # Correlation distance
  
  if(dist_type == "Corr"){
    distance_obj <- cls_means
    distance_obj <- as.matrix(distance_obj)
    distance_matrix <- as.matrix(correlation_distance(distance_obj))
  }
  
  # ###################################################################################
  # Deprecated:
  # if(dist_type == "BT"){
  #   # BHATTACHARYYA DISTR DISTANCE
  #   distance_obj <- cls_means
  #   distance_obj <- as.matrix(distance_obj)
  #   
  #   #print(dim(distance_obj))
  #   mi_distance <- as.matrix(diag(1, nrow = ncol(distance_obj), ncol = ncol(distance_obj)))
  #   
  #   
  #   n_names <- ncol(distance_obj)
  #   
  #   for(i in 1:(n_names-1)){
  #     for(j in (i+1):n_names){
  #       
  #       distr_i <- fitdistr(distance_obj[,i], "normal")
  #       distr_j <- fitdistr(distance_obj[,j], "normal")
  #       
  #       dist_ij <- bhattacharyya.dist(distr_i$estimate[[1]], distr_j$estimate[[1]], distr_i$estimate[[2]], distr_j$estimate[[2]])
  #       
  #       mi_distance[i, j] <- dist_ij#mi_calc[,1]#t(mi_calc)[1,comp_pos]
  #       mi_distance[j, i] <- dist_ij
  #     }
  #   }
  #   #print(dim(mi_distance))
  #   
  #   distance_matrix <- as.matrix(mi_distance)
  # }
  ###################################################################################
  
  # MAHALANOBIS DISTANCE
  if(dist_type == "MHL"){
    distance_obj <- cls_means
    distance_matrix <- as.matrix(dist(t(distance_obj), method = "minkowski"), upper = T, diag = T)
    
  }
  ###################################################################################
  
  # Asymmetric KL distance implemented
  
  if(dist_type == "KL2"){
    
    distance_obj <- as.matrix(cls_means)
    distance_matrix <- as.matrix(kl_cpp(distance_obj))

  }

  else if(dist_type == "ED"){
    distance_obj <- as.matrix(cls_means)
    distance_matrix <- as.matrix(ED2(cls_means))#as.matrix(dist(t(cls_means), upper = T, diag = T))
    
  }


  #print("finished calc dist")
  
  
  clust_num_HMT <- dim(distance_obj)[2]

  rownames(distance_matrix) <- paste("Clust", 1:clust_num_HMT, sep = "")
  colnames(distance_matrix) <- paste("Clust", 1:clust_num_HMT, sep = "")
  #print(distance_matrix)
  
  tourLengthLst <- c()
  highCorLst <- c()
  tourLength <- 0
  resultLst <- matrix(0, nrow = (HMT_iterations * clust_num_EM), ncol = length(cluster_result))
  norm_resultLst <- matrix(0, nrow = (HMT_iterations * clust_num_EM), ncol = length(cluster_result))
  
  forceLst <- matrix(0, nrow = 1, ncol = length(cluster_result))
  norm_forceLst <- matrix(0, nrow = 1, ncol = length(cluster_result))
  forceLength <- c()
  

  if(TSPMethod == "force"){
    
    tmp <- hamiltonian_path_force(distance_matrix)
    tourLength <- tmp[1]
    path_idx <- tmp[2: length(tmp)]
    path_labels <- paste("C", 1:clust_num_HMT, sep = "")[path_idx]
    
    ordIndex <- path_idx
    
    #print(ordIndex)
    ordResult <- (match(c(1:clust_num_EM), ordIndex))[cluster_result]
    #print(distance_matrix)
    normResult <- normalizeResult2(clust_num_EM, ordResult, tourLength, ordIndex)
    forceLength <- c(forceLength, tourLength)
    
    forceLst[1, ] <- ordResult
    norm_forceLst[1, ] <- normResult
    if(validate == T)
    {
      save(cluster_result, cls_means, ordIndex, file = paste0("hamiltonian_path_force_", clust_num_EM, ".RData"))
    } 
    
  }else if(TSPMethod == "arbitrary_insertion"){
    total_iterations <-  HMT_iterations*clust_num_EM
    if(total_iterations > 600){
      total_iterations <- 300
    }
    #for(t in 1:(HMT_iterations*clust_num_EM)){
    for(t in c(1:total_iterations)){

      
      t <- t

      path <- hamiltonian_path_as(distance_matrix)
      tourLength <- path_length(distance_matrix, path)

      
      ordIndex <- path
      
      #print(ordIndex)
      ordResult <- (match(c(1:clust_num_EM), ordIndex))[cluster_result]
      #print(distance_matrix)
      normResult <- normalizeResult2(clust_num_EM, ordResult, tourLength, ordIndex)
      tourLengthLst <- c(tourLengthLst, tourLength)
      #print(tourLengthLst)
      resultLst[t, ] <- ordResult
      norm_resultLst[t, ] <- normResult
    }
    hist(tourLengthLst)
    minLengthID <- which(tourLengthLst == min(tourLengthLst))[1]
  }
  
  
  if(TSPMethod == "force"){
    bestTour <- list(highCor = NULL, norm_result = norm_forceLst[1,], result = forceLst[1, ], tourL = forceLength, dist = distance_matrix * clust_num_EM / max(distance_matrix))

  }else{
    bestTour <- list(highCor = NULL, norm_result = norm_resultLst[minLengthID, ], result = resultLst[minLengthID, ], tourL = tourLengthLst, dist = distance_matrix * clust_num_EM / max(distance_matrix))
  }
  
  return(bestTour)

}



