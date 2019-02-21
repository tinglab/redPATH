# Modified arbitrary insertion for hamiltonian path

calc_insertion_cost_c <- function(random_dist_mat, current_order, vk){
  
  length <-length(current_order)
  costs <- integer((length - 1L))
  n <- nrow(random_dist_mat)
  vk <- vk - 1
  #t1 <- Sys.time()
  for(i in c(1:(length - 1L))){
    d_ik <- (random_dist_mat[(current_order[i]  + n * vk)])# + random_dist_mat[vk, current_order[i]]) / 2 
    #print(current_order[i+1L])
    d_kj <- (random_dist_mat[(vk + n * (current_order[i+1L] ))])# + random_dist_mat[current_order[i+1L], vk]) / 2
    d_ij <- (random_dist_mat[((current_order[i]) + n* (current_order[i+1L]))])# + random_dist_mat[current_order[i], current_order[i+1L]]) / 2
    
    #print(d_ik)
    #print(d_kj)
    #print(d_ij)
    costs[i] <- d_ik + d_kj - d_ij
  }
  #print(current_order)
  #print(vk)
  costs <- c(random_dist_mat[(vk*n + current_order[1])], costs, random_dist_mat[( current_order[length] +  n*vk)])
  print(costs)
  #print(Sys.time() - t1)
  return(costs)
  
}
calc_insertion_cost <- function(random_dist_mat, current_order, vk){
  
  length <-length(current_order)
  costs <- integer((length - 1L))
  t1 <- Sys.time()
  for(i in c(1:(length - 1L))){
    d_ik <- (random_dist_mat[current_order[i], vk])# + random_dist_mat[vk, current_order[i]]) / 2 
    #print(current_order[i+1L])
    d_kj <- (random_dist_mat[vk, current_order[i+1L]])# + random_dist_mat[current_order[i+1L], vk]) / 2
    d_ij <- (random_dist_mat[current_order[i], current_order[i+1L]])# + random_dist_mat[current_order[i], current_order[i+1L]]) / 2
    
    #print(d_ik)
    #print(d_kj)
    #print(d_ij)
    costs[i] <- d_ik + d_kj - d_ij
  }
  #print(current_order)
  #print(vk)
  costs <- c(random_dist_mat[vk, current_order[1]], costs, random_dist_mat[current_order[length], vk])
  #print(costs)
  #print(Sys.time() - t1)
  return(costs)
  
}
#############################################################################

path_length <-function(matrix, order){
  l1 <- sum(matrix[cbind(order[-nrow(matrix)], order[-1])])
  
  return(l1)
}

#dyn.load("path_insertion_cost.dll")

#' Hamiltonian Path - Modified Arbitrary Insertion
#' 
#' This function calculates the Hamiltonian path solution with a modified arbitrary insertion algorithm - designed for asymmetric cases (eg a KL distance matrix)
#'@param dist_mat A N by N distance matrix (KL distance in this case). Works on any other asymmetrical / symmetrical distance matrices
#'@keywords hamiltonian_path_as
#'@export
#'@section redPATH pseudotime:
#'@examples
#'hamiltonian_path_as(kl_cpp(data))
#'
hamiltonian_path_as <- function(dist_mat, start = NULL)
{
  dist_mat <- as.matrix(dist_mat)
  n_pts <- nrow(dist_mat)
  #result_list <- matrix(0, nrow = n_pts, ncol = n_pts)
  result_list <- rep(0, n_pts)
  
  # sets the ordering variable that sorts the random_order to give final path
  ordering <- integer(n_pts+1)
  
    
  random_order <- sample(n_pts)
  ordering[1:2] <- 1:2

  
  # proposed randomly ordered dist mat:
  random_dist <- dist_mat[random_order, random_order]

  #print("begin")
  # arbitrary insertion
  for(i in c(3:n_pts)){
    #insert_pos <- which.min(calc_insertion_cost(random_dist, ordering[1:(i-1L)], i)) + 1
    forward_cost <- .Call("path_insertion_cost", random_dist, ordering[1:(i-1L)], as.integer(i))

    insert_pos <- which.min(forward_cost)
    if(insert_pos == 1){
      insert_pos = 1
      ordering[((insert_pos):i)+1L] <- ordering[(insert_pos):i]
      ordering[insert_pos] <- i
      #print(ordering)
    }else if(insert_pos == (i)){
      insert_pos = i
      ordering[i] <- insert_pos
      #print(ordering)
    }
    else{
      insert_pos <- insert_pos
      ordering[((insert_pos):i)+1L] <- ordering[(insert_pos):i]
      ordering[insert_pos] <- i
    }

  }
  return(random_order[ordering])
}

