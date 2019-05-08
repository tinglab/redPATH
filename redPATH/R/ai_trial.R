
#' Hamiltonian Path - Modified Arbitrary Insertion
#' 
#' This function calculates the Hamiltonian path solution with a modified arbitrary insertion algorithm - designed for asymmetric cases (eg a KL distance matrix)
#'@param dist_mat A N by N distance matrix (KL distance in this case). Works on any other asymmetrical / symmetrical distance matrices
#'@keywords hamiltonian_path_as
#'@export
#'@section redPATH pseudotime:
#'@examples
#'se <- get_start_end(kl_cpp(data))
#'hamiltonian_path_ai(kl_cpp(data), se)
#'

hamiltonian_path_ai <- function(dist_mat, se)
{
  dist_mat <- as.matrix(dist_mat)
  n_pts <- nrow(dist_mat)
  #result_list <- matrix(0, nrow = n_pts, ncol = n_pts)
  result_list <- rep(0, n_pts)
  
  # sets the ordering variable that sorts the random_order to give final path
  ordering <- integer(n_pts+1)
  
  
  random_order <- c(1:n_pts)#sample(n_pts)
  #ordering[1:2] <- 1:2
  
  
  # proposed randomly ordered dist mat:
  random_dist <- dist_mat[random_order, random_order]
  #se <- get_start_end(random_dist, k = 7)
  ordering[1:2] <- 1:2#as.integer(se)
  
  total_space <- c(1:n_pts)
  search_space <- total_space[-se]
  
  search_space <- sample(search_space)
  #print(str(search_space))
  random_order <- c(se, search_space)
  random_dist <- dist_mat[random_order, random_order]
  #print("begin")
  # arbitrary insertion
  for(i in c(3:n_pts)){
    #j <- 3
    #for(i in search_space){
    #insert_pos <- which.min(calc_insertion_cost(random_dist, ordering[1:(i-1L)], i)) + 1
    forward_cost <- .Call("path_insertion_cost", random_dist, ordering[1:(i-1L)], as.integer(i))
    
    insert_pos <- which.min(forward_cost)
    #print(insert_pos)
    if(insert_pos == 1){
      insert_pos = 1
      ordering[((insert_pos):i)+1L] <- ordering[(insert_pos):i]
      ordering[insert_pos] <- i
      #print(ordering)
    }else if(insert_pos == (i)){
      insert_pos = i
      ordering[i] <- as.integer(insert_pos)
      #print(ordering)
    }
    else{
      insert_pos <- insert_pos
      ordering[((insert_pos):i)+1L] <- ordering[(insert_pos):i]
      ordering[insert_pos] <- i
    }
    # <- j + 1
  }
  #print(random_order[ordering])
  return(random_order[ordering])
}


#' Get start & end of path
#' 
#' This function estimates the initial starting and end node as input to the arbitrary insertion algorithm
#'@param dist_mat A N by N distance matrix (KL distance in this case). Works on any other asymmetrical / symmetrical distance matrices
#'@keywords get_start_end
#'@export
#'@section redPATH pseudotime:
#'@examples
#'se <- get_start_end(kl_cpp(data))
#'hamiltonian_path_ai(kl_cpp(data), se)
#'
get_start_end <- function(dist_mat, k = 6){
  
  dist_mat <- as.matrix(dist_mat)
  n_pts <- nrow(dist_mat)
  
  k_solns <- n_pts - k +1
  tmp_soln_matrix <- matrix(0, nrow = k_solns, ncol = k)
  # proposed randomly ordered dist mat:
  # random_dist <- dist_mat[random_order, random_order]
  #print(tmp_soln_matrix)
  
  for(i in c(1:k_solns)){
    selection <- c(i:(i+k-1))
    tmp_soln <- hamiltonian_path_force(dist_mat[selection, selection])[-1]
    #print(tmp_soln)
    #tmp_soln_matrix[i, ] <- c(1:n_pts)[selection][tmp_soln]#tmp_soln
    tmp_soln_matrix[i, ] <- selection[tmp_soln]
    #print(c(1:n_pts)[selection][tmp_soln])
    #print(path_length(dist_mat[selection, selection], tmp_soln))
  }
  candidate_start <- unique(tmp_soln_matrix[,1])
  
  #sequential_rank <- c(1:k_solns)
  store_tmp_start <- c()
  for(j in c(1:length(candidate_start))){
    rank_j <- apply(data.frame(tmp_soln_matrix), 1, function(x) sum(which(x == j)))
    #print(rank_j)
    seq_rank_j <- which(rank_j != 0)
    #print(seq_rank_j)
    rank_j <- rank_j[which(rank_j !=0)]
    pos_j <- sum(rank_j * seq_rank_j)
    
    tmp_value <- pos_j / (length(which(rank_j !=0)) * j)
    store_tmp_start <- c(store_tmp_start, tmp_value)
  }
  #print(store_tmp_start)
  #print(candidate_start)
  start_value <- candidate_start[which.min(store_tmp_start)]
  #print(paste0("start value ", start_value))
  
  candidate_end <- unique(tmp_soln_matrix[,k])
  store_tmp_end <- c()
  for(k in c(1:length(candidate_end))){
    rank_k <- apply(data.frame(tmp_soln_matrix), 1, function(x) sum(which(x == k)))
    
    #print(rank_k)
    seq_rank_k <- which(rank_k != 0)
    #print(seq_rank_k)
    rank_k <- rank_k[which(rank_k !=0)]
    
    pos_k <- sum(rank_k * seq_rank_k)
    
    tmp_value <- pos_k / (length(which(rank_k !=0)) * k)
    store_tmp_end <- c(store_tmp_end, tmp_value)
  }
  #print(candidate_end)
  #print(store_tmp_end)
  end_value <- candidate_end[which.max(store_tmp_end)]
  if(end_value == start_value){
    end_value <- candidate_end[-which.max(store_tmp_end)][which.min(store_tmp_end[-which.max(store_tmp_end)])]
  }
  #print(paste0("start value ", end_value))
  return(c(start_value, end_value))
}



novel_hmt <- function(dist_mat, k = 5){
  dist_mat <- as.matrix(dist_mat)
  n_pts <- nrow(dist_mat)
  #result_list <- matrix(0, nrow = n_pts, ncol = n_pts)
  result_list <- rep(0, n_pts)
  
  # sets the ordering variable that sorts the random_order to give final path
  ordering <- integer(n_pts+1)
  
  
  random_order <- sample(n_pts)
  ordering[1:2] <- 1:2
  
  print(random_order)
  k_solns <- n_pts - k +1
  tmp_soln_matrix <- matrix(0, nrow = k_solns, ncol = k)
  # proposed randomly ordered dist mat:
  random_dist <- dist_mat[random_order, random_order]
  print(tmp_soln_matrix)
  
  for(i in c(1:k_solns)){
    selection <- c(i:(i+k-1))
    tmp_soln <- hamiltonian_path_force(dist_mat[selection, selection])[-1]
    #print(tmp_soln)
    tmp_soln_matrix[i, ] <- c(1:n_pts)[selection][tmp_soln]#tmp_soln
    print(c(1:n_pts)[selection][tmp_soln])
    print(path_length(dist_mat[selection, selection], tmp_soln))
  }
  
  print(tmp_soln_matrix)
  
  # merge solutions
  
}


get_estimation <- function(dist_mat, k = 6){
  
  dist_mat <- as.matrix(dist_mat)
  n_pts <- nrow(dist_mat)
  
  k_solns <- n_pts - k +1
  tmp_soln_matrix <- matrix(0, nrow = k_solns, ncol = k)
  # proposed randomly ordered dist mat:
  # random_dist <- dist_mat[random_order, random_order]
  #print(tmp_soln_matrix)
  
  for(i in c(1:k_solns)){
    selection <- c(i:(i+k-1))
    tmp_soln <- hamiltonian_path_force(dist_mat[selection, selection])[-1]
    #print(tmp_soln)
    tmp_soln_matrix[i, ] <- selection[tmp_soln]#tmp_soln
    #print(c(1:n_pts)[selection][tmp_soln])
    #print(path_length(dist_mat[selection, selection], tmp_soln))
  }
  
  candidates <- c(1:nrow(dist_mat))
  store_soln <- c()
  for(j in c(1:length(candidates))){
    rank_j <- apply(data.frame(tmp_soln_matrix), 1, function(x) sum(which(x ==j)))
    print(rank_j)
    seq_rank_j <- which(rank_j != 0)
    
    rank_j <- rank_j[which(rank_j !=0)]
    pos_j <- sum(rank_j * seq_rank_j)
    
    tmp_value <- pos_j / (length(which(rank_j !=0 )) * j)
    store_soln <- c(store_soln, tmp_value)
    
  }

  #print(paste0("start value ", end_value))
  return(store_soln)
}
