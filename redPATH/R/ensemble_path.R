

#source("high_cor_path.R")



round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}


bestEnsembleHMT <- function(cell_gene_exprs, stageIdxSmp = NULL, HMT_iterations, beginNum = 10, 
                            endNum = 198, threads = 1, clustMethod = "Corr", dist_type = "KL2", 
                            debug = FALSE, base_path_range = c(7:10), step_size = 2, max_loop = NULL, 
                            max_num = 300)
{
  #require(RcppAlgos)
  #require(mclust)
  require(doParallel)
  message("Calculating redPATH...")
  c1 <- makeCluster(threads)
  
  #HMT_list <- matrix(0, nrow = (endNum - beginNum +1), ncol = nrow(cell_gene_exprs))
  tmp_HMT_list <- matrix(0, length(base_path_range), ncol = nrow(cell_gene_exprs))
  tmp_norm_HMT_list <- matrix(0, length(base_path_range), ncol = nrow(cell_gene_exprs))
  
  #normalized_HMT_list <- matrix(0, nrow = (endNum - beginNum + 1) , ncol = nrow(cell_gene_exprs))
  #ensemble_HMT_list <- matrix(0, nrow = (endNum - beginNum + 1), ncol = nrow(cell_gene_exprs))
  beginNum <- max(base_path_range)
  if(endNum > 500){
    endNum <- max_num
  }
  if(is.null(max_loop)){
    s1 <- (endNum - beginNum) %% step_size
    total_rows <- round2(((endNum - beginNum -s1) / step_size), 0) + 1
    HMT_list <- matrix(0, nrow = total_rows, ncol = nrow(cell_gene_exprs))
    normalized_HMT_list <- matrix(0, nrow = total_rows, ncol = nrow(cell_gene_exprs))
    ensemble_HMT_list <- matrix(0, nrow = total_rows, ncol = nrow(cell_gene_exprs))
    seq_list <- seq(from = beginNum + step_size, to = endNum, by = step_size)
  }
  else if (!is.null(max_loop))
  {
    s1 <- (endNum - max_loop) %% step_size
    max_single_step <- max_loop - beginNum
    total_rows <- (max_loop - beginNum) + round2(((endNum - max_loop - s1) / step_size), 0) + 1
    HMT_list <- matrix(0, nrow = total_rows, ncol = nrow(cell_gene_exprs))
    normalized_HMT_list <- matrix(0, nrow = total_rows, ncol = nrow(cell_gene_exprs))
    ensemble_HMT_list <- matrix(0, nrow = total_rows, ncol = nrow(cell_gene_exprs))
    seq_list <- c((beginNum+1):max_loop, seq(from = max_loop+step_size, to = endNum, by = step_size))
  }
  
  
  highCorLst <- c()
  finalHighCorLst <- c()
  
  store_distance_list <- list()
  
  # Calculate Base Score
  base_path_range <- base_path_range
  

  for (clust_num_EM in base_path_range) 
  {
    i <- clust_num_EM - (min(base_path_range)-1) 
    best_cor_path <- high_cor_path(clust_num_EM, NULL, cell_gene_exprs, stageIdxSmp, HMT_iterations = HMT_iterations, TSPMethod = "force", clustMethod = clustMethod, dist_type = dist_type, debug = debug)

    tmp_HMT_list[i, ] <- best_cor_path$result
    tmp_norm_HMT_list[i, ] <- best_cor_path$norm_result

  }
  message("Finished calculating base path...")



#############################################################################################
  tmp_HMT_list <- calc_base_path(tmp_norm_HMT_list, tmp_HMT_list, base_path_range)
###########################################################################################
  
  baseResult <- colMeans(tmp_HMT_list)
  baseResult <- normalizeResult(baseResult)

  #print(baseResult)
  #print(tmp_HMT_list[1:4,])
  
  normalized_HMT_list[1, ] <- baseResult
  ensemble_HMT_list[1, ] <- baseResult
  
  
  if (threads <= 1)
  {
    # sequential processing
    #HMT_list[1, ] <- colMeans(tmp_HMT_list)
    HMT_list <- baseResult
    i <- 2
    for (clust_num_EM in seq_list)
    {
      #i <- clust_num_EM - beginNum + 1
      #dyn.load("path_insertion_cost.dll")
      #Rcpp::sourceCpp("kl_opt.cpp")
      #source("kl_dist_optimized.R")
      require(mclust)
      require(RcppAlgos)
      best_cor_path <- high_cor_path(clust_num_EM, baseResult, cell_gene_exprs, stageIdxSmp, HMT_iterations = HMT_iterations, clustMethod = clustMethod, dist_type = dist_type, debug = debug)
      #HMT_list[i, ] <- best_cor_path$norm_result
      HMT_list <- rbind(HMT_list, best_cor_path$norm_result)
      i <- i + 1
    }
  }
  else
  {
    # parallel processing
    cl <- makeCluster(threads)
    registerDoParallel(cl)
    require(mclust)
    require(RcppAlgos)
    HMT_list <- foreach(clust_num_EM = seq_list, .combine = 'rbind', .export = c("high_cor_path","permuteGeneral", "hamiltonian_path_force", "corr_dist", "normalizeResult2", 
                                                                                 "Mclust", "mclustBIC","hamiltonian_path_as", "path_length", "ED2")) %dopar%
    {
      #dyn.load("path_insertion_cost.dll")
      #Rcpp::sourceCpp("kl_opt.cpp")
      #source("kl_dist_optimized.R")
      
      #dyn.load("../src/redPATH.dll")
      #dyn.load("../../force_trial.dll")
      #Rcpp::sourceCpp("../src/kl_opt.cpp")
      #source("kl_dist_optimized.R")
      #source("../../ai_trial.R")
      #source("high_cor_path.R")
      
      best_cor_path <- high_cor_path(clust_num_EM, baseResult, cell_gene_exprs, stageIdxSmp, HMT_iterations = HMT_iterations, clustMethod = clustMethod, dist_type = dist_type, debug = debug)
      #return(best_cor_path$result)
      return(best_cor_path$norm_result)
    }
    stopCluster(cl)
    closeAllConnections()
    HMT_list <- rbind(baseResult, HMT_list)

  }

  i <- 2
  

  message("Finished calculating all paths...")

  
  for (clust_num_EM in seq_list)
  {

    
    resulti <- HMT_list[i, ]
    moveCorLst <- 0
    revCorLst <- 0
    
    tmp <- resulti
    moveCorLst <- cor(ensemble_HMT_list[1, ], tmp, method = "spearman")
    #moveCorLst <- sum(.Call("path_length_c", kl_dist, order(tmp)))
    
    tmp <- 1 - resulti #rev(resulti)
    revCorLst <- cor(ensemble_HMT_list[1, ], tmp, method = "spearman")
    #revCorLst <- sum(.Call("path_length_c", kl_dist, order(tmp)))
    
    if(moveCorLst >= revCorLst){
    #if(moveCorLst <= revCorLst){
      normalized_HMT_list[i, ] <- resulti
    }else{
      normalized_HMT_list[i, ] <- (1)-resulti
    }
    
    i <- i + 1
  }
  
  assign("norm_res_new", normalized_HMT_list, .GlobalEnv)
  
  assign("norm_new", normalized_HMT_list, .GlobalEnv)
  i <- 2

  for (i in 2:nrow(normalized_HMT_list))
  {
    meanResults <- colMeans(matrix(normalized_HMT_list[1:i, ], nrow = i))
    ensemble_HMT_list[i, ] <- meanResults
  }


  #save(finalHighCorLst, highCorLst, HMT_list, normalized_HMT_list, ensemble_HMT_list, file = paste(getwd(), "/","best_KL_", beginNum, "-", endNum, "_", "paths.RData", sep = ""))
  #print(range(ensemble_HMT_list))
  return(list(finalHighCorLst = finalHighCorLst, ensemble_HMT_list = ensemble_HMT_list))
  
}
