# Post-Biological Analysis: Utility Functions


#' Calculate correlation of gene expression with pseudotime 
#' 
#' This function calculates the correlation of gene expression with the inferred pseudo time.
#'@param full_gene_exprs A gene (row) by cell (col) gene expression matrix with all genes. 
#'@param pseudo_time The pseudo development time / ordering of cells can be provided as input. Please leave order = F if pseudotime is provided
#'@param type The type of correlation to calculate. "dcor" and "mic" are available (Distance correlation / Maximal Information Criterion).
#'@param order Defaults to F. If pseudotime is provided, set to F. If an ordering of cells is provided, set to T.
#'@export
#'@keywords calc_correlation
#'@section Biological Analysis:
#'@examples
#'dcor_result <- calc_correlation(full_exprs_matrix, redpath_pseudotime, type = "dcor") 
#'mic_result <- calc_correlation(full_exprs_matrix, redpath_pseudotime, type = "mic")
#'selected_result <- gene_selection(full_exprs_matrix, redpath_pseudotime, dcor_result, mic_result, threshold = 0.5)

calc_correlation <- function(full_gene_exprs, pseudo_time, type = "dcor", order = F){ # dcor / mic
  require(minerva)
  require(energy)
  if(order == T){
    sorted_gene_exprs <- full_gene_exprs[,pseudo_time]
    #pt <- pseudo_time
  }else{
    sorted_gene_exprs <- full_gene_exprs[ , order(pseudo_time)]
  }
  pt <- pseudo_time[order(pseudo_time)]

  if(type == "dcor"){
    #dcor_gene_vector <- apply(sorted_gene_exprs, 1, function(x) dcor(x = x, y = c(1:dim(sorted_gene_exprs)[2])))
    dcor_gene_vector <- apply(sorted_gene_exprs, 1, function(x) dcor(x = x, y = pt))
    
    print("Summary Statistics: ")
    print(summary(dcor_gene_vector))
    return(dcor_gene_vector)
  }else if (type == "mic"){
    mic_gene_vector = apply(sorted_gene_exprs, 1, function(x) mine(x = x, y = c(1:dim(sorted_gene_exprs)[2]))$MIC)
    print("Summary Statistics: ")
    print(summary(mic_gene_vector))
    return(mic_gene_vector)
  }
}

#' Selects significant genes
#' 
#' This function selects the significant genes from either 1 or 2 correlation vectors
#'@param full_gene_exprs A gene (row) by cell (col) gene expression matrix with all genes. 
#'@param pseudo_time The pseudo development time 
#'@param correlation_vector The calculated correlation vector of each gene (from calc_correlation)
#'@param correlation_vector2 Optional vector to select the significant genes. Defaults to NULL. If provided, genes will be selected if the threshold satisfies either correlation vectors. 
#'@param threshold The threshold cut-off to retain significant genes. Defaults to 0.5
#'@export
#'@keywords gene_selection
#'@section Biological Analysis:
#'@examples
#'dcor_result <- calc_correlation(full_exprs_matrix, redpath_pseudotime, type = "dcor") 
#'mic_result <- calc_correlation(full_exprs_matrix, redpath_pseudotime, type = "mic")
#'selected_result <- gene_selection(full_exprs_matrix, redpath_pseudotime, dcor_result, mic_result, threshold = 0.5)

gene_selection <- function(full_gene_exprs, pseudo_time, correlation_vector, correlation_vector2 = NULL, threshold = 0.5){
  sorted_gene_exprs <- full_gene_exprs[ , order(pseudo_time)]
  
  if(is.null(correlation_vector2)){
    correlation_vector <- signif(correlation_vector, 2)
    selection <- which(correlation_vector >= threshold)
    selection <- selection[order(schlitz_mic[selection])]
    selected_gene_exprs <- sorted_gene_exprs[selection, ]
    return(selected_gene_exprs)
  }else{
    correlation_vector <- signif(correlation_vector, 2)
    correlation_vector2 <- signif(correlation_vector2, 2)
    selection_1 <- which(correlation_vector >= threshold)
    selection_2 <- which(correlation_vector2 >= threshold)
    union_selection <- c(selection_1, selection_2)
    non_duplicated_genes <- which(!duplicated(union_selection))
    overall_selection <- union_selection[non_duplicated_genes]
    
    combined_val <- (correlation_vector[overall_selection] + correlation_vector2[overall_selection] ) / 2
    overall_selection <- overall_selection[order(combined_val)]
    #selected_gene_exprs <- sorted_gene_exprs[union_selection[non_duplicated_genes], ]
    selected_gene_exprs <- sorted_gene_exprs[overall_selection, ]
    return(selected_gene_exprs)
  }
}


#' Selects significant genes
#' 
#' This function selects the significant genes from 2 correlation vectors (intersecting genes)
#'@param full_gene_exprs A gene (row) by cell (col) gene expression matrix with all genes. 
#'@param pseudo_time The pseudo development time 
#'@param correlation_vector The calculated correlation vector of each gene (from calc_correlation)
#'@param correlation_vector2 Genes will be selected if the threshold satisfies BOTH correlation vectors. 
#'@param threshold The threshold cut-off to retain significant genes. Defaults to 0.5
#'@export
#'@keywords gene_selection_intersection
#'@section Biological Analysis:
#'@examples
#'dcor_result <- calc_correlation(full_exprs_matrix, redpath_pseudotime, type = "dcor") 
#'mic_result <- calc_correlation(full_exprs_matrix, redpath_pseudotime, type = "mic")
#'selected_result <- gene_selection(full_exprs_matrix, redpath_pseudotime, dcor_result, mic_result, threshold = 0.5)

gene_selection_intersection <- function(full_gene_exprs, pseudo_time, correlation_vector, correlation_vector2 = NULL, threshold = 0.5){
  sorted_gene_exprs <- full_gene_exprs[ , order(pseudo_time)]
  
  if(is.null(correlation_vector2)){
    correlation_vector <- signif(correlation_vector, 2)
    selection <- which(correlation_vector >= threshold)
    selection <- selection[order(schlitz_mic[selection])]
    selected_gene_exprs <- sorted_gene_exprs[selection, ]
    return(selected_gene_exprs)
  }else{
    correlation_vector2 <- signif(correlation_vector2, 2)
    correlation_vector <- signif(correlation_vector, 2)
    selection_1 <- which(correlation_vector >= threshold)
    selection_2 <- which(correlation_vector2 >= threshold)
    union_selection <- c(selection_1, selection_2)
    duplicated_genes <- which(duplicated(union_selection))
    overall_selection <- union_selection[duplicated_genes]
    
    combined_val <- (correlation_vector[overall_selection] + correlation_vector2[overall_selection] ) / 2
    overall_selection <- overall_selection[order(combined_val)]
    #selected_gene_exprs <- sorted_gene_exprs[union_selection[non_duplicated_genes], ]
    selected_gene_exprs <- sorted_gene_exprs[overall_selection, ]
    return(selected_gene_exprs)
  }
}

calc_all_hmm_diff <- function(sorted_gene_exprs_matrix, cls_num){
  total_cells <- ncol(sorted_gene_exprs_matrix)
  total_genes <- nrow(sorted_gene_exprs_matrix)
  output <- matrix(0, nrow = total_genes, ncol = total_cells)
  for(i in c(1:total_genes)){
    input_exprs <- as.numeric(sorted_gene_exprs_matrix[i, ])
    tmp_hmm <- get_HMM_order_diff(input_exprs, c(1:total_cells), cls_num)
    output[i, ] <- tmp_hmm
    print(i)
  }
  rownames(output) <- rownames(sorted_gene_exprs_matrix)
  return(output)
  
}

#source("D:/Single Cells/Single Cell Analysis/HMM code/baum_welch_code.R")
#source("D:/Single Cells/Single Cell Analysis/HMM code/HMM.R")
#source("D:/Single Cells/Single Cell Analysis/HMM code/get_hmm_code.R")


#' HMM segmentation for each gene
#' 
#' HMM is inferred for a ON/HIGH and OFF/LOW state for each significant gene along a developmental pseudo time. Extension to a 3-state model is also provided.
#'@param sorted_gene_exprs_matrix Result returned from gene_selection function. A re-ordered gene by cell matrix
#'@param cls_num The number of states for HMM. Defaults to 2. 
#'@param threadnum The number of threads to use.
#'@export
#'@keywords gene_selection_intersection
#'@section Biological Analysis:
#'@examples
#'dcor_result <- calc_correlation(full_exprs_matrix, redpath_pseudotime, type = "dcor") 
#'mic_result <- calc_correlation(full_exprs_matrix, redpath_pseudotime, type = "mic")
#'selected_result <- gene_selection(full_exprs_matrix, redpath_pseudotime, dcor_result, mic_result, threshold = 0.5)
#'gene_state_result <- calc_all_hmm_diff_parallel(selected_result, cls_num = 2, threadnum = 4)

calc_all_hmm_diff_parallel <- function(sorted_gene_exprs_matrix, cls_num = 2, threadnum = 1){
  require(doParallel)
  total_cells <- ncol(sorted_gene_exprs_matrix)
  total_genes <- nrow(sorted_gene_exprs_matrix)
  #output <- matrix(0, nrow = total_genes, ncol = total_cells)
  seq_list <- c(1:total_genes)
  cls_num <- cls_num
  
  cl <- makeCluster(threadnum)
  registerDoParallel(cl)
  
  output <- foreach(i = seq_list, .combine = 'rbind', .export = c("get_HMM_order_diff", "normalDis", "myBW", "eln", 
                                                                  "myviterbi", "calc_lnb", "get_Ndpara", "get_A",
                                                                  "get_Pi", "eln_Xi", "eln_Gamma", "elnproduct",
                                                                  "eln_Alpha", "elnsum", "eln_Beta", "eexp",
                                                                  "calcEmissionProbs", "vmulP", "transProb")) %dopar%{
    input_exprs <- as.numeric(sorted_gene_exprs_matrix[i, ])
    #print("yes")
    tmp_hmm <- get_HMM_order_diff(exprs_val = input_exprs, ordIndex = c(1:total_cells), cls_num = cls_num)
    #output[i, ] <- tmp_hmm
    return(tmp_hmm)
  }
  stopCluster(cl)
  closeAllConnections()
  
  output <- as.matrix(output) 
  #for(i in c(1:total_genes)){
  #  input_exprs <- as.numeric(sorted_gene_exprs_matrix[i, ])
  #  tmp_hmm <- get_HMM_order_diff(input_exprs, c(1:total_cells), cls_num)
  #  output[i, ] <- tmp_hmm
  #}
  rownames(output) <- rownames(sorted_gene_exprs_matrix)
  return(output)
  
}

confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

