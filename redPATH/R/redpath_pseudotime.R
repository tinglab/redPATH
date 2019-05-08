#' redPATH Pseudotime
#' 
#' This function performs the redPATH algorithm for single cell gene expression data, obtaining the pseudo development time of each single cell.
#'@param cell_gene_exprs The gene expression input obtained from get_GO_exp
#'@param threadnum The number of threads to use
#'@param step_size The interval of paths to merge in the final solution. Defaults to 1.
#'@param base_path_range Sets the range of any consequtive 5 cluster solutions to calculate the base path. The maximum number can be 10 (ie c(6:10)). Defaults to c(3:7).
#'@param clustMethod = "Corr" This sets the clustering method used. Defaults to "Corr", which is hierarchical clustering on doubly centered matrix. "GMM" is also available for Gaussian Mixture Model.
#'@export
#'@keywords redpath
#'@section redPATH pseudotime:
#'@examples
#'redpath_pseudotime <- redpath(data, threadnum = 4) 

redpath <- function(cell_gene_exprs, threadnum, step_size = 1, base_path_range = c(3:7), dist_type="KL2", clustMethod = "Corr")
{
  max_loop = NULL
  
  result <- bestEnsembleHMT(cell_gene_exprs = cell_gene_exprs, HMT_iterations = 2, beginNum = max(base_path_range), endNum = dim(cell_gene_exprs)[1], threads = threadnum, clustMethod = "Corr", 
                            dist_type = dist_type, base_path_range = base_path_range, step_size = step_size, max_loop = max_loop, max_num = 500)
  pseudotime <- result$ensemble_HMT_list[dim(result$ensemble_HMT_list)[1], ]
  #print(pseudotime)
  
  return(pseudotime)
}