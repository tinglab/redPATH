# POS score
# Adapted from TSCAN R Package 
# Zhicheng Ji and Hongkai Ji. (2016)
# TSCAN: Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis. (2016) Nucleic Acids Research, 44(13):e117.

pos_score <- function(cell_labels, pseudotime, order = F) {
  if(order == F){
    pseudotime_ordering <- order(pseudotime)
  }else{
    pseudotime_ordering <- pseudotime
  }
  cell_labels <- as.numeric(cell_labels)
    
  score_order <- cell_labels[pseudotime_ordering]
  opt_score_order <- sort(score_order)
  opt_score <- sum(sapply(1:(length(opt_score_order)-1),function(x) {
    sum(opt_score_order[(x+1):length(opt_score_order)] - opt_score_order[x])
  })) 
  final <- sum(sapply(1:(length(score_order)-1),function(x) {
    sum(score_order[(x+1):length(score_order)] - score_order[x])
  })) / opt_score
  return(abs(final))
}      