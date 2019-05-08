# Evaluation Metrics

#' Change Index
#' 
#' It calculates the change index of the ordered cells. Evaluation in the range of [0, 1] with 1 being perfectly aligned.
#'@param cell_labels The cell type labels
#'@param pseudotime The pseudo time / order index calculated from various softwares.
#'@param order Defaults to FALSE for input of pseudo-time, TRUE for input of an order index of cells
#'@export
#'@keywords change_index
#'@section Evaluation:
#'@examples
#'ci <- change_index(cell_labels, redpath_pseudotime, order = F) 
change_index <- function(cell_labels, pseudotime, order = F){
  if(order == F){
    timeSeries <- cell_labels[order(pseudotime)]
  }else if(order == T){
    timeSeries <- cell_labels[pseudotime]
  }

  sumChange <- 0
  for(i in 1:(length(timeSeries)-1))
  {
    if(timeSeries[i]!=timeSeries[i+1])
    {
      sumChange <- sumChange+1
    }
    if(i==length(timeSeries) & timeSeries[i]!=timeSeries[1])
    {
      sumChange <- sumChange+1
    }
  }
  N <- length(timeSeries)
  
  #change_index <- 1-(sumChange - 1)/(N-2)
  unique_states <- length(unique(timeSeries))
  change_index <- 1-((sumChange - unique_states + 1)/(N-unique_states))
  print(paste("Total changes: ", sumChange, sep = ""))
  return(change_index)
}

bubblesort <- function(lst){
  length_lst <- length(lst)
  num_swaps <- 0
  for(i in 1:(length_lst)){
    for(j in 1:(length_lst-1)){
      
      if(lst[j+1] < lst[j]){
        num_swaps <- num_swaps + 1
        temp <- lst[j]
        lst[j] <- lst[j+1]
        lst[j+1] <- temp
      }
    }
    #print(lst)
  }
  #print(num_swaps)
  return(lst)
}
bubblesort_function <- function(lst){
  length_lst <- length(lst)
  num_swaps <- 0
  
  for(i in 1:(length_lst)){
    for(j in 1:(length_lst-1)){
      
      if(lst[j+1] < lst[j]){
        num_swaps <- num_swaps + 1
        temp <- lst[j]
        lst[j] <- lst[j+1]
        lst[j+1] <- temp
      }
    }
    #print(lst)
  }
  #print(num_swaps)
  return(num_swaps)
}
generate_worst_case <- function(lst){
  length_lst <- length(lst)
  unique_type <- unique(lst)
  sorted_unique_type <- sort(unique_type)
  n_type <- length(unique_type)
  
  output_lst <- c()
  for(i in c(1:n_type)){
    l_temp <- length(which(lst == sorted_unique_type[i]))
    temp <- rep(sorted_unique_type[i], l_temp)
    output_lst <- c(output_lst, temp)
  }
  return(rev(output_lst))
}

#' Bubblesort Index
#' 
#' It calculates the bubblesort index of the ordered cells. Evaluation in the range of [0, 1] with 1 being perfectly aligned.
#'@param cell_labels_int The cell type labels written in terms of int, with ordering of cell types. For example, qNSC -> aNSC -> NPC becomese 1 -> 2 -> 3.
#'@param pseudotime The pseudo time / order index calculated from various softwares.
#'@param order Defaults to FALSE for input of pseudo-time, TRUE for input of an order index of cells
#'@export
#'@keywords bubblesort_index
#'@section Evaluation:
#'@examples
#'bs <- bubblesort_index(cell_labels_int, redpath_pseudotime, order = F) 
#'
bubblesort_index <- function(cell_labels_int, pseudotime, order = F){
  if(order == F){
    lst <- cell_labels_int[order(pseudotime)]
  }else if(order == T){
    lst <- cell_labels_int[pseudotime]
  } 
  steps_used <- bubblesort_function(lst)
  worst_case <- bubblesort_function(generate_worst_case(lst))
  ratio <- steps_used / worst_case
  ratio2 = 1 - ratio
  if(ratio >= ratio2){
    return(ratio)
  }else{
    return(ratio2)
  }

}
#' Kendall Correlation
#' 
#' It calculates the highest kendall correlation of the ordered cells. Evaluation in the range of [-1, 1] with 1 being perfectly aligned.
#'@param cell_labels_int The cell type labels written in terms of int, with ordering of cell types. For example, qNSC -> aNSC -> NPC becomese 1 -> 2 -> 3.
#'@param pseudotime The pseudo time / order index calculated from various softwares.
#'@param order Defaults to FALSE for input of pseudo-time, TRUE for input of an order index of cells
#'@export
#'@keywords kendall_correlation
#'@section Evaluation:
#'@examples
#'kc <- kendall_correlation(cell_labels_int, redpath_pseudotime, order = F) 
#'
kendall_correlation <- function(cell_labels_int, pseudotime, order = F){
  if(order == F){
    timeSeries <- cell_labels_int[order(pseudotime)]
  }else if(order == T){
    timeSeries <- cell_labels_int[pseudotime]
  } 
  kc1 <- cor(as.numeric(timeSeries), rev(as.numeric(bubblesort(timeSeries))), method = "kendall")
  if(kc1 < 0){
    kc1 <- cor(as.numeric(timeSeries), (as.numeric(bubblesort(timeSeries))), method = "kendall")
  }
  print(kc1)
  return(kc1)
}

# POS score
# Adapted from TSCAN R Package 
# Zhicheng Ji and Hongkai Ji. (2016)
# TSCAN: Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis. (2016) Nucleic Acids Research, 44(13):e117.
#' POS Score
#' 
#' Adapted from Zhicheng Ji and Hongkai Ji. (2016)
#' TSCAN: Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis. (2016) Nucleic Acids Research, 44(13):e117.
#' Calculates the POS score for evaluating the pseudotime alignment. This is modified to give a value between [0, 1]. 
#'@param cell_labels_int The cell type labels written in terms of int, with ordering of cell types. For example, qNSC -> aNSC -> NPC becomese 1 -> 2 -> 3.
#'@param pseudotime The pseudo time / order index calculated from various softwares.
#'@param order Defaults to FALSE for input of pseudo-time, TRUE for input of an order index of cells
#'@export
#'@keywords pos_score
#'@section Evaluation:
#'@examples
#'pos <- pos_score(cell_labels_int, redpath_pseudotime, order = F) 
#'
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