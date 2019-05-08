# Identification of G0-like cells
######################################################################################
require(car)
require(dplyr)

#' Kmeans Clustering - returns classification result
#' 
#' This function performs imple kmeans clustering and prepares as input for identification of G0 cells
#'@param mean_scores The mean score results returned from the reCAT software, "get_score"
#'@param k The number of clusters, defaults to 5. ie G0, G1, S, G2, M cycle stages
#'@param nstart Defaults to 100. Number of times to start the kmeans clustering, ensure stability of clustering results.
#'@export
#'@keywords kmeans_clustering
#'@section G0 Identification:
#'@examples
#'cycle_mean_scores <- get_score(t(test_exp))$mean_score
#'kmean_classification <- kmeans_clustering(cycle_mean_scores)
#'statistical_test(cycle_mean_scores, kmean_classification, threshold = 0.001)
#'# Sample output
#'#"Possible group of G0-like cells is: 2"
kmeans_clustering <- function(mean_scores, k = 5, nstart = 100){
  km <- kmeans(mean_scores, centers = k, nstart = nstart)
  km_result <- km$cluster
  return(km_result)
}
#' Identification of G0-like cells
#' 
#' This function performs the simple statistical tests to determine the group of G0 like cells within single cell datasets.
#'@param mean_scores The mean score results returned from the reCAT software, "get_score"
#'@param kmean_result The classification of kmeans clustering result.
#'@param threshold Defaults to 0.001. The criterion that all p-values must be less than the threshold in order to determine G0-like cells
#'@export
#'@keywords statistical_test
#'@section G0 Identification:
#'@examples
#'cycle_mean_scores <- get_score(t(test_exp))$mean_score
#'kmean_classification <- kmeans_clustering(cycle_mean_scores)
#'statistical_test(cycle_mean_scores, kmean_classification, threshold = 0.001)
#'# Sample output
#'#"Possible group of G0-like cells is: 2"
#'

statistical_test <- function(mean_scores, kmean_result, threshold = 0.001){
  # find cluster with lowest mean
  require(car)
  require(dplyr)
  mean_val <- matrix(0, nrow = length(kmean_result), ncol = 6)
  for(i in 1:length(kmean_result)){
    mean_v <- colMeans(mean_scores[which(kmean_result == i),])
    mean_val[i, ] <- mean_v
  }
  #print(which.min(rowMeans(mean_val)))
  lowest_cluster <- which.min(rowMeans(mean_val))
  #return(mean_val)
  
  new_kmean_result <- kmean_result[-which(kmean_result == lowest_cluster)]
  comparisons <- unique(new_kmean_result)
  compare_group <- which(kmean_result == lowest_cluster)
  c_group_length <- length(compare_group)
  n_total <- nrow(mean_scores)
  new_data_table <- group_by(mean_scores, as.factor(kmean_result))
  colnames(new_data_table)[7] <- "cls_result"
  
  store_results <- matrix(0, nrow = length(comparisons), ncol = 6)
  store_levine <- matrix(0, nrow = length(comparisons), ncol = 6)
  store_oneway <- matrix(0, nrow = length(comparisons), ncol = 6)
  store_t_test <- matrix(0, nrow = length(comparisons), ncol = 6)
  store_nonparam <- matrix(0, nrow = length(comparisons), ncol = 6)
  for(i in 1:length(comparisons)){
    for(j in 1:6){
      tmp_group <- c(which(kmean_result == comparisons[i]), compare_group)
      comp_group_length <- length(which(kmean_result == comparisons[i]))
      tmp_data <- mean_scores[tmp_group, j]
      tmp_group_label <- factor(c(rep("comp", comp_group_length), rep("control", c_group_length)))
      model <- lm(tmp_data~tmp_group_label)
      levene <- leveneTest(tmp_data~tmp_group_label)
      oneway_model <- oneway.test(tmp_data~tmp_group_label)
      ttest_model <- pairwise.t.test(tmp_data, tmp_group_label)
      kruskal_model <- kruskal.test(tmp_data~tmp_group_label)
      
      levene_stats <- levene$Pr[1]
      anova_stats <- anova(model)$Pr[1]
      oneway_stats <- oneway_model$p.value
      ttest_stats <- ttest_model$p.value
      kruskal_stats <- kruskal_model$p.value
      
      store_results[i, j] <- anova_stats
      store_levine[i, j] <- levene_stats
      store_oneway[i, j] <- oneway_stats
      store_t_test[i, j] <- ttest_stats
      store_nonparam[i, j] <- kruskal_stats
    }
    #print(i)
  }
  rownames(store_results) <- make.names(comparisons)
  colnames(store_results) <- colnames(mean_scores)
  #print(round(store_nonparam, 3))
  #print(store_oneway)
  if(length(which(store_results > threshold))==0){
    print(paste("Possible group of G0-like cells is: ",lowest_cluster, sep = ""))
  }else{
    print("None identified")
  }
  return(store_results)
}
######################################################################################################

#' Visualization of G0-like distribution in mean scores
#' 
#' Provides boxplot for visualizing the cell cluster distribution across different cell cycle stages mean scores
#'@param mean_scores The mean score results returned from the reCAT software, "get_score"
#'@param kmean_result The classification of kmeans clustering result.
#'@export
#'@keywords distribution_plot
#'@section G0 Identification:
#'@examples
#'cycle_mean_scores <- get_score(t(test_exp))$mean_score
#'kmean_classification <- kmeans_clustering(cycle_mean_scores)
#'require(scater)
#'distribution_plot(cycle_mean_scores, kmean_classification)

distribution_plot <- function(mean_scores, kmean_result){
  
  
  par(mfrow=c(3, 2)) 
  # Plot for each mean score
  g1 <- ggplot(mean_scores, aes(x = kmean_result, y = G1Score)) +
    geom_boxplot(fill = "grey80", colour = "blue") +
    scale_x_discrete() + xlab("Kmeans Cluster") + ggtitle("G1") + 
    ylab("Mean score")
  g1s <- ggplot(mean_scores, aes(x = kmean_result, y = G1SScore)) +
    geom_boxplot(fill = "grey80", colour = "blue") +
    scale_x_discrete() + xlab("Kmeans Cluster") + ggtitle("G1S") + 
    ylab("Mean score")
  s <- ggplot(mean_scores, aes(x = kmean_result, y = SScore)) +
    geom_boxplot(fill = "grey80", colour = "blue") +
    scale_x_discrete() + xlab("Kmeans Cluster") + ggtitle("S") + 
    ylab("Mean score")
  g2 <- ggplot(mean_scores, aes(x = kmean_result, y = G2Score)) +
    geom_boxplot(fill = "grey80", colour = "blue") +
    scale_x_discrete() + xlab("Kmeans Cluster") + ggtitle("G2") + 
    ylab("Mean score")
  g2m <- ggplot(mean_scores, aes(x = kmean_result, y = G2MScore)) +
    geom_boxplot(fill = "grey80", colour = "blue") +
    scale_x_discrete() + xlab("Kmeans Cluster") + ggtitle("G2M") + 
    ylab("Mean score")
  m <- ggplot(mean_scores, aes(x = kmean_result, y = MScore)) +
    geom_boxplot(fill = "grey80", colour = "blue") +
    scale_x_discrete() + xlab("Kmeans Cluster") + ggtitle("M") + 
    ylab("Mean score")
  
  multiplot(g1, g1s, s, g2, g2m, m, cols = 2)
}

#' Visualization of G0-like distribution in mean scores along a cell cycle time-series. 
#' 6 plots for each of G1, S, G1/S, G2, M, G2/M mean scores.
#' 
#' Provides mean score expression changees along the cell cycle time-series with classification of different cell clusters.
#'@param mean_scores The mean score results returned from the reCAT software, "get_score"
#'@param kmean_result The classification of kmeans clustering result.
#'@param ordIndex The order obtained from reCAT, "get_ordIndex" which returns the cell cycle time-series of the data.
#'@export
#'@keywords mean_score_plot
#'@section G0 Identification:
#'@examples
#'cycle_mean_scores <- get_score(t(test_exp))$mean_score
#'kmean_classification <- kmeans_clustering(cycle_mean_scores)
#'reCAT_ordIndex <- get_ordIndex(test_exp, threadnum = 2)
#'require(scater)
#'mean_score_plot(cycle_mean_scores, kmean_classification, reCAT_ordIndex)
mean_score_plot <- function(mean_scores, kmean_result, ordIndex){
  par(mfrow=c(3, 2)) 
  n_ <- nrow(mean_scores)
  g1 <- qplot(x = c(1:n_), y = mean_scores$G1Score[ordIndex], 
              col = as.factor(kmean_result)[ordIndex])+#plot_labels[,2])+
    labs(title = "G1", x = "Pseudotime", y = "Mean Scores", colour = "Selection")+
    theme(axis.title=element_text(size=10),title=element_text(size = 10),legend.text = element_text(size = 10), 
          legend.title = element_text(size= 10), axis.text = element_text(size = 10))+geom_point(size = 2)
  g1s <- qplot(x = c(1:n_), y = mean_scores$G1SScore[ordIndex], 
               col = as.factor(kmean_result)[ordIndex])+#plot_labels[,2])+
    labs(title = "G1S", x = "Pseudotime", y = "Mean Scores", colour = "Selection")+
    theme(axis.title=element_text(size=10),title=element_text(size = 10),legend.text = element_text(size = 10), 
          legend.title = element_text(size= 10), axis.text = element_text(size = 10))+geom_point(size = 2)
  s <- qplot(x = c(1:n_), y = mean_scores$SScore[ordIndex], 
             col = as.factor(kmean_result)[ordIndex])+#plot_labels[,2])+
    labs(title = "S", x = "Pseudotime", y = "Mean Scores", colour = "Selection")+
    theme(axis.title=element_text(size=10),title=element_text(size = 10),legend.text = element_text(size = 10), 
          legend.title = element_text(size= 10), axis.text = element_text(size = 10))+geom_point(size = 2)
  g2 <- qplot(x = c(1:n_), y = mean_scores$G2Score[ordIndex], 
              col = as.factor(kmean_result)[ordIndex])+#plot_labels[,2])+
    labs(title = "G2", x = "Pseudotime", y = "Mean Scores", colour = "Selection")+
    theme(axis.title=element_text(size=10),title=element_text(size = 10),legend.text = element_text(size = 10), 
          legend.title = element_text(size= 10), axis.text = element_text(size = 10))+geom_point(size = 2)
  g2m <- qplot(x = c(1:n_), y = mean_scores$G2MScore[ordIndex], 
               col = as.factor(kmean_result)[ordIndex])+#plot_labels[,2])+
    labs(title = "G2M", x = "Pseudotime", y = "Mean Scores", colour = "Selection")+
    theme(axis.title=element_text(size=10),title=element_text(size = 10),legend.text = element_text(size = 10), 
          legend.title = element_text(size= 10), axis.text = element_text(size = 10))+geom_point(size = 2)
  m <- qplot(x = c(1:n_), y = mean_scores$MScore[ordIndex], 
             col = as.factor(kmean_result)[ordIndex])+#plot_labels[,2])+
    labs(title = "M", x = "Pseudotime", y = "Mean Scores", colour = "Selection")+
    theme(axis.title=element_text(size=10),title=element_text(size = 10),legend.text = element_text(size = 10), 
          legend.title = element_text(size= 10), axis.text = element_text(size = 10))+geom_point(size = 2)
  
  multiplot(g1, g1s, s, g2, g2m, m, cols = 2)
}