# Post-Biological Analysis: Plotting Functions

##############################################################################################
# REDPATH COLORING SCHEMES



redpath_colorset <- c("#FF6766", "#97CE68", "#FF9900", "#6BCBCA", "#F8C765", 
                      "#f2f28b", "#FA723E", "#9F80DF", "#554E44", "#F58F8F")
#hmm_colorset <- c("#F2F28B", "#F8C765", "#FA723E")

hmm_colorset <- c("#F2F28B", "#F8C765", "#f59a1c")


redpath_cols <- function(...){
  cols <- c(...)
  
  if(is.null(cols))
    return(redpath_colorset)
  
  redpath_colorset[cols]
}

redpath_palettes <- list(
  'main' = redpath_cols("light_red", "sl_desat_green", "pure_orange"),
  'all' = redpath_cols()
)
redpath_palettes_func <- function(no_col = 3){
  return(redpath_colorset[1:no_col])
}


redpath_theme <- function(){
  theme(
    plot.title = element_text(size = 20),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.key = element_rect(fill = "white"))
}

##########################################################################################################
require(plotly)
require(dplyr)
require(GOsummaries)
require(gplots)


reassign_hmm <- function(hmm_results, sorted_exprs_val){
  n_cls <- length(unique(hmm_results[1,]))
  new_hmm_results <- hmm_results
  
  
  
  for(i in c(1:nrow(sorted_exprs_val))){
    if(n_cls == 2){
      l1 <- which(hmm_results[i, ] == 1)
      l2 <- which(hmm_results[i, ] == 2)
      if(length(l1) == 0){
        c1 = 0
      }else{
        c1 <- mean(as.numeric(sorted_exprs_val[i, which(hmm_results[i, ] == 1)]))
      }
      if(length(l2) == 0){
        c2 = 0
      }else{
        c2 <- mean(as.numeric(sorted_exprs_val[i, which(hmm_results[i, ] == 1)]))
      }

      val_vec <- c(c1, c2)
      val_order <- order(val_vec)
      char_vec <- c("Low", "High")
      if(length(l1) == 0){
      }else{
        new_hmm_results[i, which(hmm_results[i, ]==1)] <- char_vec[val_order[1]]
      }
      if(length(l1) == 0){
      }else{
        new_hmm_results[i, which(hmm_results[i, ]==2)] <- char_vec[val_order[2]]#val_order[2]
      }

      
      #new_hmm_results[i, ] <- (n_cls+1) - hmm_results[i, ]
      
    }else{
      c1 <- mean(as.numeric(sorted_exprs_val[i, which(hmm_results[i, ] == 1)]))
      c2 <- mean(as.numeric(sorted_exprs_val[i, which(hmm_results[i, ] == 2)]))
      c3 <- mean(as.numeric(sorted_exprs_val[i, which(hmm_results[i, ] == 3)]))
      val_vec <- c(c1, c2, c3)
      val_order <- order(val_vec)
      char_vec <- c("Low", "Medium", "High")
      
      new_hmm_results[i, which(hmm_results[i, ]==1)] <- char_vec[val_order[1]]#val_order[1]
      new_hmm_results[i, which(hmm_results[i, ]==2)] <- char_vec[val_order[2]]#val_order[2]
      new_hmm_results[i, which(hmm_results[i, ]==3)] <- char_vec[val_order[3]]#val_order[3]
    }
  }
  return(new_hmm_results)
  
}
#' Plot function for HMM & gene expression changes
#' 
#' Various versions of plots are available here to visualize the most significant gene expression change from dCor and MIC. 
#' 
#' 
#'@param sorted_exprs_val The gene expression input obtained from gene_selection. The ranked and sorted gene_expression matrix.
#'@param pseudotime Pseudo development time
#'@param sorted_hmm_results Defaults to NULL. This is the inferred HMM matrix obtained from calc_all_hmm_diff_parallel.
#'@param color_label Defaults to NULL. The cell type labels.
#'@param num_plots Number of plots desired, ranking from the most significant.
#'@export
#'@keywords plot_hmm_full
#'@section Biological Analysis:
#'@examples
#'selected_result <- gene_selection(full_exprs_matrix, redpath_pseudotime, dcor_result, mic_result, threshold = 0.5)
#'gene_state_result <- calc_all_hmm_diff_parallel(selected_result, cls_num = 2, threadnum = 4)
#'hmm_plots <- plot_hmm_full(selected_result, redpath_pseudotime, gene_state_result, cell_type_label, num_plots = 10)
#'require(scater)
#'multiplot(plotlist = hmm_plots, layout = matrix(c(1:10), 5, 2)) # For a total number of 10 plots
plot_hmm_full <- function(sorted_exprs_val, pseudotime, sorted_hmm_results=NULL, color_label=NULL, num_plots=10){
  output_list <- list()

  n_results <- nrow(sorted_exprs_val)
  n_cells <- ncol(sorted_exprs_val)
  
  #sorted_hmm_results <- sorted_hmm_results[c(n_results:(n_results - num_plots + 1)), ]
  sorted_exprs_val <- sorted_exprs_val[c(n_results:(n_results - num_plots + 1)), ]
  if(is.null(sorted_hmm_results) & !is.null(color_label)){
    # Plot with cell type only
    for(i in c(1:num_plots)){
      n_col <- length(unique(color_label))
      headings <- rownames(sorted_exprs_val)[i]
      tmp_data <- data.frame(X = pseudotime[order(pseudotime)],#c(1:n_cells),
                             Y = as.numeric(sorted_exprs_val[i, ]),
                             col_label = as.factor(color_label)[order(pseudotime)])
      g <- ggplot(tmp_data, aes(X, Y, color = as.factor(col_label)))+geom_point(size = 3)+
        scale_color_manual(name = "", values = as.character(redpath_colorset[1:n_col])) +
        xlab("Differential Time") + ylab("Logged Exprs") + labs(colour = "Cell Type", title = headings)+
        
        redpath_theme()
      
      output_list[[i]] <- g
    }
  }else if (!is.null(sorted_hmm_results) & !is.null(color_label)){
    # Plot with HMM colorbar & Cell Type coloring
    #############################################

    n_col <- length(unique(color_label))
    n_hmm <- length(unique(sorted_hmm_results[1, ]))
    sorted_hmm_results <- sorted_hmm_results[c(n_results:(n_results - num_plots + 1)), ]
    sorted_hmm_results <- reassign_hmm(sorted_hmm_results, sorted_exprs_val)
    
    for(i in c(1:num_plots)){
      headings <- rownames(sorted_exprs_val)[i]
      col <-sorted_hmm_results[i,]
      tmp_data <- data.frame(X = pseudotime[order(pseudotime)],#c(1:n_cells),
                             Y = as.numeric(sorted_exprs_val[i, ]),
                             type_label = as.factor(color_label)[order(pseudotime)])
      tmp_data2 <- data.frame(x = pseudotime[order(pseudotime)], y = max(sorted_exprs_val[i, ])+5, color = col)
      g <- ggplot(tmp_data, aes(X, Y, color = type_label))+geom_point(size = 3)+ylim(c(0, max(tmp_data$Y)+2.5))+
        #scale_fill_brewer(palette=custom_color_set)+#(palette = "Set2")+
        scale_color_manual(name = "", labels = c(as.character(unique(col)), as.character(unique(color_label[order(pseudotime)]))),
                           values = c(as.character(hmm_colorset[1:n_hmm]), as.character(redpath_colorset[1:n_col]))) +
        #scale_color_manual(name = "", values = c(as.character(col[1:n_hmm]), as.character(redpath_colorset[1:n_col]))) +
        
        #scale_color_manual(name = "", values = as.character(hmm_colorset)) +
        xlab("Differential Time") + ylab("Logged Exprs") + labs(colour = "HMM segmentation", title = headings)+
        
        redpath_theme()
      if(n_hmm == 3){
        g <- g + 
          geom_segment(data = tmp_data2[which(col=="Low"), ], aes(x = x, xend = x, y = y, yend = y+2, #color = as.factor(col),
                                                                  colour = as.character(hmm_colorset[1])), size = 1.5, show.legend = F)+
          geom_segment(data = tmp_data2[which(col =="Medium"), ], aes(x = x, xend = x, y = y, yend = y+2, #color = as.factor(col),
                                                                      colour = as.character(hmm_colorset[2])), size = 1.5, show.legend = F)+
          
          geom_segment(data = tmp_data2[which(col =="High"), ], aes(x = x, xend = x, y = y, yend = y+2, #color = as.factor(col),
                                                                    colour = as.character(hmm_colorset[3])), size = 1.5, show.legend = F)
      }else if(n_hmm == 2){
        g <- g+
          geom_segment(data = tmp_data2[which(col=="Low"), ], aes(x = x, xend = x, y = y, yend = y+2, #color = as.factor(col),
                                                                  colour = as.character(hmm_colorset[1])), size = 1.5, show.legend = F)+
          geom_segment(data = tmp_data2[which(col =="High"), ], aes(x = x, xend = x, y = y, yend = y+2, #color = as.factor(col),
                                                                    colour = as.character(hmm_colorset[3])), size = 1.5, show.legend = F)


          
      }
      output_list[[i]] <- g
    }
  }else if (!is.null(sorted_hmm_results) & is.null(color_label)) {
    # Plot with HMM only

    # Define low / med / high
    ###############################################
    n_hmm <- length(unique(sorted_hmm_results[1,]))

    sorted_hmm_results <- sorted_hmm_results[c(n_results:(n_results - num_plots + 1)), ]
    sorted_hmm_results <- reassign_hmm(sorted_hmm_results, sorted_exprs_val)
    
    ###############################################
    for(i in c(1:num_plots)){
      headings <- rownames(sorted_exprs_val)[i]
      col <- factor(sorted_hmm_results[i, ], levels = c("Low", "Medium", "High"))
      tmp_data <- data.frame(X = pseudotime[order(pseudotime)],#c(1:n_cells),
                             Y = as.numeric(sorted_exprs_val[i, ]),
                             HMM_label = sorted_hmm_results[i, ])
      tmp_data2 <- data.frame(x = pseudotime[order(pseudotime)], y = max(sorted_exprs_val[i, ])+5, color = col)
      g <- ggplot(tmp_data, aes(X, Y, color = as.factor(HMM_label)))+geom_point(size = 2)+
        #scale_fill_brewer(palette=custom_color_set)+#(palette = "Set2")+
        #geom_segment(data = tmp_data2, aes(x = x, xend = x, y = y, yend = y+2, color = col), size = 3)+
        scale_color_manual(name = "", values = as.character(hmm_colorset[1:n_hmm])) +
        xlab("Differential Time") + ylab("Logged Exprs") + labs(colour = "HMM segmentation", title = headings)+
        
        redpath_theme()
      output_list[[i]] <- g
    }
  }else{
    # Plot top ranked gene exprs only
    for(i in c(1:num_plots)){
      headings <- rownames(sorted_exprs_val)[i]
      tmp_data <- data.frame(X = pseudotime[order(pseudotime)],#c(1:n_cells),
                             Y = as.numeric(sorted_exprs_val[i, ]))#,
                             #col_label = as.factor(color_label))
      g <- ggplot(tmp_data, aes(X, Y))+geom_point(size = 3)+
        #scale_color_manual(name = "", values = as.character(redpath_colorset[1:3])) +
        xlab("Differential Time") + ylab("Logged Exprs") + labs(title = headings)+
        
        redpath_theme()
      output_list[[i]] <- g
    }
  }

  return(output_list)
}

plot_hmm <- function(sorted_hmm_results, sorted_exprs_val, pseudotime, num_plots){
  output_list <- list()

  sorted_hmm_results <- reassign_hmm(sorted_hmm_results, sorted_exprs_val)
  n_results <- nrow(sorted_hmm_results)
  n_cells <- ncol(sorted_hmm_results)
  
  sorted_hmm_results <- sorted_hmm_results[c(n_results:(n_results - num_plots + 1)), ]
  sorted_exprs_val <- sorted_exprs_val[c(n_results:(n_results - num_plots + 1)), ]
  for(i in c(1:num_plots)){
    headings <- rownames(sorted_exprs_val)[i]
    col <-sorted_hmm_results[i,]
    tmp_data <- data.frame(X = pseudotime[order(pseudotime)],#c(1:n_cells),
                           Y = as.numeric(sorted_exprs_val[i, ]),
                           HMM_label = sorted_hmm_results[i, ])
    tmp_data2 <- data.frame(x = pseudotime[order(pseudotime)], y = max(sorted_exprs_val[i, ])+5, color = col)
    g <- ggplot(tmp_data, aes(X, Y, color = factor(HMM_label, levels = c("Low", "Medium", "High"))))+geom_point(size = 4)+
      #scale_fill_brewer(palette=custom_color_set)+#(palette = "Set2")+
      scale_color_manual(name = "", values = c(as.character(hmm_colorset), as.character(redpath_colorset[1:3]))) +
#scale_color_manual(name = "", values = as.character(hmm_colorset)) +
      xlab("Differential Time") + ylab("logged exprs") + labs(colour = "HMM segmentation", title = headings)+
      redpath_theme()
    g <- g + 
      geom_segment(data = tmp_data2[which(col=="Low"), ], aes(x = x, xend = x, y = y, yend = y+2, #color = as.factor(col),
                                                              colour = as.character(hmm_colorset[1])), size = 1.5)+
      geom_segment(data = tmp_data2[which(col =="Medium"), ], aes(x = x, xend = x, y = y, yend = y+2, #color = as.factor(col),
                                                                  colour = as.character(hmm_colorset[2])), size = 1.5)+
      
      geom_segment(data = tmp_data2[which(col =="High"), ], aes(x = x, xend = x, y = y, yend = y+2, #color = as.factor(col),
                                                                colour = as.character(hmm_colorset[3])), size = 1.5)

    output_list[[i]] <- g
  }
  return(output_list)
}


#' Plot function for specific gene expression changes
#' 
#' Plotting function for specific genes of interest, along the pseudo development time. 
#' 
#' 
#'@param selected_gene The enquired specific gene expression matrix, obtained from "get_specific_gene_exprs".
#'@param pseudotime Pseudo development time
#'@param color_label The cell type labels.
#'@param order If order = T, gene expression change will be plotted against a uniform pseudo time. If order F, it will be plotted against the actual value of the calculated pseudo time.
#'@export
#'@keywords plot_gene_exprs
#'@section Biological Analysis:
#'@examples
#'cdk_gene_exprs <- get_specific_gene_exprs(full_gene_exprs_data, gene_name = "cdk11", type = "grep", organism = "mm")
#'plot_cdk <- plot_gene_exprs(cdk_gene_exprs, redpath_pseudotime, color_label = cell_type_label, order = T)
#'require(scater)
#'multiplot(plotlist = plot_cdk, layout = matrix(c(1:2), 2, 1)) # For a total number of 2 plots

plot_gene_exprs <- function(selected_gene, pseudotime, color_label, order = T){
  output_list <- list()
  if(is.null(dim(selected_gene$exprs))){
    gene_cell_val <- data.frame(t(selected_gene$exprs))
  }else{
    gene_cell_val <- data.frame(selected_gene$exprs)
  }

  num_plots <- nrow(gene_cell_val)
  #gene_cell_val <- data.frame(gene_cell_val)
  n_cells <- ncol(gene_cell_val)
  # Plot with cell type only
  for(i in c(1:num_plots)){
    n_col <- length(unique(color_label))
    headings <- selected_gene$namesr#rownames(gene_cell_val)[i]
    if(order == F){
      tmp_data <- data.frame(X = pseudotime,#c(1:n_cells),
                             Y = as.numeric(gene_cell_val[i, ]),
                             col_label = as.factor(color_label))
    }else if(order ==T){
      tmp_data <- data.frame(X = c(1:n_cells),
                             Y = as.numeric(gene_cell_val[i, ])[order(pseudotime)],
                             col_label = as.factor(color_label)[order(pseudotime)])
    }
    g <- ggplot(tmp_data, aes(X, Y, color = as.factor(col_label)))+geom_point(size = 3)+
      scale_color_manual(name = "", values = as.character(redpath_colorset[1:n_col])) +
      xlab("Differential Time") + ylab("Logged Exprs") + labs(colour = "Cell Type", title = headings)+
      
      redpath_theme()
    
    output_list[[i]] <- g
  }
  # for(i in c(1:num_plots)){
  #   headings <- rownames(sorted_exprs_val)[i]
  #   tmp_data <- data.frame(X = pseudotime[order(pseudotime)],#c(1:n_cells),
  #                          Y = as.numeric(sorted_exprs_val[i, ]),
  #                          col_label = as.factor(color_label))
  #   g <- ggplot(tmp_data, aes(X, Y, color = as.factor(col_label)))+geom_point(size = 2)+
  #     scale_color_manual(name = "", values = as.character(redpath_colorset[1:3])) +
  #     xlab("Differential Time") + ylab("Logged Exprs") + labs(colour = "Cell Type", title = headings)
  #   output_list[[i]] <- g
  # }
  return(output_list)
}


#' Gene Clustering
#' 
#' Simple gene cluster on the significant genes sorted matrix (based on dCor and MIC). Can be either HMM matrix or the gene expression matrix
#' 
#'@param sorted_matrix This can be either the HMM matrix or the sorted gene selection matrix obtained from earlier functions.
#'@param cls_num The number of gene cluster desired
#'@export
#'@keywords get_gene_cluster
#'@section Biological Analysis:
#'@examples
#'selected_result <- gene_selection(full_exprs_matrix, redpath_pseudotime, dcor_result, mic_result, threshold = 0.5)
#'gene_state_result <- calc_all_hmm_diff_parallel(selected_result, cls_num = 2, threadnum = 4) 
#'inferred_gene_cluster <- get_gene_cluster(sorted_matrix = gene_state_result, cls_num = 2)
#'heatmap_plot <- plot_heatmap_analysis(sorted_matrix = inferred_gene_cluster, redpath_pseudotime, cell_labels = NULL, gene_cluster = inferred_gene_cluster, input_type = "hmm")

get_gene_cluster <- function(sorted_matrix, cls_num = 8){
  h_clust <- hclust(dist(sorted_matrix))
  h_gene_cluster <- cutree(h_clust, cls_num)
  
  return(h_gene_cluster)
}

#' Infer Gene Ontology summaries
#' 
#' Relying on GOsummaries package, word cloud plots are generated using this function.
#' 
#'@param gene_cluster Gene cluster classification, with gene names.
#'@param organism Type of species of the genes. Default: "hsapiens"
#'@export
#'@keywords infer_GO_summaries
#'@section Biological Analysis:
#'@examples
#'selected_result <- gene_selection(full_exprs_matrix, redpath_pseudotime, dcor_result, mic_result, threshold = 0.5)
#'gene_state_result <- calc_all_hmm_diff_parallel(selected_result, cls_num = 2, threadnum = 4) 
#'inferred_gene_cluster <- get_gene_cluster(sorted_matrix = gene_state_result, cls_num = 2)
#'GO_summary_result <- infer_GO_summaries(inferred_gene_cluster, organism = "hsapiens")
#'plot(GO_summary_result, fontsize = 8)

infer_GO_summaries <- function(gene_cluster, organism="hsapiens"){
  n_cls <- length(unique(gene_cluster))
  tmp_geneset <- list()
  name_list <- c()
  for(i in c(1:n_cls)){
    g_set <- names(gene_cluster)[which(gene_cluster == i)]
    tmp_geneset[i] <- list(g_set)
    name_list <- c(name_list, paste("C", i, sep = ""))
  }

  g1_list <- list(tmp_geneset)
  #names(g1_list) <- toupper(name_list)
  print(str(g1_list[[1]]))
  #pritn("1")
  go_summary_results <- gosummaries(g1_list[[1]], organism = organism)
  plot(go_summary_results, fontsize = 8)
  return(go_summary_results)
}

gene_cluster_colors <- 
  c(
    "1" = "#97CE68",
    "2" = "#F2F28B", #554E44
    "3" = "#F8C765",
    "4" = "#FA723E",
    "5" = "#554E44", #FF9900
    "6" = "#9F80DF",
    "7" = "#F58F8F",
    "8" = "#FF9900" #F2F28B #554E44
  )
#' Heatmap Analysis
#' 
#' Comprehensive heatmap plot function with gene cluster and cell type labels. 
#' 
#'@param sorted_matrix This can be either the HMM matrix or the sorted gene selection matrix obtained from earlier functions.
#'@param pseudo_time Pseudo development time
#'@param cell_labels Defaults to NULL. The cell type labels.
#'@param gene_cluster Defaults to NULL. It can be automatically caluclated with defualt gene cluster of 8. Or user can define the inferred gene cluster results using hierarchical clustering (get_gene_cluster) or by other means. Note that the input must match for the input matrix (ie gene / hmm)
#'@param heading Plot title
#'@param input_type "hmm" or "gene". Please specify type of sorted matrix input to the function. 
#'@export
#'@keywords plot_heatmap_analysis
#'@section Biological Analysis:
#'@examples
#'selected_result <- gene_selection(full_exprs_matrix, redpath_pseudotime, dcor_result, mic_result, threshold = 0.5)
#'gene_state_result <- calc_all_hmm_diff_parallel(selected_result, cls_num = 2, threadnum = 4) 
#'inferred_gene_cluster <- get_gene_cluster(sorted_matrix = gene_state_result, cls_num = 8)
#'heatmap_plot <- plot_heatmap_analysis(sorted_matrix = inferred_gene_cluster, redpath_pseudotime, cell_labels = NULL, gene_cluster = inferred_gene_cluster, input_type = "hmm")

plot_heatmap_analysis <- function(sorted_matrix, pseudo_time, cell_labels = NULL, gene_cluster=NULL, heading = "HEATMAP", input_type = "hmm"){
  
  if(is.null(gene_cluster)){
    gene_cluster <- get_gene_cluster(sorted_matrix)
  }
  
  gene_cluster_col <- as.character(gene_cluster_colors[as.character(gene_cluster)])
  cell_unique_lab <- unique(cell_labels)
  cell_lab_col <- redpath_cols()[c(1:length(cell_unique_lab))]
  names(cell_lab_col) <- cell_unique_lab

  cell_label_col <- as.character(cell_lab_col[as.character(cell_labels)])
  if(input_type == "hmm"){
    if(is.null(cell_labels)){
      h1 <- 
        heatmap.2(as.matrix(sorted_matrix), scale="none", col=c("#F2F2F2BF", "#FF0000BF"), 
                  #ColSideColors = cell_labels[order(pseudo_time))], 
                  RowSideColors = gene_cluster_col,  
                  key=TRUE, keysize = 1.4, symkey=F, density.info="none", trace="none", offsetRow = 0, offsetCol = 0, cexCol = 0.1, cexRow=0.1, 
                  dendrogram = "none", Colv = F, labRow = FALSE, labCol = FALSE, 
                  main = heading, xlab = "Pseudo-Time", ylab = "Genes")
    }else{
      h1 <- 
        heatmap.2(as.matrix(sorted_matrix), scale="none", col =c("#F2F2F2BF", "#FF0000BF"), 
                ColSideColors = cell_label_col[order(pseudo_time)], 
                RowSideColors = gene_cluster_col,  
                key=TRUE, keysize = 1.4, symkey=F, density.info="none", trace="none", offsetRow = 0, offsetCol = 0, cexCol = 0.1, cexRow=0.1, 
                dendrogram = "none", Colv = F, labRow = FALSE, labCol = FALSE, 
                main = heading, xlab = "Pseudo-Time", ylab = "Genes")
    }
  }else if(input_type == "gene"){
    if(is.null(cell_labels)){
      h1 <- 
        heatmap.2(as.matrix(sorted_matrix), scale="none", col = rev(heat.colors(n = 12, alpha = 1)), 
                  #ColSideColors = cell_labels[order(pseudo_time))], 
                  RowSideColors = gene_cluster_col, #as.character(mgh107_gene_cluster), 
                  key=TRUE, keysize = 1.4, symkey=F, density.info="none", trace="none", offsetRow = 0, offsetCol = 0, cexCol = 0.1, cexRow=0.1, 
                  dendrogram = "none", Colv = F, labRow = FALSE, labCol = FALSE, 
                  main = heading, xlab = "Pseudo-Time", ylab = "Genes")
    }else{
      h1 <- 
        heatmap.2(as.matrix(sorted_matrix), scale="none", col = rev(heat.colors(n = 12, alpha = 1)), 
                  ColSideColors = cell_label_col[order(pseudo_time)], 
                  RowSideColors = gene_cluster_col, #as.character(mgh107_gene_cluster), 
                  key=TRUE, keysize = 1.4, symkey=F, density.info="none", trace="none", offsetRow = 0, offsetCol = 0, cexCol = 0.1, cexRow=0.1, 
                  dendrogram = "none", Colv = F, labRow = FALSE, labCol = FALSE, 
                  main = heading, xlab = "Pseudo-Time", ylab = "Genes")
    }

    
    
  }
  return(h1)

}

# g1 <- plot_hmm(mimic_hmm[sort_mimic, ], mimic_exprs[sort_mimic, ], 10)
# par(mfrow = c(5, 2))

# multiplot(plotlist = g1, layout = matrix(c(1:10), ncol = 2, nrow = 5))
# g2 <- plot_hmm(mydcor_hmm[sort_mydcor, ], mydcor_exprs[sort_mydcor, ], 10)
# multiplot(plotlist = g2, layout = matrix(c(1:10), ncol = 2, nrow = 5))

# #terrain.colors(n = 2, alpha = 0.75) #00A600FF" "#F2F2F2FF"
# heatmap.2(mimic_hmm_1000[,order(schlitz_test15_2)], col=c("#F2F2F2BF", "#FF0000BF"), scale="none", ColSideColors=schlitz_int[order(schlitz_test15_2)],
#           #RowSideColors = as.character(gene_cluster), 
#           key=TRUE, keysize = 1.4, symkey=F, density.info="none", trace="none", cexRow=0.5, dendrogram = "row")
# 


# infer_GO_genesets <- function(hmm_exprs, cluster_number = 6){
#   hc_clust <- cutree(hclust(dist(hmm_exprs)), cluster_number)
#   
# }

# 3D PLOTS:



#' Search for specific gene expression
#' 
#' Gets the subset of genes which resemble the enquired gene name within the dataset. 
#' This prepares for the plotting function of multiple gene expressions along the pseudo development time. 
#' 
#'@param full_gene_exprs_matrix This is the full gene (row) by cell (col) matrix of your dataset.
#'@param gene_name The gene name that is desired / enquired. 
#'@param type pattern matching "grep" or exact "match"
#'@param organism "mm" for Mouse or "hs" for Human
#'@export
#'@keywords get_specific_gene_exprs
#'@section Biological Analysis:
#'@examples
#'cdk_gene_exprs <- get_specific_gene_exprs(full_gene_exprs_data, gene_name = "cdk11", type = "grep", organism = "mm")
#'plot_cdk <- plot_gene_exprs(cdk_gene_exprs, redpath_pseudotime, color_label = cell_type_label, order = T)
#'require(scater)
#'multiplot(plotlist = plot_cdk, layout = matrix(c(1:2), 2, 1)) # For a total number of 2 plots
get_specific_gene_exprs <- function(full_gene_exprs_matrix, gene_name = "cdk", type = "grep", organism = "mm"){
  gene_names <- (rownames(full_gene_exprs_matrix))
  if(length(grep("ens", gene_names[1]))!=0){
    if(organism == "mm"){
      gene_names <- ensembl_to_mgi(gene_names, mm_IDs)
    }else if(organism == "hs"){
      gene_names <- ensembl_to_mgi(gene_names, hs_IDs)
    }
    
  }
  if(type == "grep"){
    grep_gene_pos <- grep(tolower(gene_name), tolower(gene_names))
  }else if(type == "match"){
    grep_gene_pos <- match(tolower(gene_name), tolower(gene_names))
  }

  if(length(grep_gene_pos) == 0){
    print("Inquired gene not found in this dataset")
    
  }else{
    selected_exprs <- full_gene_exprs_matrix[grep_gene_pos, ]
  }
  range_data <- range(full_gene_exprs_matrix)
  return(list(range = range_data, exprs = selected_exprs, names = gene_names[grep_gene_pos]))
}
calc_unit_circle <- function(recat_input){
  meanRes <- colMeans(matrix(complex(argument = recat_input[1:nrow(recat_input), ]*2*pi), nrow = nrow(recat_input)))
  meanRes <- Arg(meanRes) / 2 / pi
  tmp <- which(meanRes < 0)
  meanRes[tmp] <- meanRes[tmp] + 1
  
  rev_complex <- complex(argument = meanRes * 2 * pi)
  return(rev_complex)
}



#' Coupling cell cycle and differentiation 3-D plots
#' 
#' This function incorporates output from reCAT and redPATH to illustrate some relationship between cell proliferation & differentiation.
#' 
#'@param recat_input This is the normalizedResultLst obtained from reCAT software. It's generally save in a .RData file started with bestEnsemble...
#'@param redpath_pseudotime Pseudo development time for differentiation
#'@param recat_symbol Cell cycle labels
#'@param color_label Cell type labels
#'@export
#'@keywords plot_cycle_diff_3d
#'@section Biological Analysis:
#'@examples
#'p1 <- plot_cycle_diff_3d(normalizedResultLst, redpath_pseudotime, cell_cycle_labels, cell_type_labels)
#'p1
plot_cycle_diff_3d <- function (recat_input, redpath_pseudotime, recat_symbol, color_label) 
{
  rev_complex <- calc_unit_circle(recat_input)

    p <- plot_ly(                  x = Re(rev_complex), y = Im(rev_complex), 
                                   z = 2*redpath_pseudotime,
                                   color = ~color_label,
                                   #legendgroup = ~color_label,
                                   #marker = list(size = ~(3*(marker_gene_exprs))),
                                   #marker = list(size = ~(rescale(marker_gene_exprs, to=c(11,28)))),
                                   #size = ~log(marker_gene_exprs),
                                   colors = redpath_colorset[1:3]) %>%
      
      add_markers(type = "scatter3d", 
                  x = Re(rev_complex), y = Im(rev_complex), 
                  z = 2*redpath_pseudotime,
                  symbol = ~recat_symbol,
                  
                  mode = 'markers',
                  text = ~paste("Cell Type: ", color_label),
                  #legendgroup=~recat_symbol, 
                  
                  #symbols = c("circle", "cross", "diamond"),
                  symbols = c('circle', 'diamond-open', 'circle-open'),
                  #color = ~(color_label),
                  color = ~(color_label), colors = redpath_colorset[1:3]) %>%
      # Add redpath
      #add_paths(x = Re(rev_complex)[order(redpath_pseudotime)], y = Im(rev_complex)[order(redpath_pseudotime)],
      #          z = 2*redpath_pseudotime[order(redpath_pseudotime)], inherit=F) %>%
      
  
    
    layout(title = gene_name,#"Cell Type Label", 
           scene = list(xaxis = list(title = "Real (Unit Circle of reCAT)"), 
                        yaxis = list(title = "Imaginary (Unit Circle of reCAT"), 
                        zaxis = list(title = "Differential Time")))#,
    #legend = list(tracegroupgap =30, y=0.9, yanchor="top"))
  


  return(p)
  
}


#' Coupling cell cycle and differentiation 3-D plots - Gene Exprs
#' 
#' This function incorporates output from reCAT and redPATH to illustrate some relationship between cell proliferation & differentiation. 
#' Observing the gradual change of specific marker genes
#' 
#'@param recat_input This is the normalizedResultLst obtained from reCAT software. It's generally save in a .RData file started with bestEnsemble...
#'@param redpath_pseudotime Pseudo development time for differentiation
#'@param color_label Cell type labels
#'@param marker_gene_exprs Marker genes matrix (genes of interest) obtained from "get_specific_gene_exprs"
#'@export
#'@keywords plot_diff_3d
#'@section Biological Analysis:
#'@examples
#'p1 <- plot_diff_3d(normalizedResultLst, redpath_pseudotime, cell_cycle_labels, cdk_gene_exprs)
#'p1
plot_diff_3d <- function (recat_input, redpath_pseudotime, color_label = NULL, 
                          marker_gene_exprs = NULL) 
{
  rev_complex <- calc_unit_circle(recat_input)
  if (is.null(marker_gene_exprs)) {
    p <- plot_ly(x = Re(rev_complex) + redpath_pseudotime, 
                 y = Im(rev_complex) + redpath_pseudotime, z = redpath_pseudotime, 
                 color = as.factor(color_label), colors = "Set1", 
                 text = ~paste("Cell Type: ", color_label)) %>% layout(title = "Cell Type Label", 
                                                                       scene = list(xaxis = list(title = "Real (Unit Circle of reCAT)"), 
                                                                                    yaxis = list(title = "Imaginary (Unit Circle of reCAT"), 
                                                                                    zaxis = list(title = "Differential Time")))
  }
  else {
    Gene_Expression <- marker_gene_exprs$exprs
    if (!is.null(nrow(marker_gene_exprs$exprs))) {
      store_p <- list()
      for (i in c(1:nrow(marker_gene_exprs$exprs))) {
        p <- plot_ly(x = Re(rev_complex), y = Im(rev_complex), 
                     z = redpath_pseudotime, text = ~paste("Cell Type: ", 
                                                           as.character(color_label)), marker = list(color = ~Gene_Expression[i, 
                                                                                                                              ], showscale = T, colorscale = "Reds", cauto = F, 
                                                                                                     cmin = marker_gene_exprs$range[1], cmax = marker_gene_exprs$range[2])) %>% 
          layout(title = marker_gene_exprs$names[i], 
                 scene = list(xaxis = list(title = "Real (Unit Circle of reCAT)"), 
                              yaxis = list(title = "Imaginary (Unit Circle of reCAT"), 
                              zaxis = list(title = "Differential Time")))
        store_p[[i]] <- p
      }
      return(store_p)
    }
    else {
      p <- plot_ly(x = Re(rev_complex), y = Im(rev_complex), 
                   z = redpath_pseudotime, text = ~paste("Cell Type: ", 
                                                         as.character(color_label)), marker = list(color = ~Gene_Expression, 
                                                                                                   showscale = T, colorscale = "Reds", cauto = F, 
                                                                                                   cmin = marker_gene_exprs$range[1], cmax = marker_gene_exprs$range[2])) %>% 
        layout(title = marker_gene_exprs$names, scene = list(xaxis = list(title = "Real (Unit Circle of reCAT)"), 
                                                             yaxis = list(title = "Imaginary (Unit Circle of reCAT"), 
                                                             zaxis = list(title = "Differential Time")))
      show(p)
      return(p)
    }
  }
}

# Deprecated:
# plot_cycle_diff_3d <- function(recat_input, redpath_pseudotime, color_label = NULL, marker_gene_exprs = NULL){
#   
#   rev_complex <- calc_unit_circle(recat_input)
#   if(is.null(marker_gene_exprs)){
#     p <- plot_ly(x = Re(rev_complex), y = Im(rev_complex), z = redpath_pseudotime,
#                  color = as.factor(color_label), colors = "Set1",
#                  text = ~paste("Cell Type: ", color_label)) %>%
#       layout(title = "Cell Type Label",
#              scene = list(xaxis = list(title = 'Real (Unit Circle of reCAT)'),
#                           yaxis = list(title = 'Imaginary (Unit Circle of reCAT'),
#                           zaxis = list(title = 'Differential Time')))
#   }else{
#     Gene_Expression <- marker_gene_exprs$exprs
#     if(!is.null(nrow(marker_gene_exprs$exprs))){
#       store_p <- list()
#       for(i in c(1:nrow(marker_gene_exprs$exprs))){
#         p <- plot_ly(x = Re(rev_complex), y = Im(rev_complex), z = redpath_pseudotime,
#                               text = ~paste("Cell Type: ", as.character(color_label)),
#                               marker = list(color = ~Gene_Expression[i, ],  showscale = T,
#                                             #colors = "Reds", 
#                                             colorscale = "Reds",
#                                             cauto = F, cmin = marker_gene_exprs$range[1], cmax = marker_gene_exprs$range[2])) %>%
#                               layout(title = marker_gene_exprs$names[i],
#                                      scene = list(xaxis = list(title = 'Real (Unit Circle of reCAT)'),
#                                                   yaxis = list(title = 'Imaginary (Unit Circle of reCAT'),
#                                                   zaxis = list(title = 'Differential Time')))#,
#         store_p[[i]] <- p
#       }
#       return(store_p)
#     }else{
#       p <- plot_ly(x = Re(rev_complex), y = Im(rev_complex), z = redpath_pseudotime,
#                    text = ~paste("Cell Type: ", as.character(color_label)),
#                    marker = list(color = ~Gene_Expression,  showscale = T,
#                                  #colors = "Reds", 
#                                  colorscale = "Reds",
#                                  cauto = F, cmin = marker_gene_exprs$range[1], cmax = marker_gene_exprs$range[2])) %>%
#                    layout(title = marker_gene_exprs$names,
#                           scene = list(xaxis = list(title = 'Real (Unit Circle of reCAT)'),
#                                        yaxis = list(title = 'Imaginary (Unit Circle of reCAT'),
#                                        zaxis = list(title = 'Differential Time')))#,
#                           #annotations = list(
#                           #  title = marker_gene_exprs$names,
#                           #  xref = 'paper',
#                           #  yref = 'paper',
#                           #  showarrow = FALSE
#                           #))
#       #scale_color_gradient(low = "white", high = "red"))
#       show(p)
#       return(p)
#     }
# 
#   }
# 
# }