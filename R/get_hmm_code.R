# get HMM segmentation
#output <- c()
#for(i in 1:nrow(schlitz_logged_GO_subset)){

#  store_cor <- cor(order(schlitz_test14), schlitz_logged_GO_subset[i, ])
#  output <- c(output, store_cor)
#}

#source("HMM.R")
#source("baum_welch_scale.R")
#source("baum_welch_scale_log.R")


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

get_region <- function(list, len)
{
  if (list[1] < list[2])
  {
    return(c(list[1]:list[2]))
  }
  else
  {
    return(c(list[1]:len, 1:list[2]))
  }
}

normalDis <- function(val, cls_num)
{
  len = dim(val)[1]
  t = dim(val)[2]
  m = apply(val, 2, mean)
  #print(m)
  v = apply(val, 2, function(x) var(x))#*(len-1)/len)
  #v = apply(val, 2, function(x) (confidence_interval(x, 0.995)[2] - m))
  if((v == 0)){
    v <- 0.0001
  }
  result <- data.frame(m = m, v = v)
  return(result)
}

get_HMM_order_diff <- function(exprs_val, ordIndex, cls_num)
{
  #pred2 = bayes_score[ordIndex, ]
  #result2 = mean_score[ordIndex, ]
  
  gene_exprs <- exprs_val[ordIndex]

  
  if (cls_num == 3)
  {

    transPro = matrix(0, 3, 3)
    transPro[1,1] = 0.5
    transPro[1,2] = 0.5
    transPro[2,2] = 0.5
    transPro[2,3] = 0.5
    transPro[3,3] = 0.5
    transPro[3,1] = 0.5
    #transPro[3,3] = 1
    
    #mypi = c(1/3, 1/3, 1/3, 1/3)
    mypi = c(1/3, 1/3, 1/3)
    
    
    total_cells <- length(gene_exprs)
    sorted_val <- sort(gene_exprs)
    set_percentage <- 0.33 * length(ordIndex)
    
    split_n = round(total_cells / cls_num)
    
    low.id = c(1:round(set_percentage))
    mid_tmp <- round(total_cells / 2)
    med.id <- c((mid_tmp - round(split_n / 2)):(mid_tmp + round(split_n / 2)))
    high.id <- c((total_cells - round(set_percentage)): total_cells)
    #low.id = c(1:split_n)
    #med.id = c((split_n+1):(2*split_n))
    #high.id = c((max(med.id)+1):total_cells)
    
    h1 <- cutree(hclust(dist(as.numeric(gene_exprs))), cls_num)
    m1 <- mean(gene_exprs[which(h1 == 1)])
    m2 <- mean(gene_exprs[which(h1 == 2)])
    m3 <- mean(gene_exprs[which(h1 == 3)])
    
    m_vec <- c(m1, m2, m3)
    l_ind <- which.min(m_vec)
    h_ind <- which.max(m_vec)
    m_ind <- c(1:3)[-c(l_ind, h_ind)]#which(!c(l_ind, h_ind)%in%c(1:3))
    
    if((length(which(h1==l_ind))==1)){
      w_m <- which.min(gene_exprs[which(h1==m_ind)])[1]
      h1[which(h1==m_ind)[w_m]] <- m_ind
    }else if(length(which(h1==h_ind))==1){
      w_m <- which.max(gene_exprs[which(h1==m_ind)])[1]
      h1[which(h1==m_ind)[w_m]] <- m_ind
    }else if(length(which(h1==m_ind))==1){
      w_h <- which.min(gene_exprs[which(h1==h_ind)])[1]
      h1[which(h1==h_ind)[w_h]] <- h_ind
    }
    low_val = data.frame(x1 = gene_exprs[which(h1==l_ind)])
    med_val = data.frame(x1 = gene_exprs[which(h1 == m_ind)])
    high_val = data.frame(x1 = gene_exprs[which(h1 == h_ind)])
    
    #low_val = data.frame(x1 = gene_exprs[low.id])
    #med_val = data.frame(x1 = gene_exprs[med.id])
    #high_val = data.frame(x1 = gene_exprs[high.id])

    
    low_nd = normalDis(low_val)
    med_nd = normalDis(med_val)
    high_nd = normalDis(high_val)
    
    # val_vec <- c(low_nd$m, med_nd$m, high_nd$m)
    # val_mat <- rbind(low_nd, med_nd, high_nd)
    # reorder_val <- order(val_vec)
    # low_nd <- val_mat[reorder_val[1],]
    # med_nd <- val_mat[reorder_val[2],]
    # high_nd <- val_mat[reorder_val[3],]
    
    ndresult <- rbind(low_nd, med_nd, high_nd)
    
    #mu_state <- rbind(mu_state, c(low_nd[1], med_nd[1], high_nd[1]))

  }
  
  
  
  
  else if (cls_num == 2)
  {
    transPro = matrix(0, 2, 2)
    transPro[1,1] = 0.5
    transPro[1,2] = 0.5
    transPro[2,2] = 0.5
    transPro[2,1] = 0.5
    #transPro[2,2] = 1
    
    mypi = c(1/2, 1/2)
    

    total_cells <- length(gene_exprs)
    sorted_val <- sort(gene_exprs)
    set_percentage <- total_cells * 0.4
    
    #low.id = c(1:round(total_cells / cls_num))
    #high.id = c((max(low.id)+1):total_cells)
    low.id = c(1:round(set_percentage))
    high.id = c((total_cells - round(set_percentage)):total_cells)
    #print(low.id)
    #print(high.id)
    h1 <- cutree(hclust(dist(as.numeric(gene_exprs))), cls_num)
    
    m1 <- mean(gene_exprs[which(h1==1)])
    m2 <- mean(gene_exprs[which(h1 == 2)])
    m_vec <- c(m1, m2)
    l_ind <- which.min(m_vec)
    h_ind <- which.max(m_vec)
    if((length(which(h1==l_ind))==1)){
      w_h <- which.min(gene_exprs[which(h1==h_ind)])[1]
      h1[which(h1==h_ind)[w_h]] <- l_ind
    }else if(length(which(h1==h_ind))==1){
      w_l <- which.max(gene_exprs[which(h1==l_ind)])[1]
      h1[which(h1==l_ind)[w_l]] <- h_ind
    }
    low_val = data.frame(x1 = gene_exprs[which(h1==l_ind)])
    high_val = data.frame(x1 = gene_exprs[which(h1 == h_ind)])
    # k1 <- kmeans(as.numeric(gene_exprs), cls_num)
    # k1_clusters <- k1$cluster
    # g1 <- gene_exprs[which(k1_clusters == 1)]
    # g2 <- gene_exprs[which(k1_clusters == 2)]
    # if(mean(g1) <= mean(g2)){
    #   low_val = data.frame(x1 = g1)
    #   high_val = data.frame(x1 = g2)
    # }else{
    #   low_val = data.frame(x1 = g2)
    #   high_val = data.frame(x1 = g1)
    # }
    
    #low_val = data.frame(x1 = gene_exprs[low.id])
    #high_val = data.frame(x1 = gene_exprs[high.id])

    #print(low_val)
    #print(high_val)
    
    low_nd = normalDis(low_val)
    high_nd = normalDis(high_val)
    
    # val_vec <- c(low_nd$m, high_nd$m)
    # val_mat <- rbind(low_nd, high_nd)
    # reorder_val <- order(val_vec)
    # low_nd <- val_mat[reorder_val[1],]
    # high_nd <- val_mat[reorder_val[2],]

    
    ndresult <- rbind(low_nd, high_nd)
    print(ndresult)
    
    #mu_state <- rbind(mu_state, c(low_nd[1], high_nd[1]))
  }
  #print("yes")
  #print(ndresult)
  #print(transProb)
  #print(mypi)
  obs_val <- as.matrix(data.frame(gene_exprs))
  
  iter_max = 30#30
  myresult = myBW(A = transPro, nd_para = ndresult, mypi = mypi, ob_value = obs_val, iter_max = iter_max)
  # decoding viterbi
  rs1 = myviterbi(obval = obs_val, transProb = myresult$transProb, ndresult = myresult$nd, mypi = myresult$mypi, M = cls_num)
  #print(rs1)
  nowr = order(rs1$rmat[length(ordIndex),])[cls_num]
  rr = c(nowr)
  for (i in (length(ordIndex)-1):1)
  {
    rr = c(rr,rs1$pmat[i,nowr])
    nowr = rs1$pmat[i,nowr]
  }
  rr <- rev(rr)
  rr[1] = order(rs1$rmat[1,])[cls_num]
  #print(ndresult)
  #print(myresult$transProb)
  #print(myresult$mypi)
  #print(myresult$nd)
  return(bw_result = rr)
  
}
