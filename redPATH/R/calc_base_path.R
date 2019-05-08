# Calculate the base path for redPATH

calc_base_path <- function(tmp_normResultLst, tmpResultLst, base_path_range){
  store_rev <- c()
  store_rev_alt <- c()
  tmpCorLst <- c()
  tmpCorLst_alt <- c()
  
  ri <- tmp_normResultLst[1, ]

  #print(tmpResultLst)
  pairwise_iter <- 2
  for(j in 1:pairwise_iter){
    if(j < pairwise_iter){
      rj <- tmpResultLst[2, ]
    }else{
      rj <- rev(tmpResultLst[2,])
    }
    for(k in 1:pairwise_iter){
      if(k < pairwise_iter){
        rk <- tmpResultLst[3, ]
      }else{
        rk <- rev(tmpResultLst[3, ])
      }
      for(l in 1:pairwise_iter){
        if(l < pairwise_iter){
          rl <- tmpResultLst[4, ]
        }else{
          rl <- rev(tmpResultLst[4, ])
        }
        for(m in 1:pairwise_iter){
          if(m < pairwise_iter){
            rm <- tmpResultLst[5, ]
          }else{
            rm <- rev(tmpResultLst[5, ])
          }
          store_rev <- rbind(store_rev, c(j, k, l, m))
          #print(store_rev)
          tmpCor <- cor(ri, rj, method = "spearman")#, method = "spearman")
          tmpCor <- tmpCor + cor(ri, rk, method = "spearman")#, method = "spearman")
          tmpCor <- tmpCor + cor(ri, rl, method = "spearman")#, method = "spearman")
          tmpCor <- tmpCor + cor(ri, rm, method = "spearman")
          tmpCor <- tmpCor + cor(rj, rk, method = "spearman")#, method = "spearman")
          tmpCor <- tmpCor + cor(rj, rl, method = "spearman")#, method = "spearman")
          tmpCor <- tmpCor + cor(rj, rm, method = "spearman")
          tmpCor <- tmpCor + cor(rl, rk, method = "spearman")#, method = "spearman")
          tmpCor <- tmpCor + cor(rl, rm, method = "spearman")
          tmpCor <- tmpCor + cor(rk, rm, method = "spearman")
          
          tmpCorLst <- c(tmpCorLst, tmpCor)
        }
      }
      #print(cor(ri, rj))
    }
  }
  

  highNum <- which((tmpCorLst) == max((tmpCorLst)))

  rev_result <- store_rev[highNum, ]
  
  if(rev_result[1] == 1){
    
    tmpResultLst[2, ] <- tmp_normResultLst[2, ]
  }else{
    tmpResultLst[2, ] <- normalizeResult2(base_path_range[2], (base_path_range[2] + 1 - tmpResultLst[2, ]), tourlength = NULL, ordIndex = NULL) #rev(tmpResultLst[2, ])
  }
  if(rev_result[2] == 1){
    tmpResultLst[3, ] <- tmp_normResultLst[3, ]
  }else{
    tmpResultLst[3, ] <- normalizeResult2(base_path_range[3], (base_path_range[3] + 1 - tmpResultLst[3, ]), tourlength = NULL, ordIndex = NULL)#rev(tmpResultLst[3, ])
  }
  if(rev_result[3] == 1){
    tmpResultLst[4, ] <- tmp_normResultLst[4, ]
  }else{
    tmpResultLst[4, ] <- normalizeResult2(base_path_range[4], (base_path_range[4] + 1 - tmpResultLst[4, ]), tourlength = NULL, ordIndex = NULL)
  }
  if(rev_result[4] == 1){
    tmpResultLst[5, ] <- tmp_normResultLst[5, ]
  }else{
    tmpResultLst[5, ] <- normalizeResult2(base_path_range[5], (base_path_range[5] + 1 - tmpResultLst[5, ]), tourlength = NULL, ordIndex = NULL)
  }
  #print(range(tmpResultLst))
  return(tmpResultLst)
}