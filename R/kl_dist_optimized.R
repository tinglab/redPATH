# KL dist optimization


kl_trial <- function(object, eps=10^-4, overlap=TRUE,...)
{
  if(!is.numeric(object))
    stop("object must be a numeric matrix\n")
  
  z <- matrix(NA, nrow=ncol(object), ncol=ncol(object))
  colnames(z) <- rownames(z) <- colnames(object)
  
  w <- object < eps
  if (any(w)) object[w] <- eps
  #object <- sweep(object, 2, colSums(object) , "/")
  object <- object / colSums(object)
  
  for(k in seq_len(ncol(object)-1)){
    for(l in 2:ncol(object)){
      #ok <- (object[, k] > eps) & (object[, l] > eps)
      #print(ok)
      #if (!overlap | any(ok)) {
        z[k,l] <- sum(object[,k] *
                        (log(object[,k]) - log(object[,l])))
        z[l,k] <- sum(object[,l] *
                        (log(object[,l]) - log(object[,k])))
      #}
    }
  }
  diag(z)<-0
  z
}
#' Fast implementation of KL distance
#' 
#' This calculates the KL-distance between cells, generating the asymmetric distance matrix
#'@param object Input matrix
#'@export
#'@keywords kl_cpp
#'@section redPATH pseudotime:
#'@examples
#'kl_mat <- kl_cpp(data)

kl_cpp <- function(object, eps=10^-4, overlap=TRUE,...)
{
  if(!is.numeric(object))
    stop("object must be a numeric matrix\n")

  w <- object < eps
  if (any(w)) object[w] <- eps
  object_c <- sweep(object, 2, colSums(object) , "/")

  z = KL(object_c, log(object_c))
  colnames(z) <- rownames(z) <- colnames(object)

  return(z)
}
kl_old <- function(object, eps=10^-4, overlap=TRUE,...)
{
  if(!is.numeric(object))
    stop("object must be a numeric matrix\n")
  
  z <- matrix(NA, nrow=ncol(object), ncol=ncol(object))
  colnames(z) <- rownames(z) <- colnames(object)
  
  w <- object < eps
  if (any(w)) object[w] <- eps
  object <- sweep(object, 2, colSums(object) , "/")
  #assign("obj_o", object, envir = .GlobalEnv)
  #print(object)
  for(k in seq_len(ncol(object)-1)){
    for(l in 2:ncol(object)){
      #print(l)
      #print(object[,l])
      #print(k)
      #print(object[,k])
      ok <- (object[, k] > eps) & (object[, l] > eps)
      if (!overlap | any(ok)) {
        z[k,l] <- sum(object[,k] *
                        (log(object[,k]) - log(object[,l])))
        z[l,k] <- sum(object[,l] *
                        (log(object[,l]) - log(object[,k])))
      }
    }
  }
  diag(z)<-0
  z
}