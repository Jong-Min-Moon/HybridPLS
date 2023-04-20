setMethod("index_sample", signature("predictor_hybrid"),
################################################################
function(W, i){
  n <- dim(W@Z)[1]
  if(i > n){stop("index out of range")}
  W@Z <- matrix(W@Z[i,], nrow = 1)
  for (j in 1 : n_predictor_functional){
    W@predictor_functional_list[[i]]@coef <- matrix(W@predictor_functional_list[[i]]@coef[j,], nrow = 1)
    }
  return(W)
  }

################################################################
)
