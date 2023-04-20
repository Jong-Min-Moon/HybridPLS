setMethod("scalar_mul", signature("predictor_hybrid", "numeric"),
################################################################
function(W, scalar) {
  for (i in 1:W@n_predictor_functional){
    W@predictor_functional_list[[i]]@coef <- scalar * W@predictor_functional_list[[i]]@coef
  }
  W@Z <- scalar*(W@Z)
  return(W)
  }
################################################################
)
