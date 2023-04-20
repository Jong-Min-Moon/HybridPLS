setMethod("hybrid_from_coef", signature("predictor_hybrid"),
#######################################

function(format, xi_star){
  M <- dim(format@predictor_functional_list[[1]]@coef)[2] #number of basis
  format@predictor_functional_list[[1]]@coef <- t(as.matrix(xi_star[1:M]))
  format@predictor_functional_list[[2]]@coef <- t(as.matrix(xi_star[(M + 1):(2 * M)]))
  format@Z <- t(as.matrix(xi_star[(2 * M + 1):length(xi_star)]))
  return(format)
  }

#######################################
)
