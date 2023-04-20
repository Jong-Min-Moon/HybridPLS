# get_J_star
#
#


setMethod("get_J_star", signature("predictor_hybrid"),
#######################################

function(W){
  p <- dim(W@Z)[2]
  J_star <- W@predictor_functional_list[[1]]@J
  for (i in 2:(W@n_predictor_functional)){
    J_star <- bdiag(J_star, W@predictor_functional_list[[i]]@J)
  }
  J_star <- bdiag(J_star, diag(p))
  return(J_star)
}

#######################################
)
