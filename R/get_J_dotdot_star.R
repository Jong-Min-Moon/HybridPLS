# get_J_dotdot_star
#
#


setMethod("get_J_dotdot_star", signature("predictor_hybrid"),
#######################################

function(W){
  p <- dim(W@Z)[2]
  J_dotdot_star <- W@predictor_functional_list[[1]]@J_dotdot
  for (i in 2:(W@n_predictor_functional)){
    J_dotdot_star <- bdiag(J_dotdot_star, W@predictor_functional_list[[i]]@J_dotdot)
    }
  J_dotdot_star <- bdiag(J_dotdot_star, matrix(0, nrow = p, ncol = p))
  return(J_dotdot_star)
  }

#######################################
)
