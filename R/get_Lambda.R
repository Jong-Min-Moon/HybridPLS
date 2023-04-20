# get_Lambda
#
#
setMethod("get_Lambda", signature("predictor_hybrid", "numeric"),
#################################################################

function(W, kappa){
  p <- dim(W@Z)[2]
  n_basis_1 <- dim(W@predictor_functional_list[[1]]@J)[2]
  Lambda <- diag(n_basis_1)
  for (i in 2:(W@n_predictor_functional)){
    n_basis_now <- dim(W@predictor_functional_list[[i]]@J)[2]
    Lambda <- bdiag(Lambda, diag(n_basis_now))
    }
  Lambda <- bdiag(
    Lambda,
    matrix(0, nrow = p, ncol = p)
  )
  }

#######################################

)
