# get_Lambda
#
# input:
# W
# lambda: K-dimensional vector
#
# transition from scalar kappa to vector lambda (created October 27th 2022)
#

setMethod("get_Lambda", signature("predictor_hybrid", "numeric"),
#################################################################

function(W, lambda){
  p <- dim(W@Z)[2]
  n_basis_1 <- dim(W@predictor_functional_list[[1]]@J)[2]
  Lambda <- lambda[1] * diag(n_basis_1)
  for (i in 2:(W@n_predictor_functional)){
    n_basis_now <- dim(W@predictor_functional_list[[i]]@J)[2]
    Lambda <- bdiag(Lambda, lambda[i] * diag(n_basis_now))
    }
  Lambda <- bdiag(
    Lambda,
    matrix(0, nrow = p, ncol = p)
  )
  }

#######################################

)
