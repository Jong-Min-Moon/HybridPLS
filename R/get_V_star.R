# get_V_star
#
# computes V*,
# which is the core part of the matrix representation of the covariance,
# which is the objective function of the optimization problem
#
# input:
# W: hybrid predictor object
# y: response vector
#
# output:
# matrix V* (from Proposition 1)
# which is a (MK+p) x (MK+p) square matrix,
# where
# M = number of basis
# K = dimension of the functional vector. In our case, K = 3(R, G, B)
# p = dimension of covariate.
#


setMethod("get_V_star", signature("predictor_hybrid", "numeric"),
#######################################

function(W, y){
  n <- W@n_sample
  K <- W@n_predictor_functional
  # building blocks
  ## 1. C multiplied by J. Scalable to arbitrary K
  predictor_first <- W@predictor_functional_list[[1]]
  C_J <- predictor_first@coef %*% predictor_first@J
  for (i in 2:K){
    predictor_now <- W@predictor_functional_list[[i]]
    C_J <- cbind(
      C_J,
      predictor_now@coef %*% predictor_now@J
    )
  }

  ## 2. easy ones
  J_Ct_y <- t(C_J) %*% y
  Zt_y <- t(W@Z) %*% y

  # V matrix
  upper_left <- J_Ct_y %*% t(J_Ct_y)
  upper_right <- J_Ct_y %*% t(Zt_y)
  lower_left <- t(upper_right)
  lower_right <- Zt_y %*% t(Zt_y)

  V_star <- cbind(
    rbind(upper_left, lower_left),
    rbind(upper_right, lower_right)
    )/(n^2)

  return(V_star)
  }

#######################################
)
