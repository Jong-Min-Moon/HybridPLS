# pls_pen
# computes hybrid pls component (loading)
# by maximizing the covariance with regularization
# this function is called at each l-th iteration of PLS fitting

#input:
# W: hyprid predictor object
# y: response vector
# L: regularization matrix

setMethod("pls_pen", "predictor_hybrid",
#######################################

function(W, y, L){
  V_star <- get_V_star(W, y)
  #Last paragraph of section 3.2. To solve (6) in practice...
  A <- t(solve(L, t(V_star))) # A = V* t(inv(L))
  E <- solve(L, A) # E = inv(L) A = inv(L) V* t(inv(L))
  eigen_result <- eigen(E)
  cat("eigenvalues: ", eigen_result$values)
  e <- eigen_result$vectors[, 1]
  #e <- e / Matrix::norm(e, "2") #normalize; omit this by request of Dr. Jang

  xi_star <- solve(t(L), e) # t(L) xi* = e
  xi_hat <- hybrid_from_coef(format = W, xi_star = xi_star)


  return(
    list(
      xi = xi_hat,
      E = E,
      V_star = V_star,
      eigen_val = eigen_result$values
      )
    )

}

#######################################
)
