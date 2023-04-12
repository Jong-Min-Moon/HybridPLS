.get_J_half <- function(J) {
  # gram matrix is positive definite,
  # so eigenvalues can be square-rooted, and eigenvectors are orthonormal.
  J_eigen <- eigen(J)
  S_half <- diag(sqrt(J_eigen$values))
  V <- J_eigen$vectors
  J_half <- V %*% S_half %*% t(V)
  return(J_half)
}
