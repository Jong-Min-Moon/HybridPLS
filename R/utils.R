get_traning_idx <- function(n_sample, training_ratio){
  return( sample(n_sample, floor(training_ratio * n_sample)) )
}

get_gram_2d <- function(argval, basis_evals){
  n_basis <- dim(basis_evals)[2]
  J <- matrix(nrow= n_basis, ncol = n_basis)
  for (i in 1:n_basis) {
    for (j in 1:n_basis){
      J[i,j] <- pracma::trapz(argval, basis_evals[,i]*basis_evals[,j])
    }}
  return(J)
}

get_gram_half <- function(J) {
  # gram matrix is positive definite,
  # so eigenvalues can be square-rooted, and eigenvectors are orthonormal.
  J_eigen <- eigen(J)
  S_half <- diag(sqrt(J_eigen$values))
  V <- J_eigen$vectors
  J_half <- V %*% S_half %*% t(V)
  return(J_half)
}

fit_spine_2d <- function(argvals, evals, n_basis){

  # create 1d b-spline basis functions
  PhiB <- splines2::bSpline(argvals, df=n_basis, degree = 3, intercept = TRUE) #Jt X M B spline Basis
  cat(paste("use", n_basis, "basis functions\n"))

  # calculate basis coefficients for observed functional data
  # idea from Dr. Jang
  C <- t(MASS::ginv(PhiB) %*% t(evals))

  # calculate gram matrices for
  # - basis functions (J)
  # - 2nd order derivative of basis functions (J'')
  J <- get_gram_2d(argvals, PhiB)
  J_dotdot <- get_gram_2d(argvals, deriv(PhiB, 2))

  return(
    list(
      C = C,
      J = J,
      J_dotdot = J_dotdot
    )
  )
}
