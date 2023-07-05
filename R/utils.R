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
  cat(paste("use", n_basis, "basis functions\n"))
  # create 1d b-spline basis functions
  #PhiB <- splines2::bSpline(argvals, df=n_basis, degree = 3, intercept = TRUE) #Jt(n.eval) X M B spline Basis #old way
  my_basis <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis)
  PhiB <- predict(my_basis, argvals)
  PhiB <- PhiB[nrow(PhiB):1, ncol(PhiB):1]

  # calculate basis coefficients for observed functional data # idea from Dr. Jang
  #C <- t( MASS::ginv(PhiB) %*% t(evals) ) #old way
  fda_obj <- Data2fd(argvals, t(evals), my_basis)
  C <- t(fda_obj$coefs)
  # calculate gram matrices for
  # - basis functions (J)
  # - 2nd order derivative of basis functions (J'')
  J <- get_gram_2d(argvals, PhiB)

  #J_dotdot <- get_gram_2d(argvals, deriv(PhiB, 2))
  deriv_PhiB <- predict(my_basis, argvals, deriv=2)
  deriv_PhiB <- deriv_PhiB[nrow(deriv_PhiB):1, ncol(deriv_PhiB):1]
  J_dotdot <- get_gram_2d(argvals, deriv_PhiB)
  return(
    list(
      C = C,
      J = J,
      J_dotdot = J_dotdot,
      fd = fda_obj
    )
  )
}

dataset_for_regression_refund <- function(response, hybrid_predictor){
  dataset_for_regression <- data.frame(hybrid_predictor@Z)
  dataset_for_regression$y <- response
  dataset_for_regression$precurve <-hybrid_predictor@predictor_functional_list[[1]]@original_X
  dataset_for_regression$postcurve <-hybrid_predictor@predictor_functional_list[[2]]@original_X
  return(dataset_for_regression)
}

dataset_for_prediction_refund <- function(hybrid_predictor){
  dataset_for_regression <- data.frame(hybrid_predictor@Z)
  dataset_for_regression$precurve <-hybrid_predictor@predictor_functional_list[[1]]@original_X
  dataset_for_regression$postcurve <-hybrid_predictor@predictor_functional_list[[2]]@original_X
  return(dataset_for_regression)
}
