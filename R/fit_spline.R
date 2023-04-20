fit_spline <- function(argvals, evals, n_basis){

  # create 1d b-spline basis functions
  PhiB <- bSpline(argvals, df=n_basis, degree = 3, intercept = TRUE) #Jt X M B spline Basis
  cat(paste("Using", n_basis, "number of basis"))

  # calculate basis coefficients for observed functional data
  C <- t(ginv(PhiB) %*% t(evals))

  # calculate gram matrices for
  # - basis functions (J)
  # - 2nd order derivative of basis functions (J'')
  J <- get_gram(argvals, PhiB)
  J_dotdot <- get_gram(argvals, deriv(PhiB, 2))

  return(
    list(
      C = C,
      J = J,
      J_dotdot = J_dotdot
    )
  )
}
