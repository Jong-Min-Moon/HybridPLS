get_gram <- function(argval, basis_evals){
  n_basis <- dim(basis_evals)[2]
  J <- matrix(nrow= n_basis, ncol = n_basis)
  for (i in 1:n_basis) {
    for (j in 1:n_basis){
      J[i,j] <- trapz(argval, basis_evals[,i]*basis_evals[,j])
    }}
  return(J)
}
