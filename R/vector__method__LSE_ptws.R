setMethod("LSE_ptws", "vector",
################################

function(input, rho) {
  rho_norm_sq <- (Matrix::norm(rho, type = "2"))^2
  nu <- as.numeric(t(rho) %*% input) / rho_norm_sq
  return(nu)
}

################################
)
