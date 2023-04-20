setMethod("hybrid_inner_prod", signature("predictor_hybrid", "predictor_hybrid"),
#####################################################################

function(xi_1, xi_2) {
if (is_one_observation(xi_2) &
    is_same_basis(xi_1, xi_2)
){
  ip_functional <- 0
  for (i in 1:(xi_1@n_predictor_functional)){
    coeff_1 <- xi_1@predictor_functional_list[[i]]@coef
    coeff_2 <- xi_2@predictor_functional_list[[i]]@coef
    J <- xi_1@predictor_functional_list[[i]]@J
    ip_functional <- ip_functional + as.numeric(coeff_1  %*% J  %*% t(coeff_2))
    }
  ip_scalar <- as.numeric(xi_1@Z %*% t(xi_2@Z))
  return(ip_functional + ip_scalar)
}else{
  stop("No observation or different bases")
}

}
###################################################################################################
)
