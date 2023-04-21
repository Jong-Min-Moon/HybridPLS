setGeneric(
  "LSE_hybrid",
  function(W, rho, tau) standardGeneric("LSE_hybrid")
)

setMethod("LSE_hybrid", "predictor_hybrid",
###########################################

function(W, rho, tau) {
  p <- dim(W@Z)[2]

  rho_t_C_star_J_star <- (W@predictor_functional_list[[1]]@coef) %*% (W@predictor_functional_list[[1]]@J)
  for (i in (2 : W@n_predictor_functional)){
    rho_t_C_star_J_star <- cbind(
      rho_t_C_star_J_star,
      (W@predictor_functional_list[[i]]@coef) %*% (W@predictor_functional_list[[i]]@J)
      )
  }
  rho_t_C_star_J_star <- t(rho) %*% cbind(rho_t_C_star_J_star, W@Z)


  # matrices in the statement of Proposition 4
  J_star <-get_J_star(W)
  J_dotdot_star <- get_J_dotdot_star(W)
  convex_J <- (norm(rho, type="2"))^2 * J_star + tau * J_dotdot_star
  d_hat_star <- t(
    solve( t(convex_J), t(rho_t_C_star_J_star) )
    )
  delta_hat <- hybrid_from_coef(W, d_hat_star)
  return(delta_hat)
}

#########################################
)
