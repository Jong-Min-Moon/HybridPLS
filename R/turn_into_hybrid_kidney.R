turn_into_hybrid_kidney <- function(variables_kidney, n_basis){
  return(
    hybrid_kidney_from_data(
      argvals_base = variables_kidney$reno_base_x,
      argvals_post = variables_kidney$reno_post_x,
      reno_base = variables_kidney$reno_base_value,
      reno_post = variables_kidney$reno_post_value,
      scalar_predictors = variables_kidney$scalar_predictors,
      n_basis = n_basis
    )
  )

}
