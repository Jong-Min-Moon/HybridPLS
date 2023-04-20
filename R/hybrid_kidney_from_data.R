################# Read functions for kidney data ###################

hybrid_kidney_from_data <- function(
    argvals_base,
    argvals_post,
    reno_base,
    reno_post,
    scalar_predictors,
    n_basis
){
  split_fit_base <- .fit_spline(argvals = argvals_base, evals = reno_base, n_basis = n_basis)
  split_fit_post <- .fit_spline(argvals = argvals_post, evals = reno_post, n_basis = n_basis)

  hybrid_predictor_dataset<-
    new("hybrid_predictors_kidney",
        basis_coeff_1 = split_fit_base$C,
        basis_coeff_2 = split_fit_post$C,
        Z = scalar_predictors,

        J_1 = split_fit_base$J,
        J_half_1 = get_J_half(split_fit_base$J),
        J_dotdot_1 = split_fit_base$J_dotdot,

        J_2 = split_fit_post$J,
        J_half_2 = get_J_half(split_fit_post$J),
        J_dotdot_2 = split_fit_post$J_dotdot
        #,raw = multiFunData(reno_base_fund, reno_post_fund)
    )
  return(hybrid_predictor_dataset)
}
