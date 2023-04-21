setClass(
  "hybrid_pls_kidney",
  representation(
    eta = "predictor_hybrid",
    xi = "list",
    nu = "list",
    rho = "list",
    delta = "list",
    E = "list",
    V_star = "list",
    eigen_val = "list",
    first_eigen_val = "vector",
    J_Lambda_Jpp = "dgCMatrix",
    L_mat = "dtCMatrix",
    resid_y = "list",
    mse_W = 'vector',
    mse_y = 'vector',
    fitted_value_W = "list",
    fitted_value_y = "list",
    W_now = "list"
  ))

setMethod("predict_test", "hybrid_pls_kidney",
          function(pls_object, W_test) {
            #centering
            return(
              y_pred = hybrid_inner_prod(W_test, pls_object@eta)
            )
          })
