setMethod("fitted_value",
          signature(
            "predictor_hybrid",
            "predictor_hybrid",
            "vector"
          ),
####################################################
function(data, coeff, pls_score) {
            #               n * 1           1 * M
            for (i in 1 : (data@n_predictor_functional) ){
              data@predictor_functional_list[[i]]@coef <- pls_score %*% coeff@predictor_functional_list[[i]]@coef
            }
            data@Z <-  pls_score %*% coeff@Z
            return(data)
          }

####################################################
)
