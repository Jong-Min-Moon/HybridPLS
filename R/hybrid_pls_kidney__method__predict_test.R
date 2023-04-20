setMethod("predict_test", "hybrid_pls_kidney",
          function(pls_object, W_test) {
            #centering
            return(
              y_pred = hybrid_inner_prod(W_test, pls_object@eta)
            )
          })
