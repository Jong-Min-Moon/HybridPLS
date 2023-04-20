setMethod("hybrid_norm",
          signature("predictor_hybrid"),
          function(xi) {
            return( sqrt(hybrid_inner_prod(xi, xi)) )}
          )
