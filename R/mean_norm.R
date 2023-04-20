setMethod("mean_norm", signature("predictor_hybrid"),
######################################################
function(W){
            norm_W = 0
            n = dim(W@Z)[1]
            print(n)
            for (i in 1:n){
              w_now <- index_sample(W, i)
              norm_W <- norm_W + hybrid_norm(w_now)
            }
            return(norm_W/n)
          }

######################################################
)
