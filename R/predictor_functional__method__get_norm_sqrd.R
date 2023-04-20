setMethod("get_norm_sqrd", "predictor_functional",
###########################################

function(input){
  n <- dim(input@coef)[1] #sample size
  norm_sqrd <- 0
  for (i in 1:n){
    coef_now <- matrix(input@coef[i,], nrow = 1)
    norm_sqrd <- norm_sqrd + sum(
      (t(coef_now) %*% coef_now) * input@J
      )
  }
return(norm_sqrd)
}

###########################################
)
