setMethod("get_norm_sqrd", "predictor_functional",
###########################################

function(input){
  n <- dim(input@coef)[1] #sample size
  norm_sqrd <- 0
  for (i in 1:n){
    coef_now <- input@coef[i,]
    coef_sqrd <- coef_now %*% coef_now
    norm_sqrd <- norm_sqrd + (
      (t(coef_now) %*% coef_now) * input@J
      )
  }
return(norm_sqrd)
}

###########################################
)
