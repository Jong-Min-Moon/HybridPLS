response_transform_min_max_logit <- function(y, min, max){
  y <- (y - min)/(max-min)
  y[y==1] <- 0.99
  y[y==0] <- 0.01
  y <- log(y/(1-y))
  return(y)
}

reponse_inverse_transform_min_max_logit <- function(y_pred_pls_old, y_train_logit_mean, y_train_max, y_train_min){
  y_pred_pls <- y_pred_pls_old + y_train_logit_mean
  y_pred_pls <- exp(y_pred_pls)/( 1 + exp(y_pred_pls))
  y_pred_pls <- y_pred_pls * (y_train_max- y_train_min)
  y_pred_pls <- y_pred_pls + y_train_min
}

