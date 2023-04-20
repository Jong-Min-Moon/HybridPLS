response_transform_min_max_logit <- function(y, min, max){
  y <- (y - min)/(max-min)
  y[y==1] <- 0.99
  y[y==0] <- 0.01
  y <- log(y/(1-y))
  return(y)
}
