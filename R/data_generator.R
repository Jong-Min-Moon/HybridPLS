generate_ar_wavelet_conv <- function(n_sample, n_eval,
                                     ar_slope,
                                     gaussian_bandwidth, fd_freq_list, scalar_freq_list
                                     ){
  eval_point <- seq(0,1, length = n_eval)

  #generate original signal using AR(1)
  signal.original <-  matrix(NA, nrow = n_sample, ncol = n_eval)
  for (i in 1:n_sample){
    signal.original[i,] <- arima.sim(model = list(ar = ar_slope), n = n_eval)
  }

  signal.filtered <- filter_signal(gaussian_bandwidth,
                                   c(fd_freq_list, scalar_freq_list),
                                   signal.original, eval_point, 3)
  scalar.predictor <- matrix(NA, nrow = n_sample, ncol = length(scalar_freq_list))
  for (ii in 1:length(scalar_freq_list)){
    signal <- signal.filtered[[ii+2]]
    scalar.predictor[,ii] <- apply(signal, 1, mean)
  }

  return(
    list(
      "x_functional" = list(
        "first" = list("value" = signal.filtered[[1]], "timestamp" = t(eval_point %*% t(rep(1, n_sample)))),
        "second" = list("value" = signal.filtered[[2]], "timestamp" = t(eval_point %*% t(rep(1, n_sample))))
      ),
      "x_scalar" = scalar.predictor
    )
  )
}






