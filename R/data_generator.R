generate_ar_wavelet_conv <- function(n_sample, n_eval,
                                     ar_slope, gaussian_bandwidth){
  #generate original signal using AR(1)
  signal.original <-  matrix(NA, nrow = n_sample, ncol = n_eval)
  for (i in 1:n_sample){
    signal.original[i,] <- arima.sim(model = list(ar = ar_slope), n = n_eval)

  }
}






time        = eval_point;
srate      = n_eval;
wavelet_half_length = 19/n_eval
time_wavelet = seq(-wavelet_half_length, wavelet_half_length, length = 2* wavelet_half_length * n_eval)
n_wavelet   = length(time_wavelet);
n_data  = length(time);
n_convolution = n_data + n_wavelet - 1;
n_half_wavelet = (n_wavelet - 1)/2

frequency_1 = 20
frequency_2 = 20


functional_data_eval_1 = filter_signal(signal.original, eval_point, frequency_1, gaussian_bandwidth)
functional_data_eval_2 = filter_signal(signal.original, eval_point, frequency_2, gaussian_bandwidth)
