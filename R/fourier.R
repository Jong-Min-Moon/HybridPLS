
filter_signal <- function(signal.original, eval_point, frequency_now, gaussian_bandwidth){
n_sample = nrow(signal.original)
n_eval = ncol(signal.original)
time        = eval_point;
srate      = n_eval;
wavelet_half_length = 19/n_eval
time = seq(-wavelet_half_length, wavelet_half_length, length = 2* wavelet_half_length * n_eval)
n_wavelet   = length(time);
n_data  = length(time);
n_convolution = n_data + n_wavelet - 1;
n_half_wavelet = (n_wavelet - 1)/2
filtered_analytic_signal = matrix(NA, nrow = n_sample, ncol = n_eval)

for ( ii in 1:n_sample){
  signal_now = signal.original[ii,]
  fft_signal_now = fft(signal_now)
  fft_signal_now_padded <- c(fft_signal_now, rep(0, n_convolution - n_data))


  wavelet_now_padded <- c(wavelet_now, rep(0, n_convolution - n_wavelet))
  convolution_result_now   = ifft(fft_signal_now_padded * fft(wavelet_now_padded));
  convolution_result_now = convolution_result_now[ (n_half_wavelet+1) : (n_convolution - n_half_wavelet)]
  filtered_analytic_signal[ii, ] = Re(convolution_result_now);
}
return(filtered_analytic_signal)
}



# gaussian_bandwidth: scalar
# frequency: vector of length = f
# time: vector of length = t
complex_morlet_wavelet <- function(gaussian_bandwidth, frequency_list, time){
  n_timepoint <- length(time) # t
  n_frequency <- length(frequency_list) # f
  time_mat <- t(matrix(time) %*% matrix(rep(1, n_frequency), nrow = 1)) # f x t matrix

  gaussian_sd_colmat <- gaussian_bandwidth / (2 * pi * matrix(frequency_list) ) # f x 1 matrix
  guassian_sd_mat <- gaussian_sd_colmat %*% matrix(rep(1,n_timepoint), nrow = 1) # f x t matrix
  gaussian_height_mat = (guassian_sd_mat * sqrt(pi))^(-1/2) # f x t matrix

                   # f x t matrix           # f x t matrix         # f x t matrix
  gaussian_part = gaussian_height_mat * exp(-time_mat ^ 2 / (2 * (guassian_sd_mat^2)));
  complex_wave_part = exp(2 * 1i * pi * matrix(frequency_list) %*% t(matrix(time)));

    return(gaussian_part * complex_wave_part)
}


