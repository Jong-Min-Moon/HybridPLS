n_eval <- 201
eval_point <- seq(0,1, length = n_eval)

n_sample <- 50
ar_slope_1 <- 0.9
ar_slope_2 <- 0.5
functional_data_eval_1 <- functional_data_eval_2 <-  matrix(NA, nrow = n_sample, ncol = n_eval)
for (i in 1:n_sample){
  functional_data_eval_1[i,] <- arima.sim(model = list(ar = ar_slope_1), n = n_eval)
  functional_data_eval_2[i,] <- arima.sim(model = list(ar = ar_slope_2), n = n_eval)
}


#%%%%%%%%%%%%% THINGS TO CHANGE %%%%%%%%

freq_array_low    = 5:45; # Hz
freq_array_high    = 10:400; # Hz
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# common parameters for both of low and high
gaussian_bandwidth = 5;
time        = eval_point;
srate      = n_eval;
wavelet_half_length = 19/n_eval
time_wavelet = seq(-wavelet_half_length, wavelet_half_length, length = 2* wavelet_half_length * n_eval)
n_wavelet   = length(time_wavelet);
n_data  = length(time);
n_convolution = n_data + n_wavelet - 1;
n_half_wavelet = (n_wavelet - 1)/2



# signal for low
filtered_analytic_signal = matrix(NA, nrow = n_sample, ncol = n_eval)

# loop for low
frequency_now = 10
for ( ii in 1:n_sample){
  signal_now = functional_data_eval_1[ii,]
  fft_signal_now = fft(signal_now)
  fft_signal_now_padded <- c(fft_signal_now, rep(0, n_convolution - n_data))

  wavelet_gaussian_sd_now = gaussian_bandwidth / (2 * pi * frequency_now)
  wavelet_gaussian_height_now = (wavelet_gaussian_sd_now * sqrt(pi))^(-1/2)
  wavelet_gaussian_part_now     = wavelet_gaussian_height_now * exp(-time_wavelet ^ 2 / (2 * (wavelet_gaussian_sd_now^2)));
  wavelet_complex_wave_part_now = exp(2 * 1i * pi * frequency_now * time_wavelet);
  wavelet_now                   = wavelet_gaussian_part_now * wavelet_complex_wave_part_now;
  wavelet_now_padded <- c(wavelet_now, rep(0, n_convolution - n_wavelet))
  convolution_result_now   = ifft(fft_signal_now_padded * fft(wavelet_now_padded));
  convolution_result_now = convolution_result_now[ (n_half_wavelet+1) : (n_convolution - n_half_wavelet)]
  filtered_analytic_signal[fi, ] = convolution_result_now;
}

head(filtered_analytic_signal_low)
plot(Re(fft(filtered_analytic_signal_low[1,])), pch=".")
plot(Re(fft(filtered_analytic_signal_low[40,])), pch=".")

plot(functional_data_eval_1[1,], pch=".")
plot(abs(filtered_analytic_signal_low[10,]), pch=".")
plot(abs(filtered_analytic_signal_low[15,]), pch=".")
