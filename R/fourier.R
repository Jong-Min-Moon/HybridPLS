filter_signal <- function(
    gaussian_bandwidth, frequency_list,
    signal.original, time,
    wavelet.time.length
    ){
  n.timepoint.original <- length(time)
  srate <- n.timepoint.original/(time[n.timepoint.original] - time[1])
  n.timepoint.wavelet = srate * wavelet.time.length
  signal.concat <- c(t(signal.original))
  n.timepoint.concat <- length(signal.concat)
  n.timepoint.conv = n.timepoint.concat + n.timepoint.wavelet - 1;


  time.wavelet <- seq(-wavelet.time.length/2, wavelet.time.length/2, length = n.timepoint.wavelet)


  wavelet.mat <- complex_morlet_wavelet(gaussian_bandwidth, frequency_list, time.wavelet)
  wavelet.mat.padded <- cbind(wavelet.mat, matrix(0, nrow = length(frequency_list), ncol = n.timepoint.conv - n.timepoint.wavelet))
  wavelet.fft <- fft.rowwize(wavelet.mat.padded)

  signal.concat.fft <- fft(
    c(signal.concat, rep(0, (n.timepoint.conv - n.timepoint.concat)))
  )
  signal.concat.fft <- t(
    signal.concat.fft %*% matrix( rep(1, length(frequency_list)), nrow = 1)
    )

  conv.result <- ifft.rowwize(wavelet.fft * signal.concat.fft)
  n.timepoint.wavelet.half = (n.timepoint.wavelet - 1)/2
  conv.result <- conv.result[ ,(n.timepoint.wavelet.half+1) : (n.timepoint.conv - n.timepoint.wavelet.half)]
  conv.result <- abs(conv.result)

  result <- list()
  for (ii in 1:length(frequency_list)){
    freq <- frequency_list[ii]
    result[[ii]] <- matrix(conv.result[ii,], ncol = n.timepoint.original, byrow = T)
  }
  return(result)
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


fft.rowwize <- function(X) t(mvfft(t(X)))
ifft.rowwize <- function(F) t(mvfft(t(F), inverse=TRUE)) / ncol(F)
