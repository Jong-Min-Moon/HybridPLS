ar_slope_1 <- 0.9
ar_slope_2 <- 0.5
indep_1 <- indep_2 <-  matrix(NA, nrow = n_sample, ncol = n_eval)
for (i in 1:n_sample){
  indep_1[i,] <- arima.sim(model = list(ar = ar_slope_1), n = n_eval)
  indep_2[i,] <- arima.sim(model = list(ar = ar_slope_2), n = n_eval)
}

n_sample <- 50
ar_slope <- 0.9
gaussian_bandwidth <- 3

freq_1 <- 6
freq_2 <- 12
scalar_freq <- seq(freq_1,freq_2,length = n_predictor_scalar+2)
predictor.value <- generate_ar_wavelet_conv(
  n_sample, n_eval, ar_slope, gaussian_bandwidth,
  c(freq_1,freq_2), scalar_freq[2:(n_predictor_scalar+1)]
)
dep_1<-predictor.value$x_functional$first$value
dep_2<-predictor.value$x_functional$second$value

cor_indep <- cor(indep_1, indep_2)
cor_dep <- cor(dep_1, dep_2)

hist(diag(cor_indep), xlim = c(-1,1))
hist(diag(cor_dep), xlim = c(-1,1))


ii <- (seq(1,10,length=1000))
plot(ii, log(ii))
plot(ii, log(log(ii)))
plot(ii, ii^{1/10})

