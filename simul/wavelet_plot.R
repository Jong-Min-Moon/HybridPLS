
n.eval <- 402
srate <- n.eval
time.length = 3
n.timepoint.wavelet = srate * time.length
time.wavelet <- seq(-time.length/2, time.length/2, length = n.timepoint.wavelet) # to ensure same srate
a<-complex_morlet_wavelet(3, scalar_freq, time.wavelet)

plot(Re(a[1,]), pch = ".")

plot(abs(fft(a[1,])), pch= "", xlim = c(0,100))
lines(abs(fft(a[1,])), pch= ".", col= "red")
lines(abs(fft(a[2,])), pch= ".",)
lines(abs(fft(a[3,])), pch= ".")
lines(abs(fft(a[4,])), pch= ".", col = "blue")
lines(abs(fft(a[5,])), pch= ".")
lines(abs(fft(a[6,])), pch= ".")
lines(abs(fft(a[7,])), pch= ".", col = "blue")

