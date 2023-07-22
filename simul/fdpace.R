set.seed(1)
n <- 50
pts <- seq(0, 1, by=0.05)
# The first n samples are for training and the rest testing
sampWiener <- Wiener(2 * n, pts)
sparsity <- 2:5
train <- Sparsify(sampWiener[seq_len(n), , drop=FALSE], pts, sparsity)
test <- Sparsify(sampWiener[seq(n + 1, 2 * n), , drop=FALSE], pts, sparsity)
res <- FPCA(train$Ly, train$Lt)
pred <- predict(res, test$Ly, test$Lt, K=3)
plot(pred$predGrid, pred$predCurves[1,])
