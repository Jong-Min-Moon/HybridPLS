library(MFPCA)
argvalsList <- list( seq(0, 1, 0.01), seq(0, 1, 0.01))
simMultiSplit <- simMultiFunData(
  N = 7,
  argvals = argvalsList,
  eFunType = "Fourier",
  eValType = "linear",
  M = 10,
  type = "split"
)
simMultiSplit

uniExpansions <- list(
  list(type = "uFPCA", npc = 4),
  list(type = "uFPCA", npc = 5)
  )
MFPCA_result <- MFPCA(simMultiSplit$simData, M = 1, uniExpansions = uniExpansions)

MFPCA_result$scores
MFPCA_result$vectors

MFPCA_result$normFactors

(MFPCA_result$functions) * simMultiSplit
MFPCA_result$functions[[1]][1]
simMultiSplit[[1]]
i=3
j=1
integrate(
  (simMultiSplit$simData[[1]][i] - MFPCA_result$meanFunction[[1]]) * (MFPCA_result$functions[[1]][j])
)+integrate(
    (simMultiSplit$simData[[2]][i]- MFPCA_result$meanFunction[[2]]) * (MFPCA_result$functions[[2]][j])
  )

MFPCA_result$scores
MFPCA_result$meanFunction[[1]]
