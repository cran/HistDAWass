## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 4.5, dpi = 72, out.width = "100%")

## ----setup--------------------------------------------------------------------
library(HistDAWass)

## -----------------------------------------------------------------------------
BLOOD

## -----------------------------------------------------------------------------
d <- BLOOD@M[[1, 1]] # cholesterol of the first group
d
meanH(d)    # mean
stdH(d)     # standard deviation
skewH(d)    # skewness
kurtH(d)    # kurtosis (3 for a Gaussian)
compQ(d, 0.5)   # the median
compQ(d, 0.25); compQ(d, 0.75)  # first and third quartiles
compP(d, 150)   # the cdf at 150

## -----------------------------------------------------------------------------
str(crwtransform(d))

## -----------------------------------------------------------------------------
d1 <- BLOOD@M[[1, 1]]
d2 <- BLOOD@M[[2, 1]]
d12 <- 0.5 * d1 + 0.5 * d2   # the Wasserstein mean of the two histograms
c(meanH(d1), meanH(d2), meanH(d12))

## -----------------------------------------------------------------------------
WassSqDistH(d1, d2)
WassSqDistH(d1, d2, details = TRUE)
rQQ(d1, d2)

## -----------------------------------------------------------------------------
get.MatH.stats(BLOOD, stat = "mean")$mat
get.MatH.stats(BLOOD, stat = "std")$mat[1:5, ]
get.MatH.stats(BLOOD, stat = "quantile", prob = 0.9)$mat[1:5, ]

## ----fig.height = 3.5---------------------------------------------------------
m_chol <- WH.vec.mean(BLOOD[, "Cholesterol"])
m_chol
plot(m_chol)

## -----------------------------------------------------------------------------
w <- c(rep(2, 7), rep(1, 7))
meanH(WH.vec.mean(BLOOD[, 1], w = w))

## -----------------------------------------------------------------------------
WH.var.covar(BLOOD)
WH.correlation(BLOOD)
WH.SSQ(BLOOD)

## -----------------------------------------------------------------------------
WH.correlation2(BLOOD[, 1], BLOOD[, 2:3])

## -----------------------------------------------------------------------------
D <- WH_MAT_DIST(BLOOD, standardize = TRUE)
round(D[1:5, 1:5], 3)

## -----------------------------------------------------------------------------
u <- upper.tri(D)
for (q in c(20, 100)) {
  Dq <- WH_MAT_DIST(BLOOD, standardize = TRUE, simplify = TRUE, qua = q)
  err <- abs(Dq[u] - D[u]) / D[u]
  cat("qua =", q, " median:", signif(median(err), 2),
      " maximum:", signif(max(err), 2), "\n")
}

