## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 4.5, dpi = 72, out.width = "100%")

## ----setup--------------------------------------------------------------------
library(HistDAWass)

## -----------------------------------------------------------------------------
get.MatH.main.info(OzoneFull)$varnames

## -----------------------------------------------------------------------------
set.seed(1)
km <- WH_kmeans(OzoneFull, k = 3, rep = 5,
                simplify = TRUE, qua = 20, standardize = TRUE)
names(km$solution)
km$solution$cardinality
km$quality   # proportion of the total sum of squares explained by the partition
head(km$solution$IDX)

## ----fig.height = 3.5---------------------------------------------------------
plot(km$solution$centers)

## ----results = "hide"---------------------------------------------------------
set.seed(1)
qual <- sapply(2:7, function(k)
  WH_kmeans(OzoneFull, k = k, rep = 5, simplify = TRUE, qua = 20,
            standardize = TRUE)$quality)

## ----fig.height = 3.5---------------------------------------------------------
round(qual, 3)
plot(2:7, qual, type = "b", xlab = "k", ylab = "quality")

## ----results = "hide"---------------------------------------------------------
set.seed(1)
akm <- WH_adaptive.kmeans(OzoneFull, k = 3, schema = 2, rep = 5,
                          simplify = TRUE, qua = 20, standardize = TRUE)

## -----------------------------------------------------------------------------
round(akm$weights[, 1], 3)   # the same for all the clusters with schema = 2
akm$cardinality
akm$quality

## ----results = "hide"---------------------------------------------------------
set.seed(1)
fcm <- WH_fcmeans(OzoneFull, k = 3, m = 1.6, rep = 5,
                  simplify = TRUE, qua = 20, standardize = TRUE)

## -----------------------------------------------------------------------------
round(head(fcm$solution$membership), 3)
fcm$solution$cardinality   # crisp assignment to the cluster of maximum membership
fcm$quality

## -----------------------------------------------------------------------------
maxmemb <- apply(fcm$solution$membership, 1, max)
sum(maxmemb < 0.6)

## ----results = "hide"---------------------------------------------------------
set.seed(1)
afcm <- WH_adaptive_fcmeans(OzoneFull, k = 3, schema = 2, m = 1.6, rep = 3,
                            simplify = TRUE, qua = 20, standardize = TRUE)

## -----------------------------------------------------------------------------
round(afcm$solution$var.weights[, 1], 3)
afcm$quality

## ----fig.height = 4-----------------------------------------------------------
hc <- WH_hclust(OzoneFull, simplify = TRUE, qua = 20, standardize = TRUE,
                method = "ward.D2")
plot(hc, labels = FALSE, main = "OzoneFull, Ward method")
table(hclust = cutree(hc, k = 3), kmeans = km$solution$IDX)

