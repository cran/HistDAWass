## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 5, dpi = 72, out.width = "100%")

## ----setup--------------------------------------------------------------------
library(HistDAWass)

## ----results = "hide", fig.height = 7-----------------------------------------
pca1 <- WH.1d.PCA(BLOOD, var = 1, quantiles = 10, listaxes = c(1, 2))

## -----------------------------------------------------------------------------
round(head(pca1$PCAout$eig, 3), 2)
round(head(pca1$PCAout$ind$coord[, 1:2]), 3)

## -----------------------------------------------------------------------------
round(pca1$PCAout$quanti.sup$cor[, 1:2], 3)

## -----------------------------------------------------------------------------
mpca <- WH.MultiplePCA(BLOOD, list.of.vars = 1:3, quantiles = 10, graph = FALSE)
round(head(mpca$eig, 3), 2)

## -----------------------------------------------------------------------------
WH.plot_multiple_Spanish.funs(mpca, axes = c(1, 2), var = 1:3)

## ----fig.height = 5-----------------------------------------------------------
WH.plot_multiple_indivs(BLOOD, mpca, axes = c(1, 2), var = 1,
                        HISTO = TRUE, stat = "mean")

## ----fig.height = 5-----------------------------------------------------------
WH.plot_multiple_indivs(BLOOD, mpca, axes = c(1, 2), var = 2,
                        strx = 0.15, stry = 0.2, HISTO = FALSE)

