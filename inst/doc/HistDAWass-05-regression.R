## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 5, dpi = 72, out.width = "100%")

## ----setup--------------------------------------------------------------------
library(HistDAWass)

## -----------------------------------------------------------------------------
get.MatH.main.info(Agronomique)$varnames
pars <- WH.regression.two.components(Agronomique, Yvar = 1, Xvars = 2:4)
pars

## -----------------------------------------------------------------------------
pred <- WH.regression.two.components.predict(Agronomique[, 2:4], parameters = pars)
pred[1:3, ]

## -----------------------------------------------------------------------------
gof <- WH.regression.GOF(observed = Agronomique[, 1], predicted = pred)
gof$RMSE_W
gof$OMEGA
gof$PSEUDOR2$index

## ----fig.height = 8-----------------------------------------------------------
plotPredVsObs(pred[1:8, ], Agronomique[1:8, 1], type = "HISTO", ncolu = 2)

## ----fig.height = 4-----------------------------------------------------------
plotPredVsObs(pred[1:4, ], Agronomique[1:4, 1], type = "CDF", ncolu = 2)

## ----fig.height = 4-----------------------------------------------------------
plot_errors(pred, Agronomique[, 1], type = "HISTO_QUA")

## -----------------------------------------------------------------------------
pars20 <- WH.regression.two.components(Agronomique, Yvar = 1, Xvars = 2:4,
                                       simplify = TRUE, qua = 20)
round(rbind(exact = pars, simplified = pars20), 4)

