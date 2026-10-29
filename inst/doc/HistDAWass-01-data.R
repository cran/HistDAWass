## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 4.5, dpi = 72, out.width = "100%")

## ----setup--------------------------------------------------------------------
library(HistDAWass)

## -----------------------------------------------------------------------------
d1 <- distributionH(x = c(1, 2, 3), p = c(0, 0.4, 1))
d1

## -----------------------------------------------------------------------------
d1@m
d1@s

## -----------------------------------------------------------------------------
distributionH(5)@x
d2 <- distributionH(x = c(0, 1, 2, 4), p = c(0, 0.5, 0.5, 1))
diff(d2@p) # the probabilities of the bins

## ----fig.height = 3.5---------------------------------------------------------
plot(d1)

## -----------------------------------------------------------------------------
set.seed(1)
x <- c(rnorm(300, mean = 10, sd = 2), rnorm(200, mean = 18, sd = 1))

h_base <- data2hist(x, algo = "base")                     # the bins of hist()
h_quan <- data2hist(x, algo = "FixedQuantiles", qua = 10) # 10 equal-probability bins
h_man  <- data2hist(x, algo = "ManualBreaks", breaks = seq(0, 25, by = 2.5))
h_pl   <- data2hist(x, algo = "PolyLine", epsilon = 0.01) # piecewise linear cdf
sapply(list(base = h_base, quantiles = h_quan, manual = h_man, polyline = h_pl),
       function(h) length(h@x) - 1) # number of bins

## ----fig.height = 3.5---------------------------------------------------------
plot(h_quan)
plot(h_quan, type = "CDF")

## -----------------------------------------------------------------------------
ListOfDist <- list(
  distributionH(c(1, 2, 3), c(0, 0.4, 1)),
  distributionH(c(7, 8, 10, 15), c(0, 0.2, 0.7, 1)),
  distributionH(c(9, 11, 20), c(0, 0.5, 1)),
  distributionH(c(2, 5, 8), c(0, 0.3, 1)),
  distributionH(c(8, 10, 15), c(0, 0.75, 1)),
  distributionH(c(20, 22, 24), c(0, 0.12, 1))
)
MyMAT <- MatH(x = ListOfDist, nrows = 3, ncols = 2,
              rownames = c("I1", "I2", "I3"), varnames = c("Var1", "Var2"))
MyMAT

## -----------------------------------------------------------------------------
groups <- levels(iris$Species)
vars <- names(iris)[1:4]
cells <- list()
for (v in vars) {           # the matrix is filled by columns
  for (g in groups) {
    cells <- c(cells, data2hist(iris[iris$Species == g, v],
                                algo = "FixedQuantiles", qua = 5))
  }
}
IRIS <- MatH(x = cells, nrows = length(groups), ncols = length(vars),
             rownames = groups, varnames = vars)
IRIS

## ----fig.height = 4-----------------------------------------------------------
plot(IRIS)                  # histograms
plot(IRIS, type = "DENS")   # smoothed densities

## ----fig.height = 4-----------------------------------------------------------
p <- plot(IRIS, type = "BOXPLOT")
p + ggplot2::ggtitle("Iris data: boxplots of the histograms")

## -----------------------------------------------------------------------------
dim(BLOOD@M)
B <- BLOOD[1:4, c("Cholesterol", "Hemoglobin")]
B
get.MatH.main.info(B)

WH.bind.row(BLOOD[1:2, ], BLOOD[13:14, ])

## -----------------------------------------------------------------------------
get.cell.MatH(BLOOD, r = 1, c = 1)
B2 <- set.cell.MatH(distributionH(c(100, 150, 200), c(0, 0.5, 1)), B, r = 1, c = 1)
B2@M[[1, 1]]

## -----------------------------------------------------------------------------
R <- registerMH(MyMAT)
is.registeredMH(MyMAT)
is.registeredMH(R)
R@M[[1, 1]]@p

## ----fig.height = 3.5---------------------------------------------------------
RetHTS
plot(subsetHTS(RetHTS, from = 1, to = 15), type = "BOXPLOT")

