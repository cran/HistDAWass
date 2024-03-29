% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Met_MatH.R
\name{WH.mat.prod}
\alias{WH.mat.prod}
\alias{WH.mat.prod,MatH,MatH-method}
\alias{WH.mat.prod,MatH-method}
\title{Method WH.mat.prod}
\usage{
WH.mat.prod(object1, object2, ...)

\S4method{WH.mat.prod}{MatH,MatH}(object1, object2, traspose1 = FALSE, traspose2 = FALSE)
}
\arguments{
\item{object1}{a \code{MatH} object}

\item{object2}{a \code{MatH} object}

\item{...}{other optional parameters}

\item{traspose1}{a logical value, default=FALSE. If TRUE trasposes object1}

\item{traspose2}{a logical value, default=FALSE. If TRUE trasposes object2}
}
\value{
a matrix of numbers
}
\description{
It is the matrix product of two \code{MatH} objects, i.e. two matrices of distributions,
by using the dot product of two histograms that is consistent with
 a set of distributions equipped with a L2 wasserstein metric.
}
\examples{

M1 <- BLOOD[1:5, ]
M2 <- BLOOD[6:10, ]
MAT <- WH.mat.prod(M1, M2, traspose1 = TRUE, traspose2 = FALSE)
}
