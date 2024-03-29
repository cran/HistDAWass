% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/All_classes.R, R/Met_MatH.R
\docType{class}
\name{MatH-class}
\alias{MatH-class}
\alias{MatH}
\alias{initialize,MatH-method}
\title{Class MatH.}
\usage{
\S4method{initialize}{MatH}(
  .Object,
  nrows = 1,
  ncols = 1,
  ListOfDist = NULL,
  names.rows = NULL,
  names.cols = NULL,
  by.row = FALSE
)

MatH(
  x = NULL,
  nrows = 1,
  ncols = 1,
  rownames = NULL,
  varnames = NULL,
  by.row = FALSE
)
}
\arguments{
\item{.Object}{the object type "MatH"}

\item{nrows}{(optional, default=1)an integer, the number of rows.}

\item{ncols}{(optional, default=1) an integer, the number of columns (aka
variables).}

\item{ListOfDist}{a vector or a list of \code{distributionH} objects}

\item{names.rows}{a vector or list of strings with thenames of the rows}

\item{names.cols}{a vector or list of strings with thenames of the columns (variables)}

\item{by.row}{(optional, default=FALSE) a logical value, TRUE the matrix is
row wise filled, FALSE the matrix is filled column wise.}

\item{x}{(optional, default= an empty \code{distributionH} object) a list of
\code{distributionH} objects}

\item{rownames}{(optional, default=NULL) a list of strings containing the
names of the rows.}

\item{varnames}{(optional, default=NULL) a list of strings containing the
names of the columns (aka variables).}
}
\value{
A \code{matH} object
}
\description{
Class \code{MatH} defines a matrix of \code{distributionH} objects

This function create a matrix of histogram data, i.e. a \code{MatH}
object
}
\examples{

## ---- create a list of six distributionH objects
ListOfDist <- vector("list", 6)
ListOfDist[[1]] <- distributionH(c(1, 2, 3), c(0, 0.4, 1))
ListOfDist[[2]] <- distributionH(c(7, 8, 10, 15), c(0, 0.2, 0.7, 1))
ListOfDist[[3]] <- distributionH(c(9, 11, 20), c(0, 0.5, 1))
ListOfDist[[4]] <- distributionH(c(2, 5, 8), c(0, 0.3, 1))
ListOfDist[[5]] <- distributionH(c(8, 10, 15), c(0, 0.75, 1))
ListOfDist[[6]] <- distributionH(c(20, 22, 24), c(0, 0.12, 1))

## create a MatH object filling it by columns
MyMAT <- new("MatH",
  nrows = 3, ncols = 2, ListOfDist = ListOfDist,
  names.rows = c("I1", "I2", "I3"), names.cols = c("Var1", "Var2"), by.row = FALSE
)

showClass("MatH")

# bulding an empty 10 by 4 matrix of histograms
MAT <- MatH(nrows = 10, ncols = 4)
}
\references{
Irpino, A., Verde, R. (2015) \emph{Basic
statistics for distributional symbolic variables: a new metric-based
approach} Advances in Data Analysis and Classification, DOI
10.1007/s11634-014-0176-4
}
\author{
Antonio Irpino
}
\keyword{classes}
