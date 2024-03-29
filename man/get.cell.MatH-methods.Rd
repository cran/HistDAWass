% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Met_MatH.R
\name{get.cell.MatH}
\alias{get.cell.MatH}
\alias{get.cell.MatH,MatH,numeric,numeric-method}
\alias{get.cell.MatH,MatH-method}
\title{Method get.cell.MatH Returns the histogram in a cell of a matrix of distributions}
\usage{
get.cell.MatH(object, r, c)

\S4method{get.cell.MatH}{MatH,numeric,numeric}(object, r, c)
}
\arguments{
\item{object}{a MatH object, a matrix of distributions.}

\item{r}{an integer, the row index.}

\item{c}{an integer, the column index}
}
\value{
A \code{distributionH} object.
}
\description{
Returns the histogram data in the r-th row and the c-th column.
}
\examples{
get.cell.MatH(BLOOD, r = 1, c = 1)
}
