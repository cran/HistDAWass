% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Met_HTS.R
\docType{methods}
\name{plot-TdistributionH}
\alias{plot-TdistributionH}
\alias{plot,TdistributionH-method}
\title{plot for a TdistributionH object}
\usage{
\S4method{plot}{TdistributionH}(x, type = "HISTO", col = "green", border = "black")
}
\arguments{
\item{x}{a \code{TdistributionH} object}

\item{type}{(optional) a string describing the type of plot, default="HISTO".\cr Other allowed types are
\cr"CDF"=Cumulative distribution function, \cr"QF"= quantile function, \cr"DENS"=a density approximation,
\cr"HBOXPLOT"=horizontal boxplot, \cr"VBOXPLOT"= vertical boxplot,}

\item{col}{(optional) a string the color of the plot, default="green".}

\item{border}{(optional) a string the color of the border of the plot, default="black".}
}
\description{
A plot function for a \code{TdistributionH} object. The function returns a representation
of the histogram.
}
