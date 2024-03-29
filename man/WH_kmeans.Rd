% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unsuperv_classification.R
\name{WH_kmeans}
\alias{WH_kmeans}
\title{K-means of a dataset of histogram-valued data}
\usage{
WH_kmeans(
  x,
  k,
  rep = 5,
  simplify = FALSE,
  qua = 10,
  standardize = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{x}{A MatH object (a matrix of distributionH).}

\item{k}{An integer, the number of groups.}

\item{rep}{An integer, maximum number of repetitions of the algorithm (default \code{rep}=5).}

\item{simplify}{A logic value (default is FALSE), if TRUE histograms are recomputed in order to speed-up the algorithm.}

\item{qua}{An integer, if \code{simplify}=TRUE is the number of quantiles used for recodify the histograms.}

\item{standardize}{A logic value (default is FALSE). If TRUE, histogram-valued data are standardized,  variable by variable,
using the Wassertein based standard deviation. Use if one wants to have variables with std equal to one.}

\item{verbose}{A logic value (default is FALSE). If TRUE, details on computations are shown.}
}
\value{
a list with the results of the k-means of the set of Histogram-valued data \code{x} into  \code{k} cluster.
}
\description{
The function implements the k-means for a set of histogram-valued data.
}
\section{Slots}{

\describe{
\item{\code{solution}}{A list.Returns the best solution among the \code{rep}etitions, i.e.
the one having the minimum sum of squares criterion.}

\item{\code{solution$IDX}}{A vector. The clusters at which the objects are assigned.}

\item{\code{solution$cardinality}}{A vector. The cardinality of each final cluster.}

\item{\code{solution$centers}}{A \code{MatH} object with the description of centers.}

\item{\code{solution$Crit}}{A number. The criterion (Sum od square deviation
from the centers) value at the end of the run.}

\item{\code{quality}}{A number. The percentage of Sum of square deviation explained by the model.
(The higher the better)}
}}

\examples{
results <- WH_kmeans(
  x = BLOOD, k = 2, rep = 10, simplify = TRUE,
  qua = 10, standardize = TRUE, verbose = TRUE
)
}
\references{
Irpino A., Verde R., Lechevallier Y. (2006). Dynamic clustering of histograms using Wasserstein
metric. In: Rizzi A., Vichi M.. COMPSTAT 2006 - Advances in computational statistics. p. 869-876,
Heidelberg:Physica-Verlag
}
