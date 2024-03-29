% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unsuperv_classification.R
\name{WH_MAT_DIST}
\alias{WH_MAT_DIST}
\title{L2  Wasserstein distance matrix}
\usage{
WH_MAT_DIST(x, simplify = FALSE, qua = 10, standardize = FALSE)
}
\arguments{
\item{x}{A MatH object (a matrix of distributionH).}

\item{simplify}{A logic value (default is FALSE), if TRUE histograms are recomputed in order to speed-up the algorithm.}

\item{qua}{An integer, if \code{simplify}=TRUE is the number of quantiles used for recodify the histograms.}

\item{standardize}{A logic value (default is FALSE). If TRUE, histogram-valued data are standardized,  variable by variable,
using the Wasserstein based standard deviation. Use if one wants to have variables with std equal to one.}
}
\value{
A matrix of squared L2 distances.
}
\description{
The function extracts the L2 Wasserstein distance matrix from a MatH object.
}
\examples{
DMAT <- WH_MAT_DIST(x = BLOOD, simplify = TRUE)
}
\references{
Irpino A., Verde R. (2006). A new Wasserstein based distance for the hierarchical clustering
of histogram symbolic data. In: Batanjeli et al. Data Science and Classification, IFCS 2006. p. 185-192,
 BERLIN:Springer, ISBN: 3-540-34415-2
}
