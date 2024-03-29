% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/principal_components.R
\name{WH.MultiplePCA}
\alias{WH.MultiplePCA}
\title{Principal components analysis  of a set of histogram variable based on Wasserstein distance}
\usage{
WH.MultiplePCA(data, list.of.vars, quantiles = 10, outl = 0)
}
\arguments{
\item{data}{A MatH object (a matrix of distributionH).}

\item{list.of.vars}{A list of  integers, the active variables.}

\item{quantiles}{An integer, it is the number of quantiles used in the analysis.
Default=10.}

\item{outl}{a number between 0 (default)  and 0.5. For each distribution, is the amount of mass removed from the tails
of the distribution. For example, if 0.1, from each distribution is cut away a left tail and a right one each containing
the 0.1 of mass.}
}
\value{
a list with the results of the PCA in the MFA format of package \pkg{FactoMineR} for function MFA
}
\description{
(Beta version) The function implements a Principal components analysis of a set of histogram variables
based on Wasserstein distance. It performs a centered (not standardized) PCA on a set of quantiles of a variable.
Being a distribution a multivalued description, the analysis performs a dimensional reduction and a visualization of distributions.
It is a 1d (one dimension) becuse it is considered just one histogram variable.
}
\details{
It is an extension of WH.1d.PCA to the multiple case.
}
