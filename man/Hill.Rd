\name{hill}

\alias{Hill}

\title{Hill function}

\description{Approximation of the shape parameter (xi) of the Generalised Pareto distribution.}

\usage{
Hill(X, trsh = hVaR(X))
}

\arguments{
\item{X}{Input matrix/sequence. Sequences are treated as one column matrices.}
\item{trsh}{Vector of NCOL(X) thresholds used to identify the tail data for the estimation.}
}

\value{ A matrix 1 by NCOL(X) of computed shape parameters}
\author{RAdamant Development Team \email{team@r-adamant.org}}



