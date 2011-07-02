\name{hvar}

\alias{hVaR}

\title{Historical Value at Risk}

\description{Compute historical VaR on each column of the input matrix}

\usage{hVaR(X, p = 0.05, centered = FALSE)}

\arguments{
\item{X}{ Input matrix/sequence. Sequences are treated as one column matrices.}
\item{p}{ vector of probabilities (DEFAULT = 0.05)}
\item{centered}{ LOGICAL. If TRUE, input data are standardised}
}

\value{A matrix length(p) by NCOL(X) of computed quantiles}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\note{TO BE COMPLETED}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
