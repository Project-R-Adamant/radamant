\name{whvar}

\alias{whVaR}

\title{Weighted Historical Value at Risk}

\description{Compute Weighted historical VaR on each column of the input matrix}

\usage{whVaR(X, p = 0.05, lambda = 0.9, centered = FALSE)}

\arguments{
\item{X}{ Input matrix/sequence. Sequences are treated as one column matrices.}
\item{p}{ vector of probabilities (DEFAULT = 0.05)}
\item{lambda}{ controls the exponential window lambda^((NROW(X)-1):0) (DEFAULT = 0.9)}
\item{centered}{ LOGICAL. If TRUE, input data are standardised}
}

\value{A matrix length(p) by NCOL(X) of computed quantiles}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\note{TO BE COMPLETED}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
