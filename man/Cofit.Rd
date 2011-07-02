\name{cofit}

\alias{cofit}

\title{Cornish Fisher Transformation}

\description{Cornish Fisher Transformation}

\usage{
cofit(X, p, k = NULL, s = NULL)
}

\arguments{
\item{X}{ Input matrix/sequence. Sequences are treated as one column matrices.}
\item{p}{ vector of probability threshold (interval [0, 1])}
\item{k}{ kurtosis (DEFAULT = NULL -> becomes kurt(X))}
\item{s}{ skewness (DEFAULT = NULL -> becomes skew(X))}
}

\value{ A matrix length(trsh) by NCOL(X) of computed quantiles}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\note{TO BE COMPLETED}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
