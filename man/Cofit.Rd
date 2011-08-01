\name{cofit}

\alias{cofit}

\title{Cornish Fisher Transformation}

\description{Estimate quantiles based on Cornish Fisher formula, which only uses skewness and kurtosis.}

\usage{
cofit(X, p, k = NULL, s = NULL)
}

\arguments{
\item{X}{Input matrix/sequence. Sequences are treated as one column matrices.}
\item{p}{Vector of probability threshold (interval [0, 1])}
\item{k}{Kurtosis (Default: NULL -> becomes kurt(X))}
\item{s}{Skewness (Default: NULL -> becomes skew(X))}
}

\value{A matrix length(p) by NCOL(X) of estimated quantiles.}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
# Estimate 5% quantile from Normal random data.
cofit(rnorm(1000), p = 0.05)
# Compare to theoretical quantile
qnorm(0.05)

# Estimate 5% quantile from Student's T random data.
cofit(rt(1000, 16), p = 0.05)
# Compare to theoretical quantile
qt(0.05, df = 16)
}


