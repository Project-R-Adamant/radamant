\name{varptf}

\alias{VaRPtf}

\title{Portfolio Value at Risk}

\description{General VaR, computed for an input portfolio}

\usage{
VaRPtf(X, p = 0.05, weights = rep(1/NCOL(X), NCOL(X)), ...)
}

\arguments{
\item{X}{Input matrix/sequence. Sequences are treated as one column matrices.}
\item{p}{Vector of probabilities (Default: 0.05)}
\item{weights}{Portfolio weigths (Default: rep(1/NCOL(X), NCOL(X)))}
\item{...}{Additional parameters passed to the 'VaR' function}
}

\value{A matrix length(p) by 1 of computed portfolio VaR values.}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\seealso{\code{\link{VaR}}.}

\examples{
# Load sample asset data
data(ex_ptf);
# Compute VaR on multiple confidence levels (GPD)
VaRPtf(ex_ptf[, -1], p = seq(0.01, 0.05, by = 0.01), probf = "GPD");
}


