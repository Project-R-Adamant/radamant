\name{whes}

\alias{whES}

\title{Weighted Historical Expected Shortfall}

\description{
Compute Weighted historical ES on each column of the input matrix.
If input is a Financial Series object (class 'fs'), then 'Close' data are processed.
}

\usage{whES(X, p = 0.05, lambda = 0.9, centered = FALSE)}

\arguments{
\item{X}{Input matrix/sequence. Sequences are treated as one column matrices.}
\item{p}{Vector of probabilities (Default: 0.05).}
\item{lambda}{Controls the exponential window lambda^((NROW(X)-1):0) (Default: 0.9).}
\item{centered}{Logical. If TRUE, input data are standardised.}
}

\value{A matrix length(p) by NCOL(X) of computed historical weighted ES.}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
# Load sample Financial series data
data(ex_fs);
# Compute Historical Weighted ES (5% confidence level) on 1-day Returns
whES(Ret(ex_fs));

# Generate some random data
X = cbind(rnorm(1000), rnorm(1000, sd = 2))
# Compute multiple Historical Weighted ES (1%, 2.5%, 5% confidence levels)
whES(X, p = c(1, 2.5, 5)/100);
}


