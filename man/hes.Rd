\name{hes}

\alias{hES}

\title{Historical Expected Shortfall}

\description{
Compute historical ES on each column of the input matrix. 
If input is a Financial Series object (class 'fs'), then 'Close' data are processed.
}

\usage{hES(X, p = 0.05, centered = FALSE)}

\arguments{
\item{X}{Input matrix/sequence. Sequences are treated as one column matrices.}
\item{p}{vector of probabilities (Default: 0.05)}
\item{centered}{Logical. If TRUE, input data are standardised prior to compute ES.}
}

\value{A matrix length(p) by NCOL(X) of computed historical VaR}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
# Load sample Financial series data
data(ex_fs);
# Compute Historical ES (5% confidence level) on 1-day Returns
hES(Ret(ex_fs));

# Generate some random data
X = cbind(rnorm(1000), rnorm(1000, sd = 2))
# Compute multiple Historical ES (1%, 2.5%, 5% confidence levels)
hES(X, p = c(1, 2.5, 5)/100);
}


