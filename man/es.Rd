\name{es}

\alias{ES}
\alias{ES.default}

\title{Expected Shortfall}

\description{
General ES, computed on each column of the input matrix.
If input is a Financial Series object (class 'fs'), then 'Close' data are processed.
}

\usage{
ES(X, ...)

\method{ES}{default}(X
    , p = 0.05
    , probf = c("Normal", "T-Student", "Cornish-Fisher", "GPD-POT")
    , df = max(4, (kurt(X)+3))
    , trsh = -hVaR(X)
    , ...
    )
}

\arguments{
\item{X}{Input matrix/sequence. Sequences are treated as one column matrices.}
\item{p}{Vector of probabilities (Default: 0.05)}
\item{probf}{Probability distribution (see details). Case insensitive, partial matching is supported.}
\item{df}{Degrees of freedom for the Student T distribution (Default: max(4, (kurt(X)+3)))}
\item{trsh}{vector of NCOL(X) thresholds used to identify the tail data for the GPD-POT method}
\item{...}{Additional parameters passed to the functions 'cofit' and 'gpd.ES'.}
}

\details{
Accepted probability distributions:
\itemize{
    \item{"Normal": Normal distribution.}
    \item{"T-Student": Student'T distribution.}
    \item{"Cornish-Fisher": Cornish-Fischer formula for quantiles estimation.}
    \item{"GPD-POT": Peak Over Threshold method, based on Generalised Pareto Distribution (EVT).}
}
}

\value{A matrix length(p) by NCOL(X) of computed ES values, based on the input distribution.}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\seealso{\code{\link{gpd.ES}}, \code{\link{mqt}}, \code{\link{cofit}}.}

\examples{
# Load sample asset data
data(ex_ptf);
# Compute ES on multiple confidence levels (Normal)
ES(ex_ptf, p = seq(0.03, 0.05, by = 0.01), probf = "Normal");

# T-Student
ES(ex_ptf, p = seq(0.03, 0.05, by = 0.01), probf = "T");

# Extreme Value Theory (GPD)
ES(ex_ptf, p = seq(0.03, 0.05, by = 0.01), probf = "GPD");

}