\name{getacfci}
\alias{get.acf.ci}
\title{Normal confidence intervals for correlation}
\description{Compute the Normal confidence intervals for correlation and partial autocorrelation data.}
\usage{
get.acf.ci(X, ci = 0.95) 
}
\arguments{
\item{X}{Instance of class 'acf' as returned by functions acf, pacf, ccf}
\item{ci}{Confidence interval required (DEFAULT: 0.95)}
}
\value{
A vector containing the two symmetrical confidence intervals.
}

\examples{
# Generate some random integrated data
x = cumsum(rnorm(30));
# The confidence intervals in mcf are calculated using get.acf.ci
res = mcf(x)
# Extract the PACF part and compute the confidence intervals
get.acf.ci(res$PACF[[1]])
# Same as
get.acf.ci(pacf(x, plot = FALSE))
}
\author{RAdamant Development Team \email{team@r-adamant.org}}



