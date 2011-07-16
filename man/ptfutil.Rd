\name{ptfutil}

\alias{PtfUtility}

\title{Portfolio Utility}

\description{Calculate utility and plot for efficient portfolio}

\usage{PtfUtility(PTF = NULL, W, R = NULL, SIGMA = NULL, 
af = 3, plot = TRUE, ...) }

\arguments{
\item{PTF}{ Matrix containing TWO series of returns, one series for each asset.}
\item{W}{ Initial vector of weights.}
\item{R}{ Vector of PTF returns.}
\item{SIGMA}{ PTF sample covariance matrix.}
\item{af}{ Numeric (range: 0,1). Adversion factor (Default = 3)}
\item{plot}{ LOGICAL. If TRUE plot is returned.}
\item{...}{Further arguments to or from other methods.}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\seealso{\code{\link{PtfFront}}, \code{\link{PtfOpt}} }

\examples{
# vector of returns for two assets A and B
R = c(A=0.021, B=0.09)
# Covariance matrix
SIGMA = matrix(c(0.101^2, 0.005,  0.005, 0.208^2),2,2)
# Calculate and show utility for the two assets
PtfUtility(PTF=NULL, R=R, SIGMA=SIGMA, W=c(0.4,0.6))
}



