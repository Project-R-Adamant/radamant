\name{croscf}

\alias{cross.ccf}

\title{Cross Correlation Function}

\description{Compute the cross correlation function for each pairs of variables (Yi Xj)}

\usage{
cross.ccf(Y, X, lag.max = 10, ci = 0.95, plot = TRUE, ...)
}


\arguments{ 
\item{Y}{Matrix of data series (one column per variable)}
\item{X}{Matrix of data series (one column per variable)}
\item{lag.max}{Max lag to be computed by the cross correlation function (Default: 10)}
\item{ci}{Confidence Interval (Default: 0.95)}
\item{plot}{Logical. If TRUE, results are plotted.}
\item{...}{Additional parameters accepted by the function plot.cross.ccf.}
}

\value{
An object of class "cross.acf". This is a list of Ny*Nx elements, where each entry is the cross correlation of the pair (Yi, Xj).
}

\examples{
# Generate two random integrated series
N = 100
X = matrix(rnorm(N), nrow = N/2, ncol=2);
# Create two series as a linear combination of X plus noise
Y = X %*% diag(c(1, -1)) + 0.1*rnorm(N);
# Perform Cross Correlation Analysis
cross.ccf(Y, X)
}

\author{ RAdamant Development Team \email{team@r-adamant.org} }


