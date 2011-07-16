\name{mcf}
\alias{mcf}

\title{Auto-Correlation and Partial Auto-Correlation}

\description{Compute auto-correlation and partial auto-correlation function on a matrix}

\usage{mcf(X, lag.max = 10, ci = 0.95, plot=TRUE, ...) }
\arguments{ 
\item{X}{ Matrix of data series (one column per variable)}
\item{lag.max}{ Max lag to be computed by the cross correlation function (DEFAULT: 10)}
\item{ci}{ Confidence Interval (DEFAULT: 0.95)}
\item{plot}{ LOGICAL. If TRUE, results are plotted.}
\item{...}{ additional parameters accepted by the function plot.cross.ccf.}
}

\value{ A list with two entries:
\item{ACF}{list of Auto-Correlation Functions (one for each column of X)}
\item{PACF}{list of Partil Auto-Correlation Functions (one for each column of X)}
}

\examples{
\dontrun{
# Plot Autocorrelation Function and Partial ACF 
mcf(RSP, lag.max = 30) 
# using another theme
theme = getTheme("Vanilla")
mcf(RDJ, lag.max = 30, theme = getTheme(2)) 
}

}

\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{ \code{\link{cross.ccf}} }


