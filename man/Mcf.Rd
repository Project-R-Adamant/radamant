\name{mcf}
\alias{mcf}

\title{Auto-Correlation and Partial Auto-Correlation}

\description{Compute auto-correlation and partial auto-correlation function on a matrix}

\usage{
mcf(X, lag.max = 10, ci = 0.95, plot=TRUE, ...)
}
\arguments{ 
\item{X}{ Matrix of data series (one column per variable)}
\item{lag.max}{ Max lag to be computed by the cross correlation function (DEFAULT: 10)}
\item{ci}{Confidence Interval (DEFAULT: 0.95)}
\item{plot}{LOGICAL. If TRUE, results are plotted.}
\item{...}{Additional parameters accepted by the function plot.cross.ccf.}
}

\value{
An object of class "mcf". This is a list with two entries:
\item{ACF}{list of Auto-Correlation Functions (one for each column of X).}
\item{PACF}{list of Partil Auto-Correlation Functions (one for each column of X).}
}

\examples{
\dontrun{
# Dow Jones
DJ = get.fs("^DJI", SName = "DowJones", from=as.Date("2008-06-01"), to=as.Date("2009-04-01"));
# Compute Returns
RDJ = Ret(DJ, na.rm = TRUE)

# Plot Autocorrelation Function and Partial ACF 
mcf(RDJ, lag.max = 30) 
# Using another theme
mcf(RDJ, lag.max = 30, theme = getTheme("vanilla")) 
}
}

\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{ \code{\link{cross.ccf}} }


