\name{Plot_MA}

\alias{plot.ma}

\title{Plot Moving Average}

\description{Plot method for object of class 'ma' (Moving Average)}

\usage{
plot.ma(MovAv, X = NULL, main = attr(MovAv, "desc")) 
}

\arguments{ 
\item{MovAv}{ instance of class 'ma'}
\item{X}{ Matrix containing the original data series (one column per variable).  For financial time series (class = 'fs'), only 'Close' column is processed.}
\item{main}{ Main title of the plot}
\item{...}{ Additional parameters accepted by the functions cplot and fin.plot}
}

\details{If the original data series is an instance of class 'fs', then the plot will have two panels:\cr
	- plot of X and MovAv on the top; \cr
	- histogram of the Volume data of the financial series X.
}

\value{VOID}

\references{}

\author{ RAdamant Development Team \email{team@r-adamant.org} }

\note{}

\seealso{ \code{\link{cplot}} }

\examples{
# Compute Exponentioal Moving Average and plot results
x = ema(rnorm(100), 10)
plot(x)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
