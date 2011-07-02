\name{plotmov}

\alias{plot.Movav}

\title{Plot Moving Average}

\description{Plot method for object of class 'Movav' (Moving Average)}

\usage{
\method{plot}{Movav}(x, fs = NULL, main = attr(x, "desc"), ...) 
}

\arguments{ 
\item{x}{ instance of class 'Movav'}
\item{fs}{ Matrix containing the original data series (one column per variable).  For financial time
series (class = 'fs'), only 'Close' column is processed.}
\item{main}{ Main title of the plot}
\item{...}{ Additional parameters accepted by the functions cplot and fin.plot}
}

\details{If the original data series is an instance of class 'fs', then the plot will have two
panels:\cr
	- plot of fs and x on the top; \cr
	- histogram of the Volume data of the financial series X.
}

\value{VOID}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{ \code{\link{cplot}} }

\examples{
# Compute Exponentioal Moving Average and plot results
x = ema(rnorm(100), 10)
plot(x)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
