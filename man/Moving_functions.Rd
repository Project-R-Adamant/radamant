\name{Moving_functions}

\alias{movMax}
\alias{movMin}
\alias{movSd}
\alias{movVar}


\title{Moving Base Functions}

\description{Applies the function "Max", "Min", "Standard Deviation" or "Variance" to a sliding window of the input data}

\usage{movMax(X, win.size = 1, ...) \cr
movMin(X, win.size = 1, ...) \cr
movSd(X, win.size = 1, ...) \cr
movVar(X, win.size = 1, ...) }

\arguments{
\item{X}{ Matrix of data series (one column per variable). }
\item{win.size}{ vector of data window sizes that will be used for the calculations (DEFAULT = 1).}
\item{...}{ Additional parameters accepted by the function movApply}
}

\details{For financial time series (class = 'fs'), only 'Close' column is processed.}

\value{A matrix of size NROW(X) by NCOL(X)*length(win.size). max is applied to each sliding window SWi (given by win.size[i])
and}

\references{}

\author{ RAdamant Development Team \email{team@r-adamant.org} }

\note{}

\seealso{ \code{\link{movApply}} }

\examples{}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
