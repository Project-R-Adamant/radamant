\name{movfunc}
\alias{movMax}
\alias{movMin}
\alias{movSd}
\alias{movVar}
\title{Moving Base Functions}
\description{Applies the function "Max", "Min", "Standard Deviation" or "Variance" to a sliding
window of the input data}
\usage{movMax(X, win.size = 1, ...) 
movMin(X, win.size = 1, ...) 
movSd(X, win.size = 1, ...) 
movVar(X, win.size = 1, ...) }
\arguments{
\item{X}{ Matrix of data series (one column per variable). }
\item{win.size}{ Vector of data window sizes that will be used for the calculations (DEFAULT = 1).}
\item{...}{ Additional parameters accepted by the function movApply}
}
\details{For financial time series (class = 'fs'), only 'Close' column is processed.}
\value{A matrix of size NROW(X) by NCOL(X)*length(win.size). max is applied to each sliding window
SWi (given by win.size[i])
and}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{ \code{\link{movApply}} }


