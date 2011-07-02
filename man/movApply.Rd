\name{movapply}
\alias{movApply}
\title{Moving Apply function}
\description{Applies a given function to a sliding window of the input data}
\usage{movApply(X, win.size = 1, padding = NA, rm.transient = FALSE, func = NULL, ...) }
\arguments{ 
\item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of data window sizes that will be passed to the given function "func"
(DEFAULT = 1). }
\item{padding}{ Padding value to fill transient of result (output data rows from 1 to win.size-1).
(DEFAULT = NA)}
\item{rm.transient}{ transient: LOGICAL. If TRUE, transient is removed, otherwise func is applied to
the transient. (DEFAULT
= FALSE)}
\item{func}{ Function to be run}
\item{...}{ Additional parameters accepted by the function func}
}
\details{For financial time series (class = 'fs'), only 'Close' column is processed.}
\value{A matrix of size NROW(X) by NCOL(X)*length(win.size). func is applied to each sliding window
SWi (given by
win.size[i]) and each column of X.}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
