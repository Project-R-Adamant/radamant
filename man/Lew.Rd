\name{lew}
\alias{lew}
\title{Moving window}
\description{Apply a given function to an extending window of the lagged data series of the input
matrix, each column separately.}
\usage{lew(X, lag = 0, padding = NA, na.rm = FALSE, 
func = NULL, is.cumulative = TRUE, ...) }
\arguments{ 
\item{X}{ Input matrix/sequence}
\item{lag}{ vector of integer lags. If lag >= 0 data are shifted to the right, else to the left.
(Default: 0)}
\item{padding}{ value used to initialise the output matrix (Default: NA)}
\item{na.rm}{Logical. If TRUE, N-lag entries are removed from the output (Default: FALSE)}
\item{func}{ function applied to the extending data window (Default: NULL)}
\item{is.cumulative}{Logical. If TRUE it the function provided must be cumulative by itself (like
cummax, cummin, etc..) (Default: TRUE)}
\item{...}{ Additional parameters accepted by the function 'func'}
}
\details{ Sequences are treated as one-column matrices}
\value{A matrix where func has been applied on increasing data windows for each column of X. Number
of rows depends on the na.rm parameter. Number of columns is NCOL(X)}
\author{RAdamant Development Team \email{team@r-adamant.org} }
\seealso{ \code{\link{cumSum}}, \code{\link{cumMin}}, 
\code{\link{cumMax}}, \code{\link{cumSd}}, \code{\link{cumVar}}}


