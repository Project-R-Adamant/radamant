\name{cumfun}

\alias{cumMax}
\alias{cumMin}
\alias{cumMean}
\alias{cumSd}
\alias{cumVar}
\alias{cumSum}
\title{Cumulative functions}
\description{
Cumulative max / min / Mean / Standard Deviation / Variance / sum on each column of the input
matrix.
}
\usage{cumMax(X, lag = 0, padding = NA, na.rm = FALSE) }
\arguments{ 
\item{X}{ Input matrix/sequence}
\item{lag}{ vector of integer lags. If lag >= 0 data are shifted to the right, else to the left.
(DEFAULT = 0)}
\item{padding}{ value used to initialise the output matrix (DEFAULT = NA)}
\item{na.rm}{ LOGICAL. If TRUE, N-lag entries are removed from the output. Also NA in the input are
replaced by -Inf (DEFAULT = FALSE)}
}
\details{Sequences are treated as one-column matrices}
\value{A matrix of cumulative maximums of X. Number of rows depends on the na.rm parameter. Number
of columns is NCOL(X)}
\author{RAdamant Development Team \email{team@r-adamant.org} }
\seealso{ \code{\link{lew}} }


