\name{scaledf}
\alias{scalApply}
\alias{scalMin}
\alias{scalMax}
\title{Apply functions on a scaled window}
\description{
\code{scalApply}: Applies a given function to the pairs (X[n, i], X[n-lag, i]). \cr
\code{scalMax}: Scaled max on each column of the input matrix.
\code{scalMin}: Scaled min on each column of the input matrix
}
\usage{
scalApply(X
        , lag = 0
        , padding = NA
        , na.rm = FALSE
        , func = NULL
        , ...
        )
        
scalMax(X
        , lag = 1
        , padding = -Inf
        , na.rm = FALSE
        , func = NULL
        ) 
        
scalMin(X
        , lag = 1
        , padding = Inf
        , na.rm = FALSE
        , func = NULL
        ) 
}
\arguments{ 
\item{X}{ Input matrix/sequence}
\item{lag}{ vector of integer lags. If lag >= 0 data are shifted to the right, else to the left.(Default: 0)}
\item{padding}{ value used to initialise the output matrix (Default: NA)}
\item{na.rm}{ Logical. If TRUE, N-lag entries are removed from the output (Default: FALSE)}
\item{func}{ function applied to the data (Default: NULL)}
\item{...}{ Additional parameters accepted by the function 'func'}
}

\details{Sequences are treated as one-column matrices.}

\value{
A matrix where func/max/min has been applied on each pair (X[n, i], X[n-lag,i]) for each column i of X. \cr
Number of rows depends on the na.rm parameter. Number of columns is NCOL(X)
}
\author{RAdamant Development Team \email{team@r-adamant.org} }


