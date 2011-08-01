\name{movapply}
\alias{movApply}
\title{Moving Apply function}
\description{Applies a given function to a sliding window of the input data}
\usage{
movApply(X
    , win.size = 1
    , padding = NA
    , rm.transient = FALSE
    , func = NULL
    , ...
    ) 
}
\arguments{ 
\item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ Vector of data window sizes that will be passed to the given function "func" (Default: 1). }
\item{padding}{ Padding value to fill transient of result (output data rows from 1 to win.size-1). Default: NA.}
\item{rm.transient}{ transient: Logical. If TRUE, transient is removed, otherwise funct is applied to the transient (Default: FALSE).}
\item{func}{ Function to be run.}
\item{...}{ Additional parameters accepted by the function func.}
}
\details{For financial time series (class = 'fs'), only 'Close' column is processed.}

\value{
A matrix of size NROW(X) by NCOL(X)*length(win.size). \cr
The function specified by 'func' is applied to each sliding window SWi (given by win.size[i]) and each column of X.}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
