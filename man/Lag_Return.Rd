\name{Lag_Return}
\alias{Ret}
\alias{Lag}
\alias{MLag}
\alias{MDiff}
\alias{Diff}
\title{Time Series Operators}
\description{
\code{Ret}: Compute N-points Returns on each column of the input matrix. \cr 
\code{Lag}: Compute lag on each column of the input matrix. \cr
\code{Diff}: Compute lagged difference on each column of the input matrix. \cr
\code{MDiff}: Compute Multiple lagged differences on each column of the input matrix. \ cr
\code{MLag}: Compute Multiple lags on each column of the input matrix
}
\usage{

Ret(X, lag = 1, log = FALSE, na.rm = FALSE, plot = FALSE, ...) 

Lag(X, lag = 1, na.rm = FALSE, padding = NA) 

Diff(X, lag = 1, padding = NA, na.rm = FALSE) 

MDiff(X, lag = 1, padding = NA, 
mode = c("auto", "range", "selected"), na.rm = FALSE) 

MLag(X, lag = 1, na.rm = FALSE, padding = NA, 
mode = c("auto", "range", "selected"), autolag.start = 1) 
}
\arguments{ \item{X}{Input data (i.e. matrix/vector of prices)}
\item{lag}{number of lags }
\item{log}{BOOLEAN: compute log-returns }
\item{na.rm}{BOOLEAN: remove NAs }
\item{plot}{BOOLEAN: return plot }
\item{padding}{value to replace removed observations}
\item{mode}{ mode}
\item{autolag.start}{ autolag.start}
\item{...}{ Further arguments to or from other methods}
}
\details{Sequences are treated as one-column matrices.}
\author{RAdamant Development Team \email{team@r-adamant.org}}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
