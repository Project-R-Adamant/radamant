\name{state_space}
\alias{ss.sym}
\title{State Space system simulation}
\description{Generic function for State Space system simulation. The system can be either linear or
non linear.}
\usage{ss.sym(X, F = NULL, G = NULL, H = NULL, D = NULL, 
init = 0, SLen = ifelse(is.function(F), NA,
NROW(F)), YLen = ifelse(is.function(H), NA, NROW(H)), ...) }
\arguments{ 
\item{X}{ Matrix of data series (one column per variable).}
\item{F}{ [State -> State] transition matrix or [(State, Input) -> State] function (F = function(S,
X, n, ...) \cr
 returning the new state vector S_new based on the current State S and the data X at time period n)
(DEFAULT = NULL)}
\item{G}{ [Input -> State] transition matrix. Only for linear models (DEFAULT = NULL)}
\item{H}{ [State -> Output] transition matrix or [(State, Input) -> Output] function (H =
function(S, X, n, ...) \cr
 returning the new output vector Y[, n] based on the new state S[, n] and the data X at time period
n) (DEFAULT = NULL -> converted in diag(SLen)) }
\item{D}{ [Input -> Output] transition matrix. Only for linear models (DEFAULT = NULL -> converted
to a zero matrix SLen by NCOL(X) )}
\item{init}{ Initial values for the state vactor S (DEFAULT = 0, recycled to length SLen if
necessary)}
\item{SLen}{ Length of the state vector S. (DEFAULT = ifelse(is.function(F), NA, NROW(F)) )}
\item{YLen}{ Number of columns of the output vector Y. (DEFAULT = ifelse(is.function(H), NA,
NROW(H)) )}
\item{...}{ Additional parameters accepted by the functions F and H}
}
\details{For financial time series (class = 'fs'), only 'Close' column is processed.}
\value{ A object of class 'ss' with attributes 'F', 'G', 'H', 'D' as given by the corresponding
input parameters: \cr
	- matrix of size NROW(X) by YLen, result of the symulation of the given dynamic system subject to
input 'X' and initial condition 'init'.
}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
