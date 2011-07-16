\name{ptfoper}

\alias{PtfBeta}
\alias{PtfRet}
\alias{PtfVar}
\title{Portfolio operators}

\description{Get portfolio Beta}

\usage{
PtfRet(PTF, w = NULL, glob = TRUE, calc.ret = FALSE, ...) 

PtfVar(PTF, w = NULL, glob = TRUE, 
vol = FALSE, calc.ret = FALSE, ...)

PtfBeta(beta, w = NULL, glob = TRUE) 
}

\arguments{
\item{PTF}{ Matrix containing one or more series of prices/returns, one time series for each asset}
\item{w}{ Vector of portfolio weights}
\item{glob}{ Logical. If TRUE return the value for the whole portfolio.}
\item{vol}{ Logical. If TRUE returns volatility (standard deviation instead of variance).}
\item{calc.ret}{ Logical. If TRUE the input matrix is considered as a matrix of prices, so returns are calculated.}
\item{beta}{ Value of the Beta coefficient or an object of class "Capm".}
\item{...}{ Further arguments to or from other methods.}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
# load example portfolio
data(ex_ptf)
# results for each series 
PtfRet(ex_ptf, glob=FALSE)
PtfVar(ex_ptf, glob=FALSE)
# results for the whole portfolio
PtfRet(ex_ptf, glob=TRUE)
PtfVar(ex_ptf, glob=TRUE)

# Example with a series of prices instead of returns
data(EuStockMarkets)
PtfRet(PTF = EuStockMarkets, w=c(0.3, 0.4, 0.2, 0.1), calc.ret=TRUE)
PtfRet(PTF = EuStockMarkets, w=c(0.3, 0.4, 0.2, 0.1), glob = FALSE, calc.ret=TRUE)

}



