\name{lagret}

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

Ret(X
    , lag = 1
    , log = FALSE
    , mode = "selected"
    , na.rm = FALSE
    , plot = FALSE
    , ...) 

Lag(X
    , lag = 1
    , na.rm = FALSE
    , padding = NA) 

Diff(X
    , lag = 1
    , padding = NA
    , na.rm = FALSE) 

MDiff(X
      , lag = 1
      , padding = NA
      , mode = c("auto", "range", "selected")
      , na.rm = FALSE) 

MLag(X
    , lag = 1
    , na.rm = FALSE
    , padding = NA
    , mode = c("auto", "range", "selected")
    , autolag.start = 1) 
}

\arguments{ 
\item{X}{Input data (i.e. matrix/vector of prices)}
\item{lag}{INTEGER or VECTOR. number of lags (it can be both positive and negative)}
\item{log}{BOOLEAN: compute log-returns }
\item{na.rm}{BOOLEAN: remove NAs }
\item{plot}{BOOLEAN: return plot }
\item{padding}{value to replace removed observations}
\item{mode}{ mode of using the vector of lags}
\item{autolag.start}{ autolag.start}
\item{...}{ Further arguments to or from other methods}
}

\details{Sequences are treated as one-column matrices. \cr
The parameter "mode" allows to control the calculation when the parameter is passed as a vector:
\itemize{
\item{auto}: only the first element is used;
\item{range}: if the lag arguments is composed of two numbers, the computation is performed for all \cr
the integers contained in the interval, ex: lag = c(4,10) allow to calculate all the lags between 4 and 10; 
\item{selected}: the computation is done only for the lag specified in the argument.
}

}

\value{ A matrix (n.obs X n.lag) containing lagged /differenced time series or returns }


\examples{

# load an example dataset containing financial daily prices
data(ex_fs)
x = ex_fs[ ,1:4]

# compute multiple multiple lags for single time series 
# different uses of the parameter "mode"
res = MLag(x[,1], lag = c(4,8), mode="range")
res[1:10, ]
res = MLag(x[,1], lag = c(4,8), mode="selected")
res[1:10, ]
res = MLag(x[,1], lag = 4, mode="auto")
res[1:10, ]

## SINGLE LAG
# calculate return for single time series 
res = Ret(x[ ,1], lag=4, log=TRUE, na.rm=TRUE)
res[1:10, ,drop=FALSE]

# calculate return for multiple time series
res = Ret(x, lag=10, log=TRUE, na.rm=TRUE)
res[1:10, ,drop=FALSE]

## MULTIPLE LAGS
# calculate return for single time series 
res = Ret(x[ ,1], lag=c(2,4,6,8), mode = "selected", log=TRUE, na.rm=TRUE)
res[1:10, ,drop=FALSE]

# calculate return for multiple time series
res = Ret(x[, 1:2], lag=c(2,4,6,8), mode = "selected", log=FALSE, na.rm=FALSE)
res[1:10, ,drop=FALSE]


## PLOT RESULTS
# calculation and plot for single series
Ret(x[,1]
    , lag = 5
    , mode = "selected"
    , plot=TRUE
    , style="bar"
    , main="Returns - 5 Lags")

# calculation and plot for multiple series
par(mfrow=c(2,2))
Ret(x
    , lag = 5
    , mode = "selected"
    , plot=TRUE
    , style="bar"
    , main="Returns - 5 Lags")

\dontrun{
# get APPLE financial series 
symbol.lookup("Apple") 
APPLE = get.fs("AAPL", from=as.Date("2008-06-01"), to=as.Date("2011-04-01"));
RAPPLE = Ret(APPLE
                , mode = "selected"
                , plot = TRUE
                , style = "bar"
                ,  ylab.fmt = .3
                , na.rm = TRUE)
RAPPLE;
}

}

\seealso{ \code{\link{plot.ret}} }


\author{RAdamant Development Team \email{team@r-adamant.org}}


