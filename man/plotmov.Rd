\name{plotmov}

\alias{plot.Movav}

\title{Plot Moving Average}

\description{Plot method for object of class 'Movav' (Moving Average)}

\usage{
\method{plot}{Movav}(x, fs = NULL, main = attr(x, "desc"), ...) 
}

\arguments{ 
\item{x}{ instance of class 'Movav'}
\item{fs}{ Matrix containing the original data series (one column per variable).  \cr 
For financial time series (class = 'fs'), only 'Close' column is processed.}
\item{main}{ Main title of the plot}
\item{...}{ Additional parameters accepted by the functions cplot and fin.plot}
}

\details{If the original data series is an instance of class 'fs', then the plot will have two
panels:
    \itemize{
    \item{plot of fs and x on the top;}
    \item{histogram of the Volume data of the financial series X.}
    }
}

\value{VOID}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{ \code{\link{cplot}} }

\examples{

# Compute Exponential Moving Average and plot results
x = ema(rnorm(100), 10)

# Plot Multiple Moving Averages together using "" plotting class 
plot(x)

## load a dataset provided by R 
data(EuStockMarkets)
# extract sample (log) time series 
x = log(EuStockMarkets[500:1000,2, drop=FALSE])
# set RAdamant theme (1 - Finance or 2 - Vanilla)
setCurrentTheme(1)
plot.Movav(cbind(kama(x)
                ,frama(x)
                ,ema(x, 10)
                ,gdema(x, 10)
                ,zlma(x, 10)
                ) 
           , x 
           )

# plot multiple moving average results from an object of class "fs"
data(ex_fs)
class(ex_fs)
x = ex_fs
# set RAdamant theme (1 - Finance or 2 - Vanilla)
setCurrentTheme(2)
plot.Movav(cbind(kama(x)
                ,frama(x)
                ,ema(x, 10)
                ,gdema(x, 10)
                ,zlma(x, 10)
                ) 
           , x 
           )

}



