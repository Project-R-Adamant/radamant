\name{drawdown}

\alias{drawdown.default}
\alias{drawdown}
\alias{summary.drawdown}
\alias{ExtremeDD}

\title{Drawdown}

\description{Drawdown risk analysis}

\usage{
drawdown(x, ...) 
\method{drawdown}{default}(x
                                    , FUN=max
                                    , relative=FALSE
                                    , plots=c("regular", "smooth", "no.plot")
                                    , ...) 
\method{summary}{drawdown}(object, show.extr=TRUE, ...)
ExtremeDD(DD, FUN, lag = 1, rolling = FALSE, plot = TRUE, ...)
}

\details{The function "ExtremeDD" is called inside "summary.drawdown".}

\arguments{ 
\item{x}{ Univariate input time series.}
\item{FUN}{ Extreme function applied for the max / min drawdown calculation (Default: max)}
\item{relative}{ Logical. If TRUE relative drawdown will be calculated.}
\item{plots}{ Character. Type of plot to be returned (De)}
\item{DD, object}{ An object of class "drawdown"}
\item{show.extr}{ Logical. if TRUE extreme drawdown will be calculated. }
\item{lag}{ Integer. Number of lag periods used for rolling calculation.}
\item{rolling}{ Logical. If TRUE extreme will be calculated on a moving window.}
\item{plot}{ Logical. If TRUE plot is returned.}
\item{...}{ Further arguments accepted by the function \code{\link{cplot}} or \code{\link{sma}}.}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
# load example time series
data(ex_ts)
x = ex_ts

# calculate drawdown - no plot
dd = drawdown(x, plots="no.plot")
# calculate drawdown - regular plot
dd = drawdown(x, plots="regular")
# calculate drawdown - smoothed plot with different color
dd = drawdown(x, plots="smooth", col="green")

# summary information and maximum drawdown
summary(dd)
# ... summary information and rolling maximum drawdown
summary(dd, rolling=TRUE, lag=10)

}


