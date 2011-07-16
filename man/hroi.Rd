\name{hroi}
\alias{hroi}
\title{Historical Returns on Investment}
\description{
Computes historical returns on investment and two-sided VaR. Analysis of the performance of the returns as a function of the holding period.
For Financial series objects (class 'fs'), Close data is processed.
}

\usage{
hroi(X 
	, lag = 1 
	, mode = c("auto", "range", "selected") 
	, autolag.start = 1 
	, range.step = 1 
	, log = TRUE 
	, VaR.type = "norm" 
	, p = 0.05 
	, ... 
	)
}

\arguments{h
\item{X}{Input matrix of data to be plotted.}
\item{lag}{The maximum lag used to compute returns (DEFAULT = 1).}
\item{mode}{Controls how the lags are computed. See details.}
\item{autolag.start}{Starting lag value for the case where mode = "auto" (DEFAULT = 1). See details.}
\item{range.step}{Lag increment used for the case where mode = "range" (DEFAULT = 1). See details.}
\item{log}{LOGICAL. If TRUE, log returns are computed. DEFAULT = TRUE.}
\item{VaR.type}{The distribution used for VaR calculation. See \code{\link{VaR}} for details.}
\item{p}{The confidence interval used for VaR calculation. (DEFAULT = 0.05)}
\item{...}{Additional parameters passed to the VaR function.}
}

\value{
An instance of the class 'roi'. This is a list of length given by the number of columns of the input X.
Each entry is a matrix with columns [Return (Avg.), VaR (Profit), VaR (Loss)] where the rows are calculated for each lag.
The following attributes are attached to the object:
\item{log}{The input log parameter.}
\item{lag}{The lags for which returns are computed.}
}

\details{
For each input time series, returns are calculated for multiple lags, hence average and two-sided Value at Risk (Profit & Loss with p% confidence interval) are computed on the returns.
The number and the way lags are computed is controlled by the mode parameter:
\itemize{
\item{auto: All lags between autolag.start and max(lag) (DEFAULT option)}
\item{range: All lags between min(lag) and max(lag) with increment given by range.step}
\item{selected: Only selected lags are calculated.}
}
}

\seealso{
\code{\link{Ret}}, \code{\link{VaR}}, \code{\link{plot.roi}}.
}

\examples{
# Load sample financial series data
data(ex_fs)

# Historical returns for all lags between 1 and 10 days
hroi(ex_fs, lag = 10)

# Historical returns for lags between 2 and 10 with increment 2
hroi(ex_fs, lag = c(2, 10), mode = "range", range.step = 2)

# Historical returns for selected lags
hroi(ex_fs, lag = c(2, 5, 10), mode = "selected")

# Analyse the performance of the returns up to 200 days and plot results
plot(hroi(ex_fs, lag = 200, log = FALSE), xlab.srt = 0)

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
