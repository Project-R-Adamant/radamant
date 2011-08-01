\name{gartest}

\alias{Archlm}
\alias{LjungBox}

\title{Garch residual tests}

\description{Compute ARCH-LM and Ljung-Box test for residual correlation}

\usage{
Archlm(x, lags, std=FALSE, plot.acf=FALSE) 
LjungBox(x, lags, plot.acf = FALSE)
}

\arguments{ 

\item{x}{ Series of residual or an object of class "Garch".}
\item{lags}{ Number of lags to calculate the autocorrelation function.}
\item{plot.acf}{ Logical. If TRUE plot of autocorrelation function is returned.}
\item{std}{ Logical. If TRUE input residual will be standardised.}
}

\seealso{ \code{\link{Garch}}  }

\author{RAdamant Development Team \email{team@r-adamant.org}}


\examples{
# load example time series 
data(ex_ts)
x = ex_ts

gg1 = Garch(x, order = c(1,1), type="garch", prob="norm")
# perform Ljung-Box test with 10 lags
LjungBox(gg1, 10)
# perform ARCH-LM test with 10 lags and show ACF plot
Archlm(gg1, 1, std=TRUE, plot.acf=TRUE)
}


