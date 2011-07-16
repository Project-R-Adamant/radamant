\name{plotroi}
\alias{plot.roi}
\title{Plot Return on Investment objects}
\description{Plot method for class 'roi'.}
\usage{
\method{plot}{roi}(x, main = "Historical Return on Investment", xtitle = "Lag", ...) 
}

\arguments{
\item{x}{Instance of class 'roi'.}
\item{main}{Title for the plot.}
\item{xtitle}{The title for the x-axis.}
\item{...}{Additional parameters passed to the cplot function.}
}

\value{Void}

\seealso{
\code{\link{cplot}}.
}

\examples{
# Load sample financial series data
data(ex_fs)

# Analyse the performance of the returns (Close data) up to 200 days and plot results
plot(hroi(ex_fs, lag = 200, log = FALSE), xlab.srt = 0)

# Analyse the performance of the returns (All data) up to 200 days and plot results
plot(hroi(ex_fs[,], lag = 200, log = FALSE), xlab.srt = 0)

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
