\name{asfs}
\alias{as.fs}
\title{Convert Yahoo! Data into Financial Series object}
\description{Converts a stock data series (dataframe) into a Financial Series (fs) object.}
\usage{as.fs(X, SName = "", Symbol = "") }

\arguments{
\item{X}{Input dataframe with columns (Open, High, Low, Close, Volume, Adj.Close).}
\item{SName}{The name assigned to the fs object.}
\item{Symbol}{The symbol assigned to the fs object.}
}

\value{
A financial Time Series object. This is a matrix with columns (Open, High, Low, Close, Volume, Adj.Close).\cr
The following attributes are attached to the object:
\item{SName}{The Name/Description of the financial series.}
\item{Symbol}{The input stock symbol.}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}


\examples{

# Load sample financial series data
data(ex_fs)
# Subset data and create another fs object
as.fs(as.data.frame(ex_fs[1:10,])
    , SName = "My Financial Series"
    , Symbol = "My Symbol"
    )

}
