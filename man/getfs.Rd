\name{getfs}
\alias{get.fs}

\title{Download Financial Series data from Yahoo!}

\description{Download Yahoo! time series data and returns a Financial Series (fs) object.}

\usage{
get.fs(symbol = NULL 
	, SName = NULL 
	, from = as.Date("1950-01-01") 
	, to = Sys.Date() 
	, strip.spaces = TRUE 
	, strip.char = "." 
	)
}

\arguments{
\item{symbol}{Stock symbol to download.}
\item{SName}{Name that will be assigned to the time series. If NULL (default) the name is retrieved from Yahoo!}
\item{from}{Date object. The start date of the time series (DEFAULT = as.Date("1950-01-01")).}
\item{to}{Date object. The end date of the time series (DEFAULT = Sys.Date()).}
\item{strip.spaces}{LOGICAL. If TRUE, spaces from SName are replaced with the value of strip.char (DEFAULT = TRUE).}
\item{strip.char}{The character used to replaces spaces in SName (DEFAULT = ".").}
}

\value{
A financial Time Series object. This is a matrix of Yahoo! daily data with columns (Open, High, Low, Close, Volume, Adj.Close).
The following attributes are attached to the object:
\item{SName}{The Name/Description of the financial series.}
\item{Symbol}{the input stock symbol.}
}

\examples{

# Get Dow Jones quotes from Jan 2010
DowJones = get.fs("^DJI", from = as.Date("2010-01-01"))
DowJones

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
