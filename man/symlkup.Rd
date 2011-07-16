\name{symlkup}
\alias{symbol.lookup}

\title{Lookup Stock Symbol from Yahoo!}

\description{Lookup stock symbols for which the symbol, name or description matches the input string value.}

\usage{symbol.lookup(what = "")}

\arguments{
\item{what}{The string to search for.}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\value{ 
A matrix containing the top 10 stock symbols that match the input, with the following columns:
\item{Symbol}{The stock symbol.}
\item{Name}{The stock name.}
\item{Exchange}{The Exchange symbol.}
\item{Type}{The Exchange Name.}
}


\examples{
# lookup the symbol for Apple 
symbol.lookup("Apple") 
# Apple 
APPLE = get.fs("AAPL", from=as.Date("2008-06-01"), to=as.Date("2011-04-01"));
}

\seealso{\code{\link{get.fs}}}
