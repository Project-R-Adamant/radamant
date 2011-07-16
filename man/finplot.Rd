\name{finplot}
\alias{fin.plot}
\title{Plot financial time series}
\description{Generic plotting for financial data. Produces a two panels plot}

\usage{
fin.plot(X 
	, top.vars = c("Close", "High", "Low") 
	, bottom.vars = "Volume" 
	, style = c("default", "candlestick") 
	, snames = attr(X, "SName") 
	, xlabels = rownames(X) 
	, main = "" 
	, main2 = "" 
	, ytitle = "" 
	, ytitle2 = "" 
	, theme.top = getCurrentTheme() 
	, overrides = list(...) 
	, theme.bottom = getCurrentTheme() 
	, overrides2 = NULL 
	, ... 
	)
}

\arguments{
\item{X}{Input matrix of data to be plotted.}
\item{top.vars}{Indices or names of the columns for the top plot.}
\item{bottom.vars}{Indices or names of the columns for the bottom plot.}
\item{style}{Not used. For future releases.}
\item{snames}{Names of the series being plotted.}
\item{xlabels}{labels for the x-axis.}
\item{main}{Main title for the top plot.}
\item{main2}{Main title for the bottom plot.}
\item{ytitle}{Title for the y-axis (top plot).}
\item{ytitle2}{Title for the y-axis (bottom plot).}
\item{theme.top}{Theme parameters list for the top plot (DEFAULT = getCurrentTheme()).}
\item{overrides}{List of parameters to override theme for the top plot. Only parameters that match those defined by the theme are overridden (DEFAULT}{list(...)).}
\item{theme.bottom}{Theme parameters list for the bottom plot.}
\item{overrides2}{List of parameters to override theme for the bottom plot. (DEFAULT = NULL).}
\item{...}{Additional parameters passed to the cplot function. Also used to quickly specify theme overrides.}
}

\seealso{
\code{\link{cplot}}.
}

\examples{
# Load sample financial series data
data(ex_fs)
# Plot the data
plot(ex_fs)
# Change the style and color of the bottom chart
plot(ex_fs, overrides2 = list(type = "l", col = "grey"))
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
