\name{funlcnt}
\alias{func.line.cnt}
\alias{plot.modularity}

\title{Modularity Analysis}

\description{
Given a package name or a list of functions, for each function X
in the package or the list it counts the lines of code, the number of subcalls
made to any other function Y of the list/package and the number of other
functions that make calls to the function X. \cr 
Results are plotted if requested.
}

\usage{
func.line.cnt(package = NULL, plot = TRUE, ...) 

\method{plot}{modularity}(x
	, qtz.type = "linear"
	, qtz.nbins = 30
	, qtz.cutoff = 30
	, theme.params = getCurrentTheme()
	, overrides = list(...)
	, border = "transparent"
	, savepng = FALSE
	, savepath = getwd()
	, save.width = 480
	, save.height = 480
	, save.resolution = 72
	, ...
	)
}

\arguments{ 
\item{x}{ Package name or a list of functions.}
\item{package}{CHARACTER. Single name of the package to load or array list of function names.}
\item{plot}{LOGICAL. If TRUE, results are plotted on bar charts.}
\item{qtz.type}{CHARACTER. qtz.type = "Linear" | "Log" | "None". Partial match on the value is attempted.}
\item{qtz.nbins}{INTEGER. Number of bins to be computed. Used only when qtz.type is "Linear" or "Log" (Default = 30).}
\item{qtz.cutoff}{Used only when qtz.type = "Log" (Default = 30). More granular binning below the cutoff point.}
\item{theme.params}{A valid RAdamant Theme. See setThemeAttr for details. (DEFAULT = getCurrentTheme())}
\item{overrides}{List of parameters used to override the theme. Only parameters that match those defined by the theme are overridden (DEFAULT = list(...))}
\item{border}{Color used for the border line of the barplot.}
\item{savepng}{LOGICAL. If true, charts are saved to png file.}
\item{savepath}{The path where png files are saved (DEFAULT = getwd()).}
\item{save.width}{The image width of the png file. See png for details.}
\item{save.height}{The image height of the png file. See png for details.}
\item{save.resolution}{The image resolution of the png file. See png for details.}
\item{...}{ Further arguments to or from other methods.} 
}

\details{
The parameter "qtz.type" controls the type of quantization used to set the bin size for the bar chart of the Code Length Distribution.
Values:
\itemize{
	\item{If "Linear", qtz.nbins equispaced intervarls are computed.}
	\item{If "Log", qtz.nbins log-spaced intervals are computed based on qtz.cutoff.}
	\item{In any other case the bin size is set to 1.}
}
The parameter "qtz.cutoff" controls how bins are computed when qtz.type = "Log": qtz.nbins equispaced intervals are computed on a log(x/qtz.cutoff) scale. 
This creates more intervals/bins in the range 0 < x < qtz.cutoff.

}

\value{
An object of the class "modularity". 
This is a data frame containing the stats for each function in the input list/package, with the following columns:
\item{fcn.name}{Name of the function.}
\item{fcn.lines}{Number of lines of code.}
\item{fcn.subcalls}{Number of distinct calls made to other functions.}
\item{fcn.called}{Number of distinct functions using this function.}

The following attribute is attached to the object:
\item{package}{The input package argument.}
}

\examples{
\dontrun{
# Modularity Analysis for the RAdamant package
rad = func.line.cnt(package = "RAdamant")
# Log quantization
plot(rad, qtz.type = "Log", qtz.cutoff = 10)
}
}
\author{RAdamant Development Team \email{team@r-adamant.org}}



