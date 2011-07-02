\name{funlcnt}
\alias{func.line.cnt}

\title{Function line counting}

\description{Given a package name or a list of functions, for each function X
in the package or the list it counts the lines of code, the number of subcalls
made to any other function Y of the list/package and the number of other
functions that make calls to the function X. \cr 
Results are plotted if requested.}

\usage{func.line.cnt(package = NULL, plot = TRUE, 
qtz.type = "NONE", qtz.nbins = 10, qtz.cutoff = 30) }

\arguments{ 
\item{package}{ CHARACTER. Single name of the package to load or array list of function names.}
\item{plot}{ LOGICAL. If TRUE, results are plotted on bar charts.}
\item{qtz.type}{ CHARACTER. qtz.type = "NONE" | "LINEAR" | "LOG".}
\item{qtz.nbins}{ INTEGER. Number of bins to be computed. Used only when qtz.type != "NONE".
(Default = 10)}
\item{qtz.cutoff}{ Used only when qtz.type = "LOG". (Default = 30)}
}

\details{
Parameter "qtz.type" is Case Insensitive. It states the yype of quantizzation to be used to set bin
size for the barchart plotting the distribution of lines of code.
Values:
	\itemize{
	\item{If "NONE", bin size is set to 1.}
	\item{If "LINEAR", qtz.nbins equispaced intervarls are computed.}
	\item{If "LOG", qtz.nbins log-spaced intervals are computed based on qtz.cutoff.}
	}
Parameter "qtz.bins": qtz.nbins equispaced intervals are computed on a log(x/qtz.cutoff) scale. This
creates more intervals in the range 0 < x < qtz.cutoff.

}

\value{
Data frame containing the stats for each function in the input list/package:
	\itemize{
	\item{fcn.name = Name of the function}
	\item{fcn.lines = Number of lines of code}
	\item{fcn.subcalls = Calls made to other functions}
	\item{fcn.called = Number of function calling the function}
	}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
