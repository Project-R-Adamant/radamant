\name{plotkit}
\alias{draw.grid}
\alias{draw.legend}
\alias{draw.x.axis}
\alias{draw.x.title}
\alias{draw.y.title}
\alias{draw.y.axis}
\title{Plotting Tools}
\description{Utility functions used for Plotting}

\usage{ 
draw.grid(X
	, base = NULL
	, theme.params = getCurrentTheme()
	, method = c("equispaced", "sampling")
	) 

draw.legend(legend = ""
	, theme.params = getCurrentTheme()
	, overrides = list(...)
	, ...
	) 

draw.x.axis(X
	, base = NULL
	, xlabels = NULL
	, theme.params = getCurrentTheme()
	, show.labels = TRUE
	, show.ticks = TRUE
	, ...
	) 
	
draw.x.title(xtitle = "", theme.params = getCurrentTheme())

draw.y.axis(X
	, ylabels = NULL
	, theme.params = getCurrentTheme()
	, side = 1
	, show.labels = TRUE
	, show.ticks = TRUE
	, ...
	) 
	
draw.y.title(ytitle = "", theme.params = getCurrentTheme(), side = 1)

}

\arguments{
\item{X}{Matrix of data series being plotted (y-values). One column per series.}
\item{base}{Corresponding x-values (common to all series) associated to the entries of X. If NULL, then base = 1:NROW(X). }
\item{theme.params}{A valid RAdamant Theme. See setThemeAttr for details. (DEFAULT = getCurrentTheme())}
\item{overrides}{List of parameters used to override the theme. Only parameters that match those defined by the theme are overridden (DEFAULT = list(...))}
\item{legend}{Vector of legend texts}
\item{xlabels}{Labels for the x-axis}
\item{ylabels}{Labels for the y-axis}
\item{xtitle}{Title for the x-axis}
\item{ytitle}{Title for the y-axis}
\item{show.labels}{LOGICAL. If TRUE, labels are showed.}
\item{show.ticks}{LOGICAL. If TRUE, tickmarks are showed.}
\item{side}{The side (1 = left, 2 = right) where the y-axis labels and title are plotted.} 
\item{method}{Controls how the x-coordinates of the grid vertical lines are computed. \cr
If method = "equispaced", N = getThemeAttr("x.ticks", exact = TRUE) points between min(base) and max(base) are computed. \cr
If method = "sampling", the N lines are drawn at the points given by base[seq(1, length(base), len = N)].
}
\item{...}{ Further arguments to or from other methods.} 
}

\details{
These are utility funtions used as building blocks for high level plotting with cplot.
Most of the behaviour is controlled by the theme options. \cr
For details on the available options, see \code{\link{setThemeAttr}}.
}


\value{Void}

\seealso{\code{\link{cplot}}}

\author{ RAdamant Development Team \email{team@r-adamant.org} }

