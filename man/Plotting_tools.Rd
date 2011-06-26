\name{Plotting_tools}
\alias{draw.grid}
\alias{draw.legend}
\alias{draw.projections}
\alias{draw.x.axis}
\alias{draw.x.title}
\alias{draw.y.title}
\alias{draw.y.axis}

\title{Plotting Tools}


\description{Utilities functions used for Plotting}

\usage{ 
draw.grid(X, base = NULL, theme.params = getCurrentTheme()) \cr
draw.legend(legend = "", theme.params = getCurrentTheme(), overrides = list(...), ...) \cr
draw.projections(X, Y, Y.fit, col = getCurrentTheme()[["projection.col"]][1]) , \cr
type = getCurrentTheme()[["projection.type"]][1], lty = getCurrentTheme()[["projection.lty"]][1]) \cr
draw.x.axis(X, base = NULL, xlabels = NULL, theme.params = getCurrentTheme(), show.labels = TRUE) \cr
draw.x.title(xtitle = "", theme.params = getCurrentTheme()) \cr
draw.y.axis(X, ylabels = NULL, theme.params = getCurrentTheme(), side = 1, show.labels = TRUE) \cr
draw.y.title(ytitle = "", theme.params = getCurrentTheme(), side = 1) 
}

\arguments{
\item{X}{}
\item{Y}{}
\item{base}{}
\item{theme.pars}{}
\item{overrides}{}
\item{legend}{}
\item{xlabels}{}
\item{ylabels}{}
\item{xtitle}{}
\item{ytitle}{}
\item{show.labels}{}
}

\details{ }

\value{ }

\references{ }

\author{ RAdamant Development Team \email{team@r-adamant.org} }

\note{ }

\seealso{ }

\examples{ }

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
