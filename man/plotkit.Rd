\name{plotkit}
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
draw.grid(X, base = NULL, theme.params = getCurrentTheme()) 

draw.legend(legend = "", theme.params = getCurrentTheme(), 
overrides = list(...), ...) 

draw.projections(X, Y, Y.fit, 
col = getCurrentTheme()[["projection.col"]][1],
type = getCurrentTheme()[["projection.type"]][1], 
lty = getCurrentTheme()[["projection.lty"]][1]) 

draw.x.axis(X, base = NULL, xlabels = NULL, 
theme.params = getCurrentTheme(), show.labels = TRUE) 

draw.x.title(xtitle = "", theme.params = getCurrentTheme()) 

draw.y.axis(X, ylabels = NULL, theme.params = getCurrentTheme(), 
side = 1, show.labels = TRUE) 

draw.y.title(ytitle = "", theme.params = getCurrentTheme(), side = 1) 
}

\arguments{
\item{X}{ X}
\item{Y}{ Y}
\item{base}{ base}
\item{theme.params}{ theme.params}
\item{overrides}{ overrides}
\item{legend}{ legend}
\item{xlabels}{ xlabels}
\item{ylabels}{ ylabels}
\item{xtitle}{ xtitle}
\item{ytitle}{ ytitle}
\item{show.labels}{ show.labels}
\item{Y.fit}{ Y.fit}  
\item{col}{ col} 
\item{type}{ type}  
\item{lty}{ lty}  
\item{side}{ side} 
\item{...}{ Further arguments to or from other methods.} 
}

\author{ RAdamant Development Team \email{team@r-adamant.org} }
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
