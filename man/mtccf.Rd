\name{Method_ccf}

\alias{print.cross.ccf}
\alias{plot.cross.ccf}

\title{Cross.ccf functions}

\description{Methods for class 'cross.ccf'}

\usage{
\method{print}{cross.ccf}(x, ...) 
\method{plot}{cross.ccf}(x, theme.params = getCurrentTheme(), xtitle = "Lag", ytitle =
expression(rho), overrides = list(...), ...) 
}

\arguments{ 
\item{x}{ Instance of class 'cross.ccf'}
\item{theme.params}{ Theme parameters (DEFAULT: getCurrentTheme())}
\item{xtitle}{ Title for the x-axis (DEFAULT: "Lag")}
\item{ytitle}{ Title for the y-axis (DEFAULT: expression(rho))}
\item{overrides}{ List of parameters to override the theme. Must match by name the parameters
defined by the theme (DEFAULT: NULL)}
\item{...}{ Further arguments to or fromother methods}
}

\value{Void}

\author{ RAdamant Development Team \email{team@r-adamant.org} }
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
