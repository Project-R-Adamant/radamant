\name{MethodMcf}

\alias{print.mcf}
\alias{plot.mcf}

\title{Mcf methods}

\description{Plot and Print method for class 'mcf'}

\usage{
\method{print}{mcf}(x, ...) 

\method{plot}{mcf}(x, theme.params = getCurrentTheme(), 
xtitle = "Lag", ytitle = expression(rho),
overrides = NULL, ...)
}

\arguments{ 
\item{x}{ Instance of class 'mcf'}
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
