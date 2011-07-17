\name{mtmcf}

\alias{print.mcf}
\alias{plot.mcf}

\title{Multi-Correlation Function methods}

\description{Plot and Print method for class 'mcf'}

\usage{
\method{print}{mcf}(x, ...) 

\method{plot}{mcf}(x
	, theme.params = getCurrentTheme()
	, xtitle = "Lag"
	, ytitle = expression(rho)
	, overrides = list(...)
	, ...)
}

\arguments{ 
\item{x}{ Instance of class 'mcf'}
\item{theme.params}{ Theme parameters (DEFAULT: getCurrentTheme())}
\item{xtitle}{ Title for the x-axis (DEFAULT: "Lag")}
\item{ytitle}{ Title for the y-axis (DEFAULT: expression(rho))}
\item{overrides}{ List of parameters to override the theme. Must match by name the parameters defined by the theme (DEFAULT: list(...))}
\item{...}{Alternative way to quickly override theme parameters}
}

\value{Void}

\examples{
\dontrun{
# Dow Jones
DJ = get.fs("^DJI", SName = "DowJones", from=as.Date("2008-06-01"), to=as.Date("2009-04-01"));
# Compute Returns
RDJ = Ret(DJ, na.rm = TRUE)

# Compute Multi Correlation Function
res = mcf(RDJ, lag.max = 30, plot = FALSE)
# Plot Autocorrelation Function and Partial ACF 
plot(res)
# Using another theme
plot(res, theme = getTheme("vanilla")) 
}
}

\seealso{\code{\link{cplot}} }
\author{ RAdamant Development Team \email{team@r-adamant.org} }
