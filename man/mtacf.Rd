\name{mtacf}

\alias{print.cool.acf}
\alias{plot.cool.acf}

\title{Cool.Acf methods}

\description{Plot and Print methods for class 'cool.acf'}

\usage{
\method{print}{cool.acf}(x, ...)

\method{plot}{cool.acf}(x
	, theme.params = getCurrentTheme()
	, xtitle = "Lag"
	, ytitle = expression(rho)
	, overrides = list(...)
	, ...
	)
}

\arguments{ 
\item{x}{Instance of class 'cool.acf'}
\item{theme.params}{Theme parameters (DEFAULT: getCurrentTheme())}
\item{xtitle}{Title for the x-axis (DEFAULT: "Lag")}
\item{ytitle}{Title for the y-axis (DEFAULT: expression(rho))}
\item{overrides}{List of parameters to override the theme. Must match by name the parameters defined by the theme (DEFAULT: list(...))}
\item{...}{Alternative way to quickly override theme parameters}
}

\value{Void}

\examples{
# Run Multi correlation analysis
X = mcf(rnorm(30), plot = FALSE);
# Extract cool.acf component
Y = X$ACF[[1]]
class(Y)
# Plot Autocorrelation function
plot(Y)

}
\author{ RAdamant Development Team \email{team@r-adamant.org} }


