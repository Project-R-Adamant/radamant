\name{mtccf}

\alias{print.cross.ccf}
\alias{plot.cross.ccf}

\title{Cross.ccf functions}

\description{Methods for class 'cross.ccf'}

\usage{
\method{print}{cross.ccf}(x, ...)

\method{plot}{cross.ccf}(x
    , theme.params = getCurrentTheme()
    , xtitle = "Lag"
    , ytitle = expression(rho)
    , overrides = list(...)
    , ...
    )
}

\arguments{ 
\item{x}{Instance of class 'cross.ccf'}
\item{theme.params}{Theme parameters (DEFAULT: getCurrentTheme())}
\item{xtitle}{Title for the x-axis (DEFAULT: "Lag")}
\item{ytitle}{Title for the y-axis (DEFAULT: expression(rho))}
\item{overrides}{List of parameters to override the theme. Must match by name the parameters defined by the theme (DEFAULT: list(...))}
\item{...}{Alternative way to quickly override theme parameters.}
}

\value{Void}

\examples{
# Generate two random integrated series
N = 100
X = matrix(rnorm(N), nrow = N/2, ncol=2);
# Create two series as a linear combination of X plus noise
Y = X %*% diag(c(1, -1)) + 0.1*rnorm(N);
# Perform Cross Correlation Analysis
Z = cross.ccf(Y, X, plot = FALSE)
plot(Z)
}

\author{ RAdamant Development Team \email{team@r-adamant.org} }


