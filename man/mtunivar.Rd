\name{mtunivar}

\alias{plot.univar}
\alias{print.univar}
\alias{summary.univar}

\title{Methods for univariate analysis}

\description{Print, Plot and Summary methods for class 'univar'}

\usage{
\method{summary}{univar}(object, ...) 
\method{plot}{univar}(x, theme.params = getCurrentTheme(), overrides = NULL, ...)
\method{print}{univar}(x, ...) 
}

\arguments{ 
\item{x, object}{ Instance of class 'univar'}
\item{theme.params}{ params: Theme parameters (DEFAULT: getCurrentTheme())}
\item{overrides}{ list of parameters to override the theme. Must match by name the parameters
defined by the theme (DEFAULT: NULL)}
\item{...}{ Further arguments to or fromother methods}
}


\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{ \code{\link{univar}}}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
