\name{crosplot}

\alias{cross.plot}
\title{Cross Plot}

\description{Plot the input dependent variable Y versus each input independent variable X}

\usage{cross.plot(Y, X, theme.params = getCurrentTheme(), 
xlabels = NULL, two.axis = TRUE,
shaded.first = FALSE, overrides = NULL) 
}

\arguments{ 
\item{Y}{ serie of the dependent variable}
\item{X}{ Matrix containing all independent variables (one column per variable)}
\item{theme.params}{ theme parameters (DEFAULT: getCurrentTheme())}
\item{xlabels}{ serie of the lables associated to the rows of X (i.e. Time  libels)(DEFAULT: NULL)}
\item{two.axis}{ LOGICAL. If TRUE, series are plotted on two axis (two scales).}
\item{shaded.first}{ LOGICAL. If TRUE, the variable Y is shaded.}
\item{overrides}{ list of parameters to override the theme. Must match by name the parameters
defined by the theme (DEFAULT: NULL)}
}

\value{VOID}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
