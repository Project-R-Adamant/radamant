\name{Cross_ccf}

\alias{cross.ccf}

\title{Cross correlation function}

\description{Compute the cross correlation function for each pairs of variables (Yi Xj)}

\usage{cross.ccf(Y, X, lag.max = 10, ci = 0.95, plot = TRUE, ...) }

\arguments{ 
\item{Y}{ Matrix of data series (one column per variable)}
\item{X}{ Matrix of data series (one column per variable)}
\item{lag.max}{ Max lag to be computed by the cross correlation function (DEFAULT: 10)}
\item{ci}{ Confidence Interval (DEFAULT: 0.95)}
\item{plot}{ LOGICAL. If TRUE, results are plotted.}
\item{...}{ additional parameters accepted by the function plot.cross.ccf.}
}

\details{}

\value{A list of Ny*Nx cross correlation objects of the class "cool.acf"}

\references{}

\author{ RAdamant Development Team \email{team@r-adamant.org} }

\note{}

\seealso{}

\examples{}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
