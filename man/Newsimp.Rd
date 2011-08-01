\name{newsimp}

\alias{newsimp.default}
\alias{newsimp.Garch}
\alias{newsimp}

\title{News impact curve }

\description{Compute News impact curve for Garch models}

\usage{
newsimp(x, ...)
\method{newsimp}{Garch}(x, plot = TRUE, ...)
\method{newsimp}{default}(x
	, theta
	, order
	, type=c("garch","mgarch", "egarch","tgarch")
	, plot=FALSE
	, ...) 
}

\arguments{ 
    \item{x}{ A vector of innvations (x axis of the plot) or an object of class "Garch".}
    \item{theta}{ Vector of Garch model parameters.}
    \item{order}{ Vector of integers. Arch and Garch parameters order. (Default: 1,1)}
    \item{type}{ Type of Garch to be estimated: "garch", "mgarch", "tgarch", "egarch". (Default: "garch"). }
    \item{plot}{Logical. If TRUE plot of the NIC is returned.}
    \item{...}{ Further arguments to or from other methods}
}

\value{
The function returns the NIC curve plus a matrix containing: Sigma values (y axis) and Innovations (x axis). \cr
The plot is made by the cplot function, for more information about the graphical parameters take a look here \code{\link{cplot}}.
}

\seealso{ \code{\link{Garch}}, \code{\link{cplot}}}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
    # load example time series 
    data(ex_ts)
    x = ex_ts
    
    # Symmetric NIC - GARCH example 
    gg1 = Garch(x, order = c(2,1), type="garch", prob="g")
    newsimp(gg1)
    
    # Asymmetric NIC - EGARCH and TGARCH example
    gg2 = Garch(x, type="egarch", prob="g")
    newsimp(x=gg2)
    gg3 = Garch(x, type="tgarch")
    newsimp(x=gg3)
}

