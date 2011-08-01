\name{plotret}
\alias{plot.ret}
\title{Plot Returns}
\description{Plot method for class "ret"}

\usage{
\method{plot}{ret}(x
    , style = c("line"
    , "bar")
    , xlabels = rownames(x)
    , theme.params = getCurrentTheme()
    , ...
    )
}

\arguments{ 
\item{x}{ an objekt of class "ret"}
\item{style}{ plot style, "line" plot or "bar" plot}
\item{xlabels}{ xlabels}
\item{theme.params}{theme.params}
\item{...}{ Further arguments to or from other methods}}

\value{Void}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\seealso{ \code{\link{Ret}} }

\examples{
# load an example dataset containing financial daily prices
data(ex_fs)
x = ex_fs[ ,1:4]

# calculation and plot for single series
Ret(x[,1]
    , lag = 5
    , plot = TRUE
    , mode = "selected"
    , style="bar"
    , main="Returns - 5 Lags"
    )
    
# calculation and plot for multiple series
par(mfrow = c(2,2))
Ret(x
    , lag = 5
    , mode = "selected"
    , plot = TRUE
    , style = "bar"
    , main = "Returns - 5 Lags"
    )
}



