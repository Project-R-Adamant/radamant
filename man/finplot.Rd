\name{finplot}
\alias{fin.plot}
\title{Plot financial time series}
\description{Plot financial time series}

\usage{fin.plot(X, top.vars = c("Close", "High", "Low"), 
bottom.vars = "Volume", style = c("default", "candlestick"), 
snames = attr(X, "SName"), xlabels = rownames(X), 
main = "", main2 = "",
ytitle = "", ytitle2 = "", 
theme.top = getCurrentTheme(), 
overrides = list(...), 
theme.bottom = getCurrentTheme(), 
overrides2 = NULL, ...) }

\arguments{\item{X}{ X}\item{top.vars}{ top.vars}\item{bottom.vars}{ bottom.vars}\item{style}{
style}\item{snames}{ snames}\item{xlabels}{ xlabels}\item{main}{ main}\item{main2}{
main2}\item{ytitle}{ ytitle}\item{ytitle2}{ ytitle2}\item{theme.top}{ theme.top}\item{overrides}{
overrides}\item{theme.bottom}{ theme.bottom}\item{overrides2}{ overrides2}\item{...}{ Further
arguments to or from other methods.}}
\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
