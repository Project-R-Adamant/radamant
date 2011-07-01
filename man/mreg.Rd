\name{mreg}

\alias{mreg}
\title{Multiple regressions}

\description{Multiple regressions}

\usage{
mreg(Y, X, xlabels = NULL, tick.step = 1, backtest = 0, 
stress.idx = c(), type = "simple",
model = "lm", ci = 0.95, max.vars = NCOL(X), 
intercept = TRUE, family = gaussian, weights = NULL,
plot = TRUE, scope = NULL, trace = FALSE, ...) 
}

\arguments{ \item{Y}{ Y}
\item{X}{ X}
\item{xlabels}{ xlabels}
\item{tick.step}{ tisck.step}
\item{backtest}{backtest}
\item{stress.idx}{ stress.idx}
\item{type}{ type}
\item{model}{ model}
\item{ci}{ ci}
\item{max.vars}{ max.vars}
\item{intercept}{ intercepts}
\item{family}{ family}
\item{weights}{ weights}
\item{plot}{ LOGICAL. If TRUE plot is returned.}
\item{trace}{ trace}
\item{scope}{ scope}
\item{...}{ Further arguments to or from other methods}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
