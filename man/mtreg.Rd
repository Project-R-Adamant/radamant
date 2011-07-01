\name{MethodReg}

\alias{summary.reg}
\alias{print.reg}
\alias{plot.reg}

\title{Methods for reg}

\description{Plot, Print ND Summary method for "reg"}

\usage{
\method{print}{reg}(x, ...) 
\method{summary}{reg}(object, ...) 
\method{plot}{reg}(x, mode = c("response", "link"), 
title = ifelse(x$model.type == "lm", "LS Regression","GLM Regression"), 
theme.params = getCurrentTheme(), 
overrides = list(...), ...) 
}


\arguments{ \item{x, object}{ x}\item{mode}{ mode}\item{title}{ title}
\item{theme.params}{ theme.params}\item{overrides}{ overrides}\item{...}{Further arguments to or
from other methods}}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\note{TO BE COMPLETED}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
