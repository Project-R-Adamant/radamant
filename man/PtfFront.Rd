\name{PtfFront}
\alias{PtfFront}

\title{Portfolio efficient frontier}

\description{Compute / Simulate portfolio mean-variance efficient frontier}

\usage{PtfFront(PTF = NULL, ret = NULL, SIGMA = NULL, mi = NULL, n_sim = 10,
volatility = TRUE, plot = TRUE, main = paste("Frontier Simulation:",
ifelse(is.null(mi), n_sim, length(mi)), "points"), xtitle = ifelse(volatility,
expression(sigma), expression(sigma^2)), ytitle = expression(mu), xlab.srt =
0, ytitle.srt = 0, type = "o", legend = "Mean-Variance Frontier", ...) }

\arguments{ \item{PTF}{ }\item{ret}{ }\item{SIGMA}{ }\item{mi}{ }\item{n_sim}{ }\item{volatility}{ }\item{plot}{ }\item{main}{ }\item{xtitle}{ }\item{ytitle}{ }\item{xlab.srt}{ }\item{ytitle.srt}{ }\item{type}{ }\item{legend}{ }\item{...}{ }}

\details{}

\value{}

\references{}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\note{TO BE COMPLETED}

\seealso{}

\examples{}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
