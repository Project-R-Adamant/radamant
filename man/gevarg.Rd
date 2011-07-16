\name{gevarg}
\alias{gev.VaR.range}
\title{GEV - VaR range grid for contour calculation}
\description{GEV - VaR range grid for contour calculation}

\usage{gev.VaR.range(Xbmax, 
VaR = sum(gev.VaR.constraint(parms = c(0, xi, sigma), 
type = "both", Xbmax = Xbmax, prob = prob))/2, xi = 0.1, 
sigma = 1, alpha = 0.01, df = 3, prob = alpha[1], ...) }
\arguments{\item{Xbmax}{ Xbmax}\item{VaR}{ VaR}\item{xi}{ xi}\item{sigma}{ sigma}\item{alpha}{
alpha}\item{df}{ df}\item{prob}{ prob}\item{...}{ Further arguments to or from other methods.}}
\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}


