\name{ptfront}
\alias{PtfFront}
\title{Portfolio efficient frontier}
\description{Compute/Simulate portfolio mean-variance efficient frontier}

\usage{
PtfFront(PTF = NULL
    , ret = NULL
    , SIGMA = NULL
    , mi = NULL
    , n_sim = 10
    , volatility = TRUE
    , plot = TRUE
    , main = paste("Frontier Simulation:"
                   , ifelse(is.null(mi)
                            , n_sim
                            , length(mi)
                            )
                   , "points"
                   )
    , xtitle = ifelse(volatility
                      , expression(sigma)
                      , expression(sigma^2)
                      )
    , ytitle = expression(mu)
    , xlab.srt = 0
    , ytitle.srt = 0
    , type = "o"
    , legend = "Mean-Variance Frontier"
    , ...
    ) 
}

\arguments{
\item{ PTF}{ PTF}
\item{ ret}{ ret}
\item{ SIGMA}{ SIGMA}
\item{ mi}{ mi}
\item{ n_sim}{n_sim}
\item{ volatility}{ volatility}
\item{ plot}{ plot}
\item{ main}{ main}
\item{ xtitle}{xtitle}
\item{ ytitle}{ ytitle}
\item{ xlab.srt}{ xlab.srt}
\item{ ytitle.srt}{ ytitle.srt}
\item{ type}{type}
\item{ legend}{ legend}
\item{ ...}{ Further arguments to or from other methods}}

\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}
