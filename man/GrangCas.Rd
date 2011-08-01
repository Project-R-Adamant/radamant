\name{grangcas}

\alias{GrangCas.VecAr}
\alias{print.GrangCas}

\title{Granger Causality test}

\description{Perform Granger causality test for parameters of VAR model}

\usage{
\method{GrangCas}{VecAr}(X, cause = colnames(coef(X)), digits = 3, ...) 
}

\arguments{ 
\item{X}{ An object of class "VecAr"}
\item{cause}{ Vector of character. Name of the variables to be used as "cause". By default all the variables are tested.}
\item{digits}{ number of digits to be printed.}
\item{...}{ Further arguments to or from other methods.}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}


