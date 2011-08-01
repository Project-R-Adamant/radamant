\name{fmlmreg}

\alias{formula.mreg}
\alias{formula.reg}

\title{Extract formula from regression object}

\description{Extract formula from regression ("reg" / "mreg") object}

\usage{
\method{formula}{reg}(x, ...) 

\method{formula}{mreg}(x, ...) 
}

\arguments{ 
\item{x}{An object of class "reg" / "mreg"}
\item{...}{Further arguments passed to or from other methods.}
}

\value{
A formula if input x is an object of class "reg". \cr
A list of formulas if x is an object of class "mreg".
}

\seealso{\code{\link{mreg}}}

\author{RAdamant Development Team \email{team@r-adamant.org}}



