\name{irsvecar}

\alias{IRS.VecAr}
\alias{PHI.VecAr}

\title{VAR Impulse response}

\description{Compute Impulse response function and Wold decomposition for VAR model}

\usage{
IRS.VecAr(X, imp, resp = NULL, steps = 5, cum = TRUE, ortho = FALSE, ...) 
PHI.VecAr(X, steps, ortho = FALSE, ...)
}

\arguments{
    \item{X}{ An object of class "VecAr".}
    \item{imp}{ Vector of characters. Impulse variable(s).}
    \item{resp}{ Vector of characters. Response variable(s).}
    \item{steps}{ Integer. Number of forward steps.}
    \item{cum}{ Logical. If TRUE cumulated impulse will be returned.}
    \item{ortho}{ Logical. If TRUE orthogonal impulse will be returned.}
    \item{...}{ Further arguments to or from other methods.}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
data(ex_ptf)
colnames(ex_ptf)
X = ex_ptf[ ,1:4]
# estimate VAR(2) model
var = VecAr(X, ar.lags=1:2, type="const", exog = NULL)

PHI.VecAr(var, steps=10, ortho=TRUE)

# Impulse response function - single inpulse
imp = "Asset_1"
resp = c("Fund", "Asset_1", "Asset_2", "Asset_3")
im = IRS.VecAr(var, imp=imp, resp=resp, steps=10, ortho=TRUE)
im
# view plots
cplot(im[[1]], lwd=2)

}

