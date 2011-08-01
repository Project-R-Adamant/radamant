\name{strvar}

\alias{Strvar.VecAr}

\title{Structural Vector Autoregressive model}

\description{Estimate Structural Vector Autoregressive model}

\usage{
Strvar.VecAr(X, A = "diag", B = NULL, inter = FALSE, ...) 
}

\arguments{ 
\item{X}{ An object of class "VecAr"}
\item{A}{ Restriction matrix A.}
\item{B}{ Restriction matrix B.}
\item{inter}{ Logical. If TRUE restrictions matrix will be manually edited.}
\item{...}{ Further arguments to or from other methods}
}

\value{An object list containing the following elements: 
\item{EST_Matrix}{ List of 2 elements: 
    \itemize{
        \item{Estimated A parameters}
        \item{Estimated B parameters}
        }
    }
\item{SE}{ List of 2 elements: 
    \itemize{
        \item{Standard errors of A parameters}
        \item{Standard errors of B parameters}
        }
    }
\item{LogLik}{Log-Likelihood value.}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\seealso{ \code{\link{optim}} ,\code{\link{VecAr}}}

\examples{
# load example data sete
data(ex_ptf)
X = ex_ptf[ ,1:4]
# estimate VAR(2) model
vecar = VecAr(X, ar.lags=1:2, type="const")

## Estimate Structural VAR models
# EX. 1
# Default constraints provided by the function:
# A =      [,1] [,2] [,3] [,4]
#        [1,]    C1    0    0    0
#        [2,]    0    C2    0    0
#        [3,]    0    0    C3    0
#        [4,]    0    0    0    C4
# B =      [,1] [,2] [,3] [,4]
#        [1,]    1    0    0    0
#        [2,]    0    1    0    0
#        [3,]    0    0    1    0
#        [4,]    0    0    0    1

Strvar.VecAr(vecar)

# EX. 2
# Different constraints for A matrix:
# A =      [,1] [,2] [,3] [,4]
#        [1,]    C1    0    0    0
#        [2,]    C2    C3    0    0
#        [3,]    C4    0    C5    0
#        [4,]    C6    0    0    C6
# B =      [,1] [,2] [,3] [,4]
#        [1,]    1    0    0    0
#        [2,]    0    1    0    0
#        [3,]    0    0    1    0
#        [4,]    0    0    0    1

A = diag(NA, 4)
A[ ,1] = NA
Strvar.VecAr(vecar, A=A)
}


