\name{stepmat}

\alias{StepMat}

\title{Step matrix for binomial tree}

\description{Simulate binomial path of a binomial tree}

\usage{StepMat(init, n_step, up, down) }

\arguments{
\item{init}{ Initial price, step number 0 in the matrix.}
\item{n_step}{ Integer. Number of steps.}
\item{up}{ Up movement factor}
\item{down}{ Down movement factor}
}

\value{
Create Step probability matrix of \emph{(n_step+1) x (n_stpe+1)} dimensions
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
# simulate binomial path for 10 steps
StepMat(init = 0.5, n_step = 10, up = 0.8, down = 0.6)
}

