\name{preder}

\alias{pred_error}
\alias{av_er}
\alias{abs_avdi}
\alias{mse}
\alias{sde}
\alias{track_sign}
\alias{track_sign_exp}

\title{Prediction error}

\description{Measures for model evaluation}

\usage{
pred_error(target, pred, pc = FALSE) 

av_er(target, pred, pc=FALSE)

abs_avdi(target, pred, pc=FALSE)

mse(target, pred)

sde(target, pred)

track_sign(target, pred)

track_sign_exp(target, pred)
}

\arguments{
\item{target}{Vector. Observed target value}
\item{pred}{Vector. Predicted values}
\item{pc}{ Logical. If TRUE return results in percentage}
}

\details{
\itemize{
\item{pred_error}{: Prediction error}
\item{av_er}{: Average error}
\item{abs_avdi}{: Absolute average discard}
\item{mse}{: Mean squared error}
\item{sde}{: Error standard deviation}
\item{track_sign}{: Error track signal}
\item{track_sign_exp}{: Exponential track signal}
}
}


\author{RAdamant Development Team \email{team@r-adamant.org}}




