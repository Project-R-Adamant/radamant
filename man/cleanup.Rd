\name{cleanup}
\alias{cleanup}
\title{Clean memory}
\description{Cleanup environment and (optionally) performs Garbage Collection}

\usage{cleanup(keep = c(), env = parent.frame(), gc = FALSE) }

\arguments{
\item{keep}{ CHARACTER. Vector of variables to keep in memory.}
\item{env}{ Environment from which objects are removed. Defaults to the environment from which this
function is called.}
\item{gc}{ Logical. If TRUE, garbage collection is performed to release memory. (Default = TRUE)}
}

\value{Void}

\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}


