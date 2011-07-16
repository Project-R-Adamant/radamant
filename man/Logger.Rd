\name{logger}

\alias{Logger}

\title{Main logging function}

\description{Create Log for the functions contained in the package}

\usage{Logger(message = "", from = deparse(sys.call(sys.parent())), 
level = 1, line = NA, env = getOption("RAdamant"), 
console = getConsoleLogging(env = env), 
logfile = getLogFile(env = env)) }

\arguments{ 
\item{message}{ Message printed.}
\item{from}{ from}
\item{level}{ Log depth level, minimum = 1}
\item{line}{ line}
\item{env}{ environment}
\item{console}{ console logging}
\item{logfile}{ log file}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}


