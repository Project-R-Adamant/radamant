\name{logger}

\alias{Logger}

\title{Main logging function}

\description{Send the input message to console and log file.}

\usage{
Logger(message = ""
      , from = deparse(sys.call(sys.parent()))
      , level = 1
      , line = NA
      , env = getOption("RAdamant")
      , console = getConsoleLogging(env = env)
      , logfile = getLogFile(env = env)
      )
}

\arguments{ 
\item{message}{Message to be logged.}
\item{from}{The level in the call stack from which the log message was generated.}
\item{level}{The debug level (importance) of the input message (level >= 1).}
\item{line}{The code line number that the message refers to.}
\item{env}{The environment where the logging options are stored.}
\item{console}{Logical. If TRUE, the message is sent to console.}
\item{logfile}{The filename where the log information is saved.}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\note{This is an internal logging function. It is supposed to be called from other functions.}

