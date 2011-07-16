\name{runlog}
\alias{write.log}
\alias{error.handling}
\title{Error Handling and Log with runner}
\description{
\code{write.log}: Simple function to write/append log to file (csv format).\cr
\code{error.handling}: Error handling function
}
\usage{
write.log(log = matrix(NA, nrow = 0, ncol = 0), logfile = "runlog.log") 
error.handling(err)
}
\arguments{ 
\item{log}{ Matrix containing logging information.}
\item{logfile}{ Filename of the log}
\item{err}{ List containing the status code of the error.}
}
\details{
Function \code{error.handling} is to be called ONLY inside a tryCatch statement. \cr
It assigns three variables:\cr
	- log.status = "Failed": the status of the execution is set to "Failed" \cr
	- log.message: The error message generated inside the tryCatch \cr
	- res = NA: the result is set to NA
}
\value{VOID}
\author{RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{run}}, \code{\link{multirun}}}


