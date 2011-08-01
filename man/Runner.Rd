\name{runner}
\alias{run}
\alias{multirun}

\title{Runner and Multirunner}
\description{Wrapper function to execute any function. \cr
Run single or multiple functions and provide a list of results.
}

\usage{
run(func = NULL
    , args = list()
    , writelog = TRUE
    , logfile = "runlog.log"
    , check.input = TRUE
    , output = c("console", "sing.file")
    ) 

multirun(func.array = character(0)
    , args.list = list()
    , writelog = TRUE
    , logfile = "runlog.log"
    , output = c("console", "sing.file", "multi.file")
    ) 
}

\arguments{ 
\item{func}{ Name of the function to run}
\item{func.array}{ Array of function names to execute}
\item{args}{ Named list of parameters of the function.\cr
    Each entry is of the form: args[["PARAM.NAME"]] = VALUE.
}
\item{args.list}{ Array of named list of parameters of the function.\cr
    Each entry is a list of parameters, as required by the wrapper function
"run".
}
\item{writelog}{Logical. If TRUE, execution log is written to file.}
\item{logfile}{ Filename of the log}
\item{check.input}{Logical. If TRUE, basic checks are performed on input
data, and stop code execution in case of wrong data.
}
\item{output}{ Choose wether to return the results in the console or export the to text file.}
}
\details{ When called the function \code{multirun} the elements of the argument \code{args.list} can
be 
specified with or without names. If the names are specified the arguments can be put in a different
order from
the array function. \cr
If \code{writelog = TRUE} a log containing information about submitted
computation is saved in the current working directory.
If \code{output} = "sing.file", a text file containing all the results is saved in current working
directory. \cr
    The file will be named "Run_time_date.txt"
If \code{output} = "sing.file", a text for each called function is saved in a text file. \cr
    The files will be named "Function Name_time_date.txt"
}
\value{The object returned depends on the function being called. \cr
\code{multirun} returns a list of results, one entry for each function being
executed.}
\author{RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{write.log}}, \code{\link{error.handling}}}
\examples{
# Run Exponential Moving Average and Simple Moving Average.
# For each function a list of parameters has been specified
multirun(c("ema","sma")
	,list(list(rnorm(150), 5)
		  , list(rnorm(100), 10)
		  )
	, writelog = TRUE
	)
# Specifies names in the list of arguments
multirun(func.array = c("ema","sma")
	,args.list = list(sma = list(rnorm(150), 5)
					 , ema = list(rnorm(100), 30)
					 )
	, TRUE
	)
# Output to text file
multirun(func.array = c("ema","sma")
	, args.list = list(sma = list(rnorm(150), 5)
					  , ema = list(rnorm(100), 30)
					  )
	, output = "multi.file"
	)
}


