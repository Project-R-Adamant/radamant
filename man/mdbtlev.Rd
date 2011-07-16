\name{mdbtlev}
\alias{getDebugTraceLevel}
\alias{setDebugTraceLevel}
\title{Manage Debug Trace Level}
\description{
Set and retrieve the level of function nesting for which logging is performed. 
Controls how much information is sent to the log about the execution of each function executed inside the call stack.
}

\usage{
setDebugTraceLevel(level = 1, env = getOption("RAdamant"))
getDebugTraceLevel(env = getOption("RAdamant"))
}

\arguments{
\item{level}{The level of nesting (level >= 1). See details.}
\item{env}{The environment where the info is stored (DEFAULT = getOption("RAdamant")).}
}

\value{
The current value of debug trace level.
}

\details{
The amount of information sent to log depends on the debug trace level:
\itemize{
\item{level = 1: Only top level function calls are logged.}
\item{level = 2: Top and second level function calls (function within a function) are logged.}
\item{level = N: All functions in the call stack up to level N are logged.}
}
}

\examples{
# Retrieve current debug level
getDebugTraceLevel();

# Enable logging to console
setConsoleLogging(TRUE);

# Set minimal level of trace debugging
setDebugTraceLevel(1);
cplot(1:10);

# Set high level of trace debugging (up the 10th level of inner function call)
setDebugTraceLevel(5);
cplot(1:10);

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
