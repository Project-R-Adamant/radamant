\name{funcomx}
\alias{func.comment.idx}
\title{Function comment}

\description{Given an input file, this functions created an index based commented version of the
file.}

\usage{
func.comment.idx(control.df = 
	data.frame(FNAME = c(), FCODE = c(), 
	AREA = c(), SECTION = c(), CLASS = c()), 
infile = NULL, incode = NULL, outfile = NULL, max.dgt = 3) }

\arguments{
\item{control.df}{ List of function names. See Details}
\item{infile}{ Input file (Full path: Mandatory).}
\item{incode}{ Input code array (Alternative to infile: Mandatory). Each entry is considered to be a
line of code.}
\item{outfile}{ Output commented file (Full path: Optional). If provided, an output file is
generated.}
\item{max.dgt}{ Controls the number of digits to be used on each section of the comment.}
}

\details{
This data frame is a list of function names:
	\itemize{
		\item{FNAME = Name of the function}
		\item{FCODE = code identifier for the function. (a-Z)(0-9).}
		\item{AREA = Macro area (Description) classification for the function.}
		\item{SECTION = Section (Description) classification for the function (Sub-AREA)}
		\item{CLASS = The class of the returned object.}
	}
}

\value{
String array where every entry is a line of code. Each original line of the input code is preceded
by a special comment.
}

\examples{
#### EXAMPLE #####
tst = data.frame(FNAME = c("sd", "lm")
				, FCODE = c("SD", "LM")
				, AREA = c("s5", "s2")
				, SECTION = c("s1", "s1")
				, CLASS = c("c1", "c2")
				);
				
incode = rbind(paste("sd =", as.character(deparse(args(sd)))[1])
				, as.matrix(deparse(body(sd)))
				, ""
				, ""
				, paste("lm =", as.character(deparse(args(lm)))[1])
				, as.matrix(deparse(body(lm)))
			   )
func.comment.idx(tst, incode = incode, max.dgt=3)
}


\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}


