\name{decimals}
\alias{decimals}

\title{Count Decimals}
\description{Count the number of digits of the decimal part (mantissa) of a number }
\usage{decimals(x, max.digits = 10, ...)}

\arguments{
\item{x}{The number for which the count of decimals is required.}
\item{max.digits}{Controls the resolution. See details.}
\item{...}{Not used, for future releases.}
}

\details{
The number x is first converted into a string, where the decimal part is truncated after max.digits.
The number of significant digits of the decimal part are hence calculated. 
The truncation allows to remove the artifacts introduced by the finite resolution of the numbers representation. 
}

\value{The number of digits of the mantissa}

\examples{

# Print the mantissa of the number 0.01 with 80 digits.
sprintf("\%.80f", 0.01)

# Real number of digits is 2
decimals(0.01, max.digits = 10)

# Number of digits of the mantissa of the computer representation of 0.01
decimals(0.01, max.digits = 100)

}

\author{RAdamant Development Team \email{team@r-adamant.org}}

