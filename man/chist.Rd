\name{chist}
\alias{chist}
\title{Custom histogram function}
\description{Custom histogram function}

\usage{
chist(x
	, nclass = min(max(round(NROW(x)/10), 10), NROW(x))
	, density = c("kernel", "normal")
	, kernel = c("gaussian", "epanechnikov", "rectangular"
				, "triangular", "biweight", "cosine", "optcosine")
	, theme.params = getCurrentTheme()
	, main = "Histogram and Kernel Density Estimation"
	, xtitle = NULL
	, ytitle = NULL
	, legend = NULL
	, show.legend = TRUE
	, normalised = FALSE
	, ...
	) 
}

\arguments{
\item{x}{ x}
\item{nclass}{ nclass}
\item{density}{ density}
\item{kernel}{kernel}
\item{theme.params}{ theme.params}
\item{main}{ main}
\item{xtitle}{ xtitle}
\item{ytitle}{ytitle}
\item{legend}{ legend}
\item{show.legend}{ show.legend}
\item{normalised}{normalised}
\item{...}{ ...}}
\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}


