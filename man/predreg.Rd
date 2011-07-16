\name{predreg}

\alias{predict.reg}

\title{Predict method for regression}

\description{Predict method for class 'reg'}

\usage{

\method{predict}{reg}(object
	, newdata = NULL
	, ci = 0.95
	, mode = c("response", "link")
	, plot = FALSE
	, shaded = FALSE
	, xlabels = NULL
	, main = "Linear Model Prediction"
	, col = getThemeAttr("col", exact = TRUE)[c(1, 2, 2)]
	, shade.stripes = 1
	, shade.col = getThemeAttr("col", exact = TRUE)[2]
	, shade.density = 40
	, shade.angle = 30
	, legend = NULL
	, ...
	) 
}

\arguments{ \item{object}{ OBJECT of class "reg".}
\item{newdata}{ newdata}
\item{ci}{ ci}
\item{mode}{mode}
\item{plot}{ LOGICAL. If TRUE plot is returned.}
\item{shaded}{ shaded}
\item{xlabels}{ xlabels}
\item{main}{ main}
\item{col}{ color}
\item{shade.stripes}{ shade.stripes}
\item{shade.col}{ shade.col}
\item{shade.density}{ shade.density}
\item{shade.angle}{ shade.angle}
\item{legend}{ legend}
\item{...}{Further arguments to or from other methods}}

\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}


