#######################################################################################################################
# Copyright (C) 2011  RAdmant Development Team
# email: team@r-adamant.org
# web: http://www.r-adamant.org
#
# This library is free software;
# you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation;
# either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program;
# if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA.
#######################################################################################################################
#### DRAWDOWNS / DRAWUP  ANALYSIS ###
drawdown = function(x, ...) UseMethod("drawdown")
drawdown.default = function(x, FUN=max, relative=FALSE, plots=c("regular", "smooth", "no.plot"), ...){
	# check extreme function
	Logger(message = "check extreme function", from = "drawdown.default", line = 2, level = 1);
	if(!is.function(FUN)){
		cat("'~' FUNction must be either 'max' or 'min' \n")
		return(NULL)
	}	
	# convert x to matrix
	Logger(message = "convert x to matrix", from = "drawdown.default", line = 7, level = 1);
	if(!(is.matrix(x)))
		x = as.matrix(x)
	# check for NAs
	Logger(message = "check for NAs", from = "drawdown.default", line = 10, level = 1);
	if(any(is.na(x)))
		x = x[-is.na(x)]
	if(is.null(colnames(x)))
		colnames(x) = deparse(substitute(x))
	# series length
	Logger(message = "series length", from = "drawdown.default", line = 15, level = 1);
	lx = length(x)
	# cumulative returns
	Logger(message = "cumulative returns", from = "drawdown.default", line = 17, level = 1);
	cx = cummax(x)
	# initialise vector of results
	Logger(message = "initialise vector of results", from = "drawdown.default", line = 19, level = 1);
	i = 1
	res = rep(0, lx)
	## calculate drawdown
	Logger(message = "calculate drawdown", from = "drawdown.default", line = 22, level = 1);
	if(relative){
		res = (cx - x) / max(cx)
	} else {
		# normal drawdown
		Logger(message = "normal drawdown", from = "drawdown.default", line = 26, level = 1);
		res = (cx - x)
	}
	# list of results
	Logger(message = "list of results", from = "drawdown.default", line = 29, level = 1);
	Results = list(Series_info = cbind(Mi=mean(x, na.rm=TRUE), Sigma=sd(x, na.rm=TRUE), Time = lx)
			, Drawdown = res );
	# select plotting options
	Logger(message = "select plotting options", from = "drawdown.default", line = 32, level = 1);
	pt = match.arg(plots)
	switch(pt,
		"regular" =	cplot(res, ...)
		,
		"smooth" = cplot(sma(res, ...), ...)
		)	
	# assign class	
	Logger(message = "assign class	", from = "drawdown.default", line = 39, level = 1);
	class(Results) = "drawdown";
	Results;
}
## Summary drawdown
summary.drawdown = function(object, show.extr = TRUE, ...){
	if(class(object) != "drawdown"){
		cat("DD must be an objeckt of class 'drawdown' \n")
		return(NULL)
	}
	sums = c(Mean = mean(object[[2]], na.rm=TRUE)
			, Std.Dev = sd(object[[2]], na.rm=TRUE)
			, Max = max(object[[2]], na.rm=TRUE)
			, Min = min(object[[2]], na.rm=TRUE)
			);
	# asyntotic expected MDD
	Logger(message = "asyntotic expected MDD", from = "summary.drawdown", line = 11, level = 1);
	if(object[[1]][,1] < 0){
		edd = (object[[1]][,1] * object[[1]][,3]) + (object[[1]][,2]^2 / object[[1]][,1])
	} else if (object[[1]][,1] > 0){
		edd = (object[[1]][,2]^2 / object[[1]][,1]) * (0.63519 + 0.5*log(object[[1]][,3]) + log(object[[1]][,1] / object[[1]][,2]))
	} else {
		edd = sqrt(0.5*pi)*object[[1]][,2]*sqrt(object[[1]][,3])
	}
	# expected max drawdown per unit of sigma
	Logger(message = "expected max drawdown per unit of sigma", from = "summary.drawdown", line = 19, level = 1);
	edds = edd / object[[1]][,2]
	Results = vector("list", 2)
	Results[[1]] = list(summary=sums, Expected_DD=edd, Expected_DD_sigma=edds);
	if(show.extr)
		Results[[2]] = ExtremeDD(DD = object, ...)
	else 
		Results[[2]] = NULL
	cat("!Expected DD is calculated with asyntotic formulas! \n")		
	Results;
}
## max / min drawdown
ExtremeDD = function(DD, FUN=max, lag=1, rolling=FALSE, plot=TRUE, ...){
	if(class(DD) != "drawdown"){
		cat("DD must be an objeckt of class 'drawdown' \n")
		return(NULL)
	}
	# calculate max or min
	Logger(message = "calculate max or min", from = "ExtremeDD", line = 6, level = 1);
	ext = FUN(DD[[2]], na.rm=TRUE);
	# perform rolling extreme (max or min)
	Logger(message = "perform rolling extreme (max or min)", from = "ExtremeDD", line = 8, level = 1);
	if(rolling){
		sext = movApply(DD[[2]], win.size=lag, func=FUN, ...)
	} else {
		sext = NULL
	}
	# list of results
	Logger(message = "list of results", from = "ExtremeDD", line = 14, level = 1);
	res = list(ext, sext);
	names(res) = c(deparse(substitute(FUN)), "Rolling")
	if(rolling & plot)
		cplot(sext)
	res;
}
