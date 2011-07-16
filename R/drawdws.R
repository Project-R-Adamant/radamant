#### DRAWDOWNS / DRAWUP  ANALYSIS ###
drawdown = function(x, ...) UseMethod("drawdown")
drawdown.default = function(x, FUN=max, relative=FALSE, plot=FALSE, ...){
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
	# series length
	Logger(message = "series length", from = "drawdown.default", line = 13, level = 1);
	lx = length(x)
	# cumulative returns
	Logger(message = "cumulative returns", from = "drawdown.default", line = 15, level = 1);
	cx = cumsum(x)
	# initialise vector of results
	Logger(message = "initialise vector of results", from = "drawdown.default", line = 17, level = 1);
	i = 1
	res = rep(0, lx)
	## calculate drawdown
	Logger(message = "calculate drawdown", from = "drawdown.default", line = 20, level = 1);
	if(relative){
		# relative drawdown
		Logger(message = "relative drawdown", from = "drawdown.default", line = 22, level = 1);
		while(i <= lx){
			res[i] = FUN(cx[i] - x[i]) / FUN(cx[i])
			i = i + 1
		}
	} else {
		# normal drawdown
		Logger(message = "normal drawdown", from = "drawdown.default", line = 28, level = 1);
		while(i <= lx){
			res[i] = FUN(cx[i] - x[i]) 
			i = i + 1
		}
	}
	if(plot){
		cplot(res, ...)
	}
	
	Results = list(Series_info = cbind(Mi=mean(x, na.rm=TRUE), Sigma=sd(x, na.rm=TRUE), Time = lx) , Drawdown = res );
	
	# assign class	
	Logger(message = "assign class	", from = "drawdown.default", line = 38, level = 1);
	class(Results) = "drawdown";
	Results;
}
## Summary drawdown
SummaryDD = function(DD){
	
	if(class(DD) != "drawdown"){
		cat("DD must be an objeckt of class 'drawdown' \n")
		return(NULL)
	}
	sums = c(mean(DD[[2]], na.rm=TRUE),sd(DD[[2]], na.rm=TRUE),max(DD[[2]], na.rm=TRUE),min(DD[[2]], na.rm=TRUE)
			);
	
	# asyntotic expected MDD
	Logger(message = "asyntotic expected MDD", from = "SummaryDD", line = 8, level = 1);
	if(DD[[1]][,1] < 0){
		
		edd = (DD[[1]][,1] * DD[[1]][,3]) + (DD[[1]][,2]^2 / DD[[1]][,1])
		
	} else if (DD[[1]][,1] > 0){
		
		edd = (DD[[1]][,2]^2 / DD[[1]][,1]) * (0.63519 + 0.5*log(DD[[1]][,3]) + log(DD[[1]][,1] / DD[[1]][,2]))
	} else {
		
		edd = sqrt(0.5*pi)*DD[[1]][,2]*sqrt(DD[[1]][,3])
	}
	
	# expected max drawdown per unit of sigma
	Logger(message = "expected max drawdown per unit of sigma", from = "SummaryDD", line = 16, level = 1);
	edds = edd / DD[[1]][,2]
	
	Results = list(summary=sums, Expected_DD=edd, Expected_DD_sigma=edds);
	
	cat("!Expected DD is calculated with asyntotic formulas! \n")		
	Results;
}
## max / min drawdown
ExtremeDD = function(DD, FUN, lag=1, rolling=FALSE, plot=TRUE, ...){
	
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
