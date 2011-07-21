#### DRAWDOWNS / DRAWUP  ANALYSIS ###
drawdown = function(x, ...) UseMethod("drawdown")
drawdown.default = function(x, FUN=max, relative=FALSE, plots=c("regular", "smooth", "no.plot"), ...){
	# check extreme function
	if(!is.function(FUN)){
		cat("'~' FUNction must be either 'max' or 'min' \n")
		return(NULL)
	}	
	# convert x to matrix
	if(!(is.matrix(x)))
		x = as.matrix(x)
	# check for NAs
	if(any(is.na(x)))
		x = x[-is.na(x)]
	if(is.null(colnames(x)))
		colnames(x) = deparse(substitute(x))
	# series length
	lx = length(x)
	# cumulative returns
	cx = cummax(x)
	# initialise vector of results
	res = rep(0, lx)
	## calculate drawdown
	if(relative){
		res = (cx - x) / max(cx)
	} else {
		# normal drawdown
		res = (cx - x)
	}
	# list of results
	Results = list(Series_info = cbind(Mi=mean(x, na.rm=TRUE), Sigma=sd(x, na.rm=TRUE), Time = lx)
			, Drawdown = res );
	# select plotting options
	pt = match.arg(plots)
	switch(pt,
		"regular" =	cplot(res, ...)
		,
		"smooth" = cplot(sma(res, ...), ...)
		)	
	# assign class	
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
	if(object[[1]][,1] < 0){
		edd = (object[[1]][,1] * object[[1]][,3]) + (object[[1]][,2]^2 / object[[1]][,1])
	} else if (object[[1]][,1] > 0){
		edd = (object[[1]][,2]^2 / object[[1]][,1]) * (0.63519 + 0.5*log(object[[1]][,3]) + log(object[[1]][,1] / object[[1]][,2]))
	} else {
		edd = sqrt(0.5*pi)*object[[1]][,2]*sqrt(object[[1]][,3])
	}
	# expected max drawdown per unit of sigma
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
	ext = FUN(DD[[2]], na.rm=TRUE);
	# perform rolling extreme (max or min)
	if(rolling){
		sext = movApply(DD[[2]], win.size=lag, func=FUN, ...)
	} else {
		sext = NULL
	}
	# list of results
	res = list(ext, sext);
	names(res) = c(deparse(substitute(FUN)), "Rolling")
	if(rolling & plot)
		cplot(sext)
	res;
}