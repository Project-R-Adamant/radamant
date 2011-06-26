#### DRAWDOWNS / DRAWUP  ANALYSIS ###

## DRAWDOWNS ## 
drawdown = function(x, ...) UseMethod("drawdown")

drawdown.default = function(x, FUN=max, relative=FALSE, plot=FALSE, ...)
{
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

	# series length
	lx = length(x)

	# cumulative returns
	cx = cumsum(x)

	# initialise vector of results
	i = 1
	res = rep(0, lx)

	## calculate drawdown
	if(relative){
		# relative drawdown
		while(i <= lx){
			res[i] = FUN(cx[i] - x[i]) / FUN(cx[i])
			i = i + 1
		}
	} else {
		# normal drawdown
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
	if(DD[[1]][,1] < 0){
		
		edd = (DD[[1]][,1] * DD[[1]][,3]) + (DD[[1]][,2]^2 / DD[[1]][,1])
		
	} else if (DD[[1]][,1] > 0){
		
		edd = (DD[[1]][,2]^2 / DD[[1]][,1]) * (0.63519 + 0.5*log(DD[[1]][,3]) + log(DD[[1]][,1] / DD[[1]][,2]))

	} else {
		
		edd = sqrt(0.5*pi)*DD[[1]][,2]*sqrt(DD[[1]][,3])
	}
	
	# expected max drawdown per unit of sigma
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
		cool.plot(sext)

	res;

}


