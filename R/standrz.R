#### STANDARDIZATION ####
# Decimal scale
	# x series
	# scale 
Decscal=function(x,scale=0.1){
	# numeric input
	Logger(message = "numeric input", from = "Decscal", line = 2, level = 1);
	if(!is.matrix(x))
		x=as.matrix(x)
	# NAs 
	Logger(message = "NAs ", from = "Decscal", line = 5, level = 1);
	if(any(is.na(x)))
		naid = !is.na(x)
	# scale series
	Logger(message = "scale series", from = "Decscal", line = 8, level = 1);
	res=x[naid]/(10)^scale
	# clean memory
	Logger(message = "clean memory", from = "Decscal", line = 10, level = 1);
	cleanup(res)
	# results
	Logger(message = "results", from = "Decscal", line = 12, level = 1);
	res
}
# Min-Max scale
	# x
	# tmin
	# tmax
Minmaxscal=function(x,tmin=0,tmax=1){
	M = FALSE
	# numeric input
	Logger(message = "numeric input", from = "Minmaxscal", line = 3, level = 1);
	if(length(dim(x)) > 1L){
		M = TRUE
		x=as.matrix(x)
		}
	# min value
	Logger(message = "min value", from = "Minmaxscal", line = 8, level = 1);
	mn = if(M) apply(x,2,min,na.rm=TRUE) else
		min(x,na.rm=TRUE)
	# max value
	Logger(message = "max value", from = "Minmaxscal", line = 11, level = 1);
	mx = if(M) apply(x,2,max,na.rm=TRUE) else
		max(x,na.rm=TRUE)
	# rescale series
	Logger(message = "rescale series", from = "Minmaxscal", line = 14, level = 1);
	if(M){
		res = matrix(0,nrow(x),ncol(x))
			for(i in 1:ncol(x)){
				res[,i] = ((x[,i]- mn[i])/(mx[i]-mn[i]))*(tmax-tmin)+tmin
				}
		} else {
		res = ((x-mn)/(mx-mn))*(tmax-tmin)+tmin
		}
	# clean memory
	Logger(message = "clean memory", from = "Minmaxscal", line = 23, level = 1);
	#cleanup(res)
	Logger(message = "cleanup(res)", from = "Minmaxscal", line = 24, level = 1);
	# results
	Logger(message = "results", from = "Minmaxscal", line = 25, level = 1);
	res
}
# Z-Index
	# x
	# sigma
	# mu
Zind=function(x,sigma=1,mi=2){
	M = FALSE
	# numeric input
	Logger(message = "numeric input", from = "Zind", line = 3, level = 1);
	if(length(dim(x)) > 1L){
		M = TRUE
		x=as.matrix(x)
		}
	# get sigma
	Logger(message = "get sigma", from = "Zind", line = 8, level = 1);
	if (missing(sigma)){
		if(M){
		sigma = apply(x,2,sd,na.rm=TRUE)
		} else {
		sigma = sd(x,na.rm=TRUE)
		}
	}
	# get mu
	Logger(message = "get mu", from = "Zind", line = 16, level = 1);
	if (missing(mi)){
		if(M){
		mi = apply(x,2,mean,na.rm=TRUE)
		} else {
		mi = mean(x,na.rm=TRUE)
		}
	}
	# return originale series
	Logger(message = "return originale series", from = "Zind", line = 24, level = 1);
	if (mi == 0 & sigma ==1) 
		return(res)
	# calc z-index
	Logger(message = "calc z-index", from = "Zind", line = 27, level = 1);
	if(M){
		res = matrix(0,nrow(x),ncol(x))
			for(i in 1:ncol(x)){
				res[,i] = (x[,i]- mi[i])/sigma[i]
				}
	} else {
		res = (x - mi) / sigma
	}
	# clean memory
	Logger(message = "clean memory", from = "Zind", line = 36, level = 1);
	#cleanup(res)
	Logger(message = "cleanup(res)", from = "Zind", line = 37, level = 1);
	# results
	Logger(message = "results", from = "Zind", line = 38, level = 1);
	res
}
