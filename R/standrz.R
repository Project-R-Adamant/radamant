#### STANDARDIZATION ####
# Decimal scale
	# x series
	# scale 
Decscal=function(x,scale=0.1){
	# numeric input
	if(!is.matrix(x))
		x=as.matrix(x)
	# NAs 
	if(any(is.na(x)))
		naid = !is.na(x)
	# scale series
	res=x[naid]/(10)^scale
	# clean memory
	cleanup(res)
	# results
	res
}
# Min-Max scale
	# x
	# tmin
	# tmax
Minmaxscal=function(x,tmin=0,tmax=1){
	M = FALSE
	# numeric input
	if(length(dim(x)) > 1L){
		M = TRUE
		x=as.matrix(x)
		}
	# min value
	mn = if(M) apply(x,2,min,na.rm=TRUE) else
		min(x,na.rm=TRUE)
	# max value
	mx = if(M) apply(x,2,max,na.rm=TRUE) else
		max(x,na.rm=TRUE)
	# rescale series
	if(M){
		res = matrix(0,nrow(x),ncol(x))
			for(i in 1:ncol(x)){
				res[,i] = ((x[,i]- mn[i])/(mx[i]-mn[i]))*(tmax-tmin)+tmin
				}
		} else {
		res = ((x-mn)/(mx-mn))*(tmax-tmin)+tmin
		}
	# clean memory
	#cleanup(res)
	# results
	res
}
# Z-Index
	# x
	# sigma
	# mu
Zind=function(x,sigma=1,mi=2){
	M = FALSE
	# numeric input
	if(length(dim(x)) > 1L){
		M = TRUE
		x=as.matrix(x)
		}
	# get sigma
	if (missing(sigma)){
		if(M){
		sigma = apply(x,2,sd,na.rm=TRUE)
		} else {
		sigma = sd(x,na.rm=TRUE)
		}
	}
	# get mu
	if (missing(mi)){
		if(M){
		mi = apply(x,2,mean,na.rm=TRUE)
		} else {
		mi = mean(x,na.rm=TRUE)
		}
	}
	# return originale series
	if (mi == 0 & sigma ==1) 
		return(res)
	# calc z-index
	if(M){
		res = matrix(0,nrow(x),ncol(x))
			for(i in 1:ncol(x)){
				res[,i] = (x[,i]- mi[i])/sigma[i]
				}
	} else {
		res = (x - mi) / sigma
	}
	# clean memory
	#cleanup(res)
	# results
	res
}
