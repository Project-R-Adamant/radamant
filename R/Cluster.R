###
TSClust = function(x, ...) UseMethod("TSClust")

TSClust.default = function(X,y=NULL, n_clust=5, bk.type=c("quantile","volatility","uniform","custom"), pc_vol=0.1, win.size=10,  custom_breaks=NULL, lab.dig=0){
	
	# check for NAs values
	if(any(is.na(X)))
		X = X[-is.na(X)] 
	# convert input to matrix
	if(!is.matrix(X))
		X = as.matrix(X)

	# base to split the time series
	if (is.null(y)) 
		y = (1:NROW(X))

	n = length(y)

	# select method for getting the breaks
	bk.type = match.arg(bk.type)
	
	if(bk.type=="volatility")
		warning("If cluster type is 'volatility', the 'n_clust' is ignored!")

	switch(bk.type,
	
		# quantile method
		quantile = (bk = quantile(y, probs=seq(0,1,1/n_clust), type=8)) ,
		
		# uniform method
		uniform = 
			{		
				st=round(n/n_clust, 0)
				bk= c( y[order(y)][seq(1,n,st)], max(y))			},
				
		# user defined
		custom = 
			{
				if (is.null(custom_breaks)){ 
					cat("Breaks needed!")
					return(NULL) 
				} else {
					bk = as.numeric(custom_breaks)
				}
			},
			
		# volatility clusters
		volatility = 
			{
				# calculate moving standard deviation
				msd <- movSd(X, win.size, na.rm=TRUE)

				# calculate cumulative standard deviation
				csd <- cumSd(X, 0, na.rm=TRUE)

				# set limit above average volatility
				avs <- mean(msd,na.rm=TRUE) + (mean(msd,na.rm=TRUE) * pc_vol)

				# check volatility value above the limit
				check = c(0, which(msd>=avs), 0)
				if(any(is.na(check)))
					check = check[!is.na(check)]

				# get position of the extreme of the segments
				pos = rep(FALSE, length(check))
				i = 2
				while(i < length(pos)){
	
					if(check[i] != (check[i+1]-1) || check[i] != (check[i-1]+1) )
					pos[i] = TRUE
		
					i = i + 1
				}
				# set breaks
				bk = check[pos]
			}
				
	)
	
	# round breaks
	bk = round(bk, lab.dig)
	
	# calculate cluster for time series
	bb = as.vector(findInterval(y, bk+1, rightmost.close=TRUE, all.inside=TRUE))

	# allocate list of clusters
	clust = vector("list", length(bk)-1)
	sep = unique(bb)
	sp = 1:(length(sep)-1)
	
	# put each cluster of the time series in a list
	i = 1
	while(i <= (length(bk)-1)){
		clust[[i]] = X[which(bb == sep[i]), ]
		i = i + 1
	}

	# create and assign names to each element of the list
	lab = rep(0, length(bk)-1)
	for(i in 1:(length(bk)-1)){
		if(i == 1)
			lab[i] = paste( bk[i], bk[i+1], sep=" - " )
		else
			lab[i] = paste( bk[i]+1, bk[i+1], sep=" - " )
	}
	# assign names
	names(clust) = lab
	
	# assig class and attributes
	if(bk.type=="volatility")
		attr(clust, "Vol_info") = cbind(MovSD = msd, CumSd = csd)
	attr(clust, "Type") = bk.type
	attr(clust, "Dim") = dim(X)
	attr(clust, "Breaks") = bk
	class(clust) = "TSClust"
	clust	
}


summary.TSClust = function(cluster, funs = summary, ...){
	
	# check object class
	if(class(cluster) != "TSClust"){
		cat("'~' Provide an object of class \"TSClust\" \n")
		return(NULL)
		}

	# number of clusters
	nc = length(cluster)
	
	# apply specified function to each cluster
	res = 	(lapply(cluster, funs, ...))
	
	# assign class and attributes 
	attr(res, "Call") = deparse(substitute(funs))
	class(res) = "TSClust"
	
	# return results
	res;
}


### plot TS custers ###
plot.TSClust = function(cluster, smooth=FALSE, ...){

	# check object class
	if(class(cluster) != "TSClust"){
		cat("'~' Provide an object of class \"TSClust\" \n")
		return(NULL)
		}
	
	# smooth series plot
	if(smooth){
		X = sma(unlist(cluster), na.rm=TRUE, ...) 
	} else {
		X = unlist(cluster)
	}
		
	# get breaks
	bk = attr(cluster, "Breaks")
	
	# special plot for volatility cluster
	if(attr(cluster, "Type") == "volatility"){
		par(mfrow=c(2,1))
		plot(X,type="l")
		abline(v = bk, col="red")
		plot(attr(cluster,"Vol_info")[,1], type="l")
		abline(v = bk, col="red")
	
	} else {
		# plot for regular clusters
		plot(X, type="l")
		abline(v = bk, col="red")
	} 
}



