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
TSClust = function(x, ...) UseMethod("TSClust")
TSClust.default = function(x, y=NULL, n_clust=5, bk.type=c("quantile","volatility","uniform","custom"), pc_vol=0.1, win.size=10,  custom_breaks=NULL, lab.dig=0, ...){
	# check for NAs values
	Logger(message = "check for NAs values", from = "TSClust.default", line = 2, level = 1);
	if(any(is.na(x)))
		x = x[-is.na(x)] 
	# convert input to matrix
	Logger(message = "convert input to matrix", from = "TSClust.default", line = 5, level = 1);
	if(!is.matrix(x))
		x = as.matrix(x)
	# base to split the time series
	Logger(message = "base to split the time series", from = "TSClust.default", line = 8, level = 1);
	if (is.null(y)) 
		y = (1:NROW(x))
	n = length(y)
	# select method for getting the breaks
	Logger(message = "select method for getting the breaks", from = "TSClust.default", line = 12, level = 1);
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
				msd <- movSd(x, win.size, na.rm=TRUE)
				# calculate cumulative standard deviation
				csd <- cumSd(x, 0, na.rm=TRUE)
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
	Logger(message = "round breaks", from = "TSClust.default", line = 49, level = 1);
	bk = round(bk, lab.dig)
	# calculate cluster for time series
	Logger(message = "calculate cluster for time series", from = "TSClust.default", line = 51, level = 1);
	bb = as.vector(findInterval(y, bk+1, rightmost.close=TRUE, all.inside=TRUE))
	# allocate list of clusters
	Logger(message = "allocate list of clusters", from = "TSClust.default", line = 53, level = 1);
	clust = vector("list", length(bk)-1)
	sep = unique(bb)
	sp = 1:(length(sep)-1)
	# put each cluster of the time series in a list
	Logger(message = "put each cluster of the time series in a list", from = "TSClust.default", line = 57, level = 1);
	i = 1
	while(i <= (length(bk)-1)){
		clust[[i]] = x[which(bb == sep[i]), ]
		i = i + 1
	}
	# create and assign names to each element of the list
	Logger(message = "create and assign names to each element of the list", from = "TSClust.default", line = 63, level = 1);
	lab = rep(0, length(bk)-1)
	for(i in 1:(length(bk)-1)){
		if(i == 1)
			lab[i] = paste( bk[i], bk[i+1], sep=" - " )
		else
			lab[i] = paste( bk[i]+1, bk[i+1], sep=" - " )
	}
	# assign names
	Logger(message = "assign names", from = "TSClust.default", line = 71, level = 1);
	names(clust) = lab
	# assig class and attributes
	Logger(message = "assig class and attributes", from = "TSClust.default", line = 73, level = 1);
	if(bk.type=="volatility")
		attr(clust, "Vol_info") = cbind(MovSD = msd, CumSd = csd)
	attr(clust, "Type") = bk.type
	attr(clust, "Dim") = dim(x)
	attr(clust, "Breaks") = bk
	class(clust) = "TSClust"
	clust	
}
### summary TS custers ###
summary.TSClust = function(object, funs = summary, ...){
	# check object class
	Logger(message = "check object class", from = "summary.TSClust", line = 2, level = 1);
	if(class(object) != "TSClust"){
		cat("'~' Provide an object of class \"TSClust\" \n")
		return(NULL)
		}
	# number of clusters
	Logger(message = "number of clusters", from = "summary.TSClust", line = 7, level = 1);
	nc = length(object)
	# apply specified function to each x
	Logger(message = "apply specified function to each x", from = "summary.TSClust", line = 9, level = 1);
	res = 	(lapply(object, funs, ...))
	# assign class and attributes 
	Logger(message = "assign class and attributes ", from = "summary.TSClust", line = 11, level = 1);
	attr(res, "Call") = deparse(substitute(funs))
	class(res) = "TSClust"
	# return results
	Logger(message = "return results", from = "summary.TSClust", line = 14, level = 1);
	res;
}
### plot TS custers ###
plot.TSClust = function(x, smooth=FALSE, ...){
	# check object class
	Logger(message = "check object class", from = "plot.TSClust", line = 2, level = 1);
	if(class(x) != "TSClust"){
		cat("'~' Provide an object of class \"TSClust\" \n")
		return(NULL)
		}
	# smooth series plot
	Logger(message = "smooth series plot", from = "plot.TSClust", line = 7, level = 1);
	if(smooth){
		X = sma(unlist(x), na.rm=TRUE, ...) 
	} else {
		X = unlist(x)
	}
	# get breaks
	Logger(message = "get breaks", from = "plot.TSClust", line = 13, level = 1);
	bk = attr(x, "Breaks")
	# special plot for volatility x
	Logger(message = "special plot for volatility x", from = "plot.TSClust", line = 15, level = 1);
	if(attr(x, "Type") == "volatility"){
		par(mfrow=c(2,1))
		plot(X,type="l")
		abline(v = bk, col="red")
		plot(attr(x,"Vol_info")[,1], type="l")
		abline(v = bk, col="red")
	} else {
		# plot for regular clusters
		Logger(message = "plot for regular clusters", from = "plot.TSClust", line = 23, level = 1);
		plot(X, type="l")
		abline(v = bk, col="red")
	} 
}
