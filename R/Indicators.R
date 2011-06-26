#######################################################################################################################
# FUNCTION: Price Volume Trend
#
# SUMMARY:
# This function computes the Price Volume Trend used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - Volume: vector of volumes
# - lag: lags
#
# RETURNS:
#  List of results with index and trend
#######################################################################################################################	
pvt = function(Close, Volume, lag=5, plot=FALSE, ...){

	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		Volume = X[, "Volume", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}

	adV = Volume[-(1:lag)]
	l = length(adV)
	res = rep(0,l)
	ret = Ret(Close, lag, FALSE, TRUE)
	res[1] = Volume[1]

	i = 2
	while(i <= l){
		res[i] = res[i-1] + (adV[i] * ret[i])
		i = i+1
		}

	Results = list(Price_Volume_index = res,
			Price_Volume_trend = cumsum(res)
			)
			
	class(Results) = "oscil";
	attr(Results, "type") = "PVT";

	if(plot){
		name = deparse(substitute(Close))
		main = paste("Price_Volume_index: ", name, " - ", "Lag_", lag, sep="")
		plot.oscil(Results[[1]], Close, main = main, ...)
	}	
	
	Results
}


#######################################################################################################################
# FUNCTION: Polarized Fractal Efficency
#
# SUMMARY:
# This function computes the Polarized Fractal Efficency used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - corr_fact: correction factor
# - lag: lags
#
# RETURNS:
#  Vector of results 
#######################################################################################################################	
pfe = function(X, lag=9, corr_fact=200, plot=FALSE, ...){
	
	if(class(Close) == "fs") {	
		X = Y;
		Y = X[, "Close", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}

	# series name
	name = deparse(substitute(X))
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	}

	lc = length(X)
	
	signal = ifelse(Diff(X, 1) < 0, -1, 1)

	pf = signal * 100 * sqrt( (Diff(X, lag)^2 + lag^2)) / X

	res = ema(pf[!is.na(pf)], lag)

	class(res) = "oscil";
	attr(res, "type") = "PFE";
	
	if(plot){
		main = paste("Polarized_Fractal_Efficency: ", name, " - ", "Lag_", lag, sep="")
		plot.oscil(res, Close, main = main, ...)
	}	

	res
}


#######################################################################################################################
# FUNCTION: Buying Pressure
#
# AUTHOR: FM
#
# SUMMARY:
# This function computes the Buying Pressure indicator used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - Low: vector of low prices
# - lag: lags
#
# RETURNS:
#  Vector of results 
#######################################################################################################################	
buypre=function(Close, Low,  lag=5, plot=FALSE, ...){
	
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}

	res = Close - tlow(Low, Close, lag)
	
	class(res) = "oscil";
	attr(res, "type") = "BUYPRE";
	
	if(plot){
		name = deparse(substitute(Close))
		main = paste("Buying_Pressure: ", name, " - ", "Lag_", lag, sep="")
		plot.oscil(res, Close, main = main, ...)
	}	

	res
}


#######################################################################################################################
# FUNCTION: Absolute Relative Strenght
#
# SUMMARY:
# This function computes the Absolute Relative Strenght used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - lag: lags
#
# RETURNS:
#  Vector of results 
#######################################################################################################################	
absrs = function(X, lag=14, na.rm=FALSE, plot=FALSE, ...){
	
	if(class(Close) == "fs") {	
		X = Y;
		Y = X[, "Close", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}

	# series name
	name = deparse(substitute(X))
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	}	

	# get vector of performances
	ret = Perf(X, 1, FALSE)
	
	# calculate moving averages on two conditional vectors
	greater = ema(ret[ret>=0], lag)
	lower = ema(ret[ret<0], lag)

	# fill the gap length between vectors with NAs
	cle = length(greater) - length(lower)
	if(cle<0){
		greater = c( rep(NA,abs(cle)), greater)
	} else {  
		lower = c( rep(NA,abs(cle)), lower)
	}	

	if(na.rm)
		res = (greater / (-lower))[-(1:lag),,drop=FALSE]
	else
		res = greater / (-lower)
	
	class(res) = "oscil";
	attr(res, "type") = "ABSRS";

	if(plot){
		main = paste("Absolute_Relative_Strenght: ", name, " - ", "Lag_", lag, sep="")
		plot.oscil(res, X, main = main, ...)
	}	

	res
}


#######################################################################################################################
# FUNCTION: Relative Strenght Index
#
# SUMMARY:
# This function computes the Relative Strenght Index used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
#
# RETURNS:
#  Vector of results 
#######################################################################################################################	
rsi = function(X, plot=FALSE, ...){

	if(class(Close) == "fs") {	
		X = Y;
		Y = X[, "Close", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}

	# series name
	name = deparse(substitute(X))
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	}	

	# calculate index based on absolute relative strenght
	res = 100 * ( 1 / (1 + rs(x)) )

	class(res) = "oscil";
	attr(res, "type") = "ABSRS";

	if(plot){
		main = paste("Relative_Strenght_index ", name, sep="")
		plot.oscil(res, X, main = main, ...)
	}	

	res
}


#######################################################################################################################
# FUNCTION: Mass Index
#
# SUMMARY:
# This function computes the Mass Index used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - Low: vector of low prices
# - lag: lags
#
# RETURNS:
#  Vector of results
#######################################################################################################################	
mass = function(High, Low, Close=NULL, lag=9, plot=FALSE, ...){

	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}

	hl = High - Low;

	res = ema(hl, lag) / ema((ema(hl,lag)), lag);

	class(res) = "oscil";
	attr(res, "type") = "MASS";
	
	if(plot){
		name = deparse(substitute(Close))
		main = paste("Mass_index: ", name, " - ", "Lag_", lag, sep="")
		plot.oscil(res, Close, main = main, ...)
	}	

	res
}

#######################################################################################################################
# FUNCTION: Full Price
#
# SUMMARY:
# This function computes the Full Price used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - Low: vector of low prices
# - High: vector of high prices
# - Open: vector of open prices
# - lag: lags
#
# RETURNS:
#  Vector of results
#######################################################################################################################	
fullP=function(Close, Open, High, Low, plot=FALSE, ...){

	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		Open = X[, "Open", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}

	res = (Open + High + Low + Close)/4;

	if(plot){
		name = deparse(substitute(Close))
		main = paste("Full_Price: ", name, sep="")
		cplot(res, main = main, ...)
	}	

	res
}


#######################################################################################################################
# FUNCTION: Typical Price
#
# SUMMARY:
# This function computes the Typical Price used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - Low: vector of low prices
# - High: vector of high prices
#
# RETURNS:
#  Vector of results
#######################################################################################################################	
tyP=function(Close, High, Low, plot=FALSE, ...){

	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	};

	res = (High + Low + Close) / 3;

	if(plot){
		name = deparse(substitute(Close))
		main = paste("Typical_Price: ", name, sep="")
		cplot(res, main = main, ...)
	};	

	res;
}


#######################################################################################################################
# FUNCTION: Heikin Ashi Technique
#
# SUMMARY:
# This function computes the Heikin Ashi Technique used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - Low: vector of low prices
# - High: vector of high prices
# - Open: vector of open prices
#
# RETURNS:
#  Matrix of results with "New Open", "New Close", "New Max" and "New Min"
#######################################################################################################################	
he_as = function(Close, Open, High, Low, plot=FALSE, ...){
	
	if(class(Close) == "fs") {
		X = Close	
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		Open = X[, "Open", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}

	
	res = matrix(NA, length(Close), 4);
	colnames(res)=c("New Open","New Close","New Max","New Min");

	res[,2] = fullP(Open, High, Low, Close)	
	res[,1] = Lag(Open + Close, 1) / 2		
	res[,3] = apply(cbind(High, res[,1], res[,2]), 1, max)
	res[,4] = apply(cbind(Low, res[,1], res[,2]), 1, min)
	
	class(res) = "oscil";
	attr(res, "type") = "HEAS";
	
	if(plot){
		name = deparse(substitute(Close))
		main = paste("Heikin_Ashi:", name)
		cplot(cbind(res, Close), main = main, ...)
	};
	
	res;
}



#######################################################################################################################
# FUNCTION: Keltner Channel
#
# SUMMARY:
# This function computes the Keltner Channel used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - Low: vector of low prices
# - High: vector of high prices
#
# RETURNS:
#  Matrix of results with Lower band, Upper band and Close price
#######################################################################################################################	
kelt = function(Close, High, Low, mult=2, plot=FALSE, ...){
	
	if(class(Close) == "fs") {
		X = Close	
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}

	res = matrix(0, length(Close), 3);
	colnames(res) = c("Lower", "EMA_20 Close_Price", "Upper");

	# compute average true range
	av_tr = trf(Close, Low, High, 10)
	# compute expoential moving average
	smooth = ema(Close, 20)	

	res[,1] = smooth - (mult * av_tr) 
	res[,2] = smooth
	res[,3] = smooth + (mult * av_tr)
	
	class(res) = "oscil";
	attr(res, "type") = "kelt";
	
	if(plot){
		name = deparse(substitute(Close))
		main = paste("keltner_Channel: ", name, " - ", "Mult_", mult, sep="")
		cplot(cbind(res, Close), main = main, ...)
	};
	
	res;
}


#######################################################################################################################
# FUNCTION: Elder-Ray Force
#
# SUMMARY:
# This function computes the Elder-Ray Force index used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - Low: vector of low prices
# - High: vector of high prices
# - lag: lag period
#
# RETURNS:
#  Matrix of results with "Bear_Power", "Close_Price" and "Bull_Power"
#######################################################################################################################	
erf = function(Close, High=NULL, Low=NULL, lag=13, plot=FALSE, ...){
	
	if(class(Close) == "fs") {
		X = Close	
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	
	res = matrix(0,length(Close),3);
	colnames(res)=c("Bear_Power", "Close_Price", "Bull_Power");

	# Bull Power
	res[,3] = High - ema(Close, lag)
	# Bear_Power
	res[,1] = Low - ema(Close, lag)
	# Closing price
	res[,2] = Close

	class(res) = "oscil";
	attr(res, "type") = "kelt";
	
	if(plot){
		name = deparse(substitute(Close))
		main = paste("Elde_Ray_Force_index: ", name, " - ", "Lags_", lag, sep="")
		plot.oscil(res[,c(1,3)], Close, main = main, ...)
	};
	
	res;
}


#######################################################################################################################
# FUNCTION: Elder Force Index
#
# SUMMARY:
# This function computes the Elder Force Index used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - Volume: vector of volumes
# - lag: lag period
#
# RETURNS:
#  Vector of results 
#######################################################################################################################	
erfi = function(X, Volume, lag=13, plot=FALSE, ...){

	if(class(X) == "fs") {
		Y = X	
		Close = Y[, "Close", drop = FALSE];
		Volume = Y[, "Volume", drop = FALSE];
		colnames(Close) = attr(Y, "SName");
	}

	# force index
	res = Volume * Diff(X, lag);

	class(res) = "oscil";
	attr(res, "type") = "erfv";
	
	if(plot){
		name = deparse(substitute(Close))
		main = paste("Elde_Force_Volume_index: ", name, " - ", "Lags_", lag, sep="")
		plot.oscil(res, Close, main = main, ...)
	};
	
	res;
}


#######################################################################################################################
# FUNCTION: Chaikin Money Flow
#
# SUMMARY:
# This function computes the Chaikin Money Flow used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - Low: vector of low prices
# - High: vector of high prices
# - Volume: vector of open prices
#
# RETURNS:
#  Vector of results 
#######################################################################################################################	
cmf=function(Close, Low, High, Volume, plot=FALSE, ...){
	
	if(class(Close) == "fs") {
		X = Close	
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}

	res = (clv(Close, Low, High) * Volume) / sma(Volume);
	
	class(res) = "oscil";
	attr(res, "type") = "CMF";
	
	if(plot){
		name = deparse(substitute(Close))
		main = paste("Chaikin_Money_Flow: ", name, sep="")
		plot.oscil(res, Close, main = main, ...)
	};
	
	res;
}


#######################################################################################################################
# FUNCTION: Chaikin Volatility
#
# SUMMARY:
# This function computes the Chaikin Volatility used for technical analysis
#
# PARAMETERS:
# - Low: vector of low prices
# - High: vector of high prices
# - lag: lags period
#
# RETURNS:
#  Vector of results 
#######################################################################################################################	
Ch.vol=function(High, Low, Close, lag=5, plot=FALSE, ...){

	if(class(High) == "fs") {
		X = High	
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		close = X[, "close", drop = FALSE];
		colnames(High) = colnames(Low) = attr(X, "SName");
	}

	# differnce High Low
	dhl = High - Low
	# average rate of change
	aroc = dhl - Perf(dhl, 10,cut=FALSE)
	
	res = 100*( (ema(dhl, lag)/aroc)-1 );
	
	class(res) = "oscil";
	attr(res, "type") = "CHCOL";
	
	if(plot){
		name = deparse(substitute(Close))
		main = paste("Chaikin_Volatility: ", name, " - ", "Lags_", lag, sep="")
		plot.oscil(res, Close, main = main, ...)
	};
	
	res;
}


#######################################################################################################################
# FUNCTION: Commodity Channel Index
#
# SUMMARY:
# This function computes the Commodity Channel Index used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - Low: vector of low prices
# - High: vector of high prices
# - lag: lags period
#
# RETURNS:
#  Vector of results 
#######################################################################################################################	
cci=function(High, Low, Close, lag=5, plot=FALSE, ...){
	
	if(class(Close) == "fs"){
		X = Close	
		Close = X[, "Close", drop = FALSE]
		High = X[, "High", drop = FALSE]
		Low = X[, "Low", drop = FALSE]
		colnames(Close) = attr(X, "SName")
	};
	
	# typical price
	tp = tyP(High, Low, Close);

	# difference
	dif = (tp - sma(tp, lag));

	# average difference
	avdif = sma(abs(dif), lag);
	
	res = (dif / avdif) * 0.67;

	class(res) = "oscil";
	attr(res, "type") = "CCI";
	
	if(plot){
		name = deparse(substitute(Close))
		main = paste("Commodity_Channel_index: ", name, " - ", "Lags_", lag, sep="")
		plot.oscil(res, Close, main = main, ...)
	};
	
	res;
}


#######################################################################################################################
# FUNCTION: Commodity Channel Index - Version2
#
# SUMMARY:
# This function computes a different version of the Commodity Channel Index used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - Low: vector of low prices
# - High: vector of high prices
# - lag: lags period
#
# RETURNS:
#  Vector of results 
#######################################################################################################################	
cci.v2 = function(High, Low, Close, lag=5, plot=FALSE, ...){

	if(class(Close) == "fs"){
		X = Close	
		Close = X[, "Close", drop = FALSE]
		High = X[, "High", drop = FALSE]
		Low = X[, "Low", drop = FALSE]
		colnames(Close) = attr(X, "SName")
	};

	# typical price
	tp = tyP(High, Low, Close);

	m = sma(tp,lag);
	m2 = sum(abs(tp - sma(tp,lag))) / lag;

	# calculate index
	res = ((tp - m) / (0.015 * m2))[-(1:lag)]

	class(res) = "oscil";
	attr(res, "type") = "CCI2";
	
	if(plot){
		name = deparse(substitute(Close))
		main = paste("Commodity_Channel_V2_index: ", name, " - ", "Lags_", lag, sep="")
		plot.oscil(res, Close, main = main, ...)
	};
	
	res;
}


#######################################################################################################################
# FUNCTION: Chande Momentum Oscillator
#
# SUMMARY:
# This function computes the Chande Momentum Oscillator used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - Low: vector of low prices
# - High: vector of high prices
# - lag: lag period
#
# RETURNS:
#  Vector of results
#######################################################################################################################	
cmof = function (X, lag=5, plot=FALSE, ...){ 

	if(class(Close) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(Close) = attr(Y, "SName");
	}

	# series name
	name = deparse(substitute(X))
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	}	

	# lagged difference
	d = Diff(X, lag ,na.rm=TRUE)
	
	cmo1 = as.numeric(d >= 0);  #ifelse(d>=0,d,0)
	cmo2 = as.numeric(d < 0);	 #ifelse(d>=0,0,-d)	
	
	# up movement
	up = cumsum(d * cmo1);
	# down movement
	down = cumsum(-d * cmo2);

	# calculate index
	res = ((up - down) / (up + down)) * 100
	
	class(res) = "oscil";
	attr(res, "type") = "CMO";
	
	if(plot){
		main = paste("Chande_Momentum_Oscillator: ", name, " - ", "Lags_", lag, sep="")
		plot.oscil(res, X, main = main, ...)
	};

	# return results
	res;
}


#######################################################################################################################
# FUNCTION: Variable Chande Momentum Oscillator
#
# SUMMARY:
# This function computes the Variable Chande Momentum Oscillator used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - lag: lag period
#
# RETURNS:
#  Vector of results
#######################################################################################################################	
vcmof = function(X, lag=5, plot=FALSE, ...){

	if(class(Close) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(Close) = attr(Y, "SName");
	}

	# series name
	name = deparse(substitute(X))
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	}	
 	
	vr = (cmof(X, lag) / 100 )
	
	ma = sma(X, lag)[-(1:lag)]
	malag = Lag(ma, 1)
	sm.fac = 2/(lag+1)
	
	res = sm.fac * vr * X[-(1:lag)] + (1-sm.fac * vr) * malag;

	class(res) = "oscil";
	attr(res, "type") = "VCMO";
	
	if(plot){
		main = paste("Chande_Variable_Momentum_Oscillator: ", name, " - ", "Lags_", lag, sep="")
		plot.oscil(res, X, main = main, ...)
	};

	# return results
	res;
}


#######################################################################################################################
# FUNCTION: Variable Index Dynamic Average
#
# SUMMARY:
# This function computes the Variable Index Dynamic Average used for technical analysis
#
# PARAMETERS:
# - X: vector of prices
# - lag: lag window
#
# RETURNS:
#  Vector of results
#######################################################################################################################	
vidyaf = function (X, lag=5, plot=FALSE, ...){

	if(class(Close) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(Close) = attr(Y, "SName");
	};

	# series name
	name = deparse(substitute(X));
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	};	
	
	l = length(X);

	# vector of results
	res = rep(0, l-lag);
	sm.fac = 2/(lag+1);
	cmo = abs(cmof(X, lag));
	
	# starting value
	res[1] = X[1]
	
	# calculate index
	i = 2
	while(i <= l-lag)
	{		
		res[i] = X[i-1] * sm.fac * cmo[i] + vidya[i-1] * (1 - (sm.fac * cmo[i])) 
		i = i + 1
	}
	
	class(res) = "oscil";
	attr(res, "type") = "VIDA";
	
	if(plot){
		main = paste("Variable_Index_Dynamic_Average: ", name, " - ", "Lags_", lag, sep="")
		plot.oscil(res, X, main = main, ...)
	};

	# return results
	res;
}


#######################################################################################################################
# FUNCTION: Vertical Horizontal Filter
#
# SUMMARY:
# This function computes the Vertical Horizontal Filter used for technical analysis
#
# PARAMETERS:
# - x: vector of prices
# - lag: lag window
#
# RETURNS:
#  Vector of results
#######################################################################################################################	
vhff<- function (X, lag=9, plot=FALSE, ...){ 	

	if(class(Close) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(Close) = attr(Y, "SName");
	};

	# series name
	name = deparse(substitute(X));
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	};	

	dab = abs(Diff(X, 1)) 
	
	# scaled Max/Min
	hh = scalMax(x, lag)
	ll = scalMin(x, lag)

	# rate of change
	den = roc(x, lag, plot=FALSE)	
		 
	res = (hh-ll)/den;

	class(res) = "oscil";
	attr(res, "type") = "VHF";
	
	if(plot){
		main = paste("Vertical_Horizontal_filter: ", name, " - ", "Lags_", lag, sep="")
		plot.oscil(res, X, main = main, ...)
	};

	# return results
	res;
}



#######################################################################################################################
# FUNCTION: DeMarker Indicator
#
# SUMMARY:
# This function computes the DeMarker Indicator used for technical analysis
#
# PARAMETERS:
# - Low: vector of low prices
# - High: vector of high prices
# - lag: lag period
#
# RETURNS:
#  Vector of results
#######################################################################################################################	
demark = function(High, Low, Close, lag=5, plot=FALSE, ...){

	if(class(Close) == "fs") {	
		Y = Close;
		Close = Y[, "Close", drop = FALSE];
		High = Y[, "High", drop = FALSE];
		Low = Y[, "Low", drop = FALSE];
		colnames(High) = colnames(Low) = colnames(Close) = attr(Y, "SName");
	};

	demax = ifelse(High>Lag(High,1), Diff(High, 1), 0)
	demin = ifelse(Low>Lag(Low,1), Diff(Low, 1), 0)

	res = sma(demax, lag) / (sma(demax, lag) + sma(demin, lag));

	class(res) = "oscil";
	attr(res, "type") = "DEM";
	
	if(plot){
		name = deparse(substitute(Close))
		main = paste("Demark_index: ", name, " - ", "Lags_", lag, sep="")
		plot.oscil(res, Close, main = main, ...)
	};
	
	res;
}


#######################################################################################################################
# FUNCTION: Tirone Levels
#
# SUMMARY:
# This function computes the Tirone Levels used for technical analysis
#
# PARAMETERS:
# - Low: vector of low prices
# - High: vector of high prices
# - lag: lag period
#
# RETURNS:
#  Matrix of results with "Lower" , "Center" and "Upper" bands
#######################################################################################################################	
tirLev = function(High, Low, lag=5, plot=FALSE, ...){

	if(class(High) == "fs") {	
		Y = High;
		High = Y[, "High", drop = FALSE];
		Low = Y[, "Low", drop = FALSE];
		colnames(High) = colnames(Low) = attr(Y, "SName");
	};
	
	# declare output matrix
	res = matrix(0, length(High), 3)	
	colnames(res)=c("Lower","Center","Upper")

	res[,2] = (movMax(High, lag) - movMin(Low, lag)) / 2;
	res[,3] = res[,2] + 0.67*res[,2];
	res[,1] = res[,2] - 0.33*res[,2];
	
	class(res) = "oscil";
	attr(res, "type") = "TIR";
	
	if(plot){
		name = deparse(substitute(High))
		main = paste("Tirone_levels: ", name, " - ", "Lags_", lag, sep="")
		plot.oscil(res, Close, main = main, ...)
	};
	
	res;
} 
 
 #######################################################################################################################
# FUNCTION: Parabolic Stop-and-Reverse (PSAR)
#
# SUMMARY:
# This function computes the Parabolic Stop-and-Reverse used for technical analysis
#
# PARAMETERS:
# - Low: vector of low prices
# - High: vector of high prices
# - accel: acceleration factor
#
# RETURNS:
#  Vector of results
#######################################################################################################################	
prbsar = function(Close, High, Low, accel=c(0.02,0.2), plot=FALSE, ...){

	if(class(High) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	}


	## initialize variables
	# stop-and-reverse (sar) vector
	sar = rep(0, length(High))
	
	# signals
	signal1 = 0
	signal0 =1

	# extreme prices
	extP0 = High[1]
	extP1 = 0

	# first value for sar
	sar[1] = Low[1]-0.01

	# acceleration factors
	acc_fact0 = accel[1]
	acc_fact1 = 0

	# step minimum
	lmin = ifelse(Low < Lag(Low), Lag(Low), Low)
	# step maximum
	lmax = ifelse(High < Lag(High), Lag(High), High)

	i = 2
	while(i <= length(High)){

		signal1 = signal0
		extP1 = extP0
		acc_fact1 = acc_fact0

		if(signal1 == 1){
			if(Low[i] > sar[i-1])
				signal0 = 1 else
					signal0 = -1
			if(High[i] > extP1)
				extP0 = High[i] else
					extP0 = extP1
			} else {
			
			if(High[i] < sar[i-1])
				signal0 = -1 else
					signal0 = 1
			if(Low[i] < extP1)
				extP0 = Low[i] else
					extP0 = extP1
			}

		if(signal0 == signal1){

			sar[i] = sar[i-1] + (extP1 - sar[i-1]) * acc_fact1

			if(signal0 == 1){
				if(extP0 > extP1){
					if(acc_fact1 == accel[2])
						acc_fact0 = accel[2] else
							acc_fact0 = accel[1] + acc_fact1	
				} else {
					acc_fact0 = acc_fact1
				}

			if(sar[i] > lmin[i])
				sar[i] = lmin[i]
		} else {

			if(extP0 < extP1){
				if(acc_fact1 == accel[2])
					acc_fact0 = accel[2] else
						acc_fact0 =	accel[1] + acc_fact1
				} else {
				acc_fact0 = acc_fact1				
				}
			
			if(sar[i] < lmax[i])
				sar[i] = lmax[i]
		}
		
		} else {
			acc_fact0 = accel[1]
			sar[i] = extP1
		}
	i = i + 1
	};

	class(sar) = "oscil";
	attr(sar, "type") = "PRBSAR";
	
	if(plot){
		name = deparse(substitute(Close))
		main = paste("Parabolic_SAR: ", name, " - ", "Acc.fact_", accel, sep="")
		plot.oscil(Osc = sar, X = Close, ...)
	};

	sar;
		
}	


 #######################################################################################################################
# FUNCTION: Mass Index
#
# SUMMARY:
# This function computes the Mass Index used for technical analysis
#
# PARAMETERS:
# - Low: vector of low prices
# - High: vector of high prices
# - accel: acceleration factor
#
# RETURNS:
#  Vector of results
#######################################################################################################################	
mass.cum = function(High, Low, Close=NULL, lag=9, plot=FALSE, ...){  

	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}

	res = rep(NA,length(High))
	emr = ema(High - Low, lag) / ema( ema(High - Low, lag), lag)
	
	i = length(High)

	# calculated sum
	while(i>(lag-1)){
		res[i]=sum(emr[i:(i-lag)])
		i=i-1
	}

	class(res) = "oscil";
	attr(res, "type") = "MASSC";
	
	if(plot){
		name = deparse(substitute(Close))
		main = paste("Mass_Cumulative_index: ", name, " - ", "Lags_", lag, sep="")
		plot.oscil(Osc = res, X = Close, ...)
	};

	res;

}


 #######################################################################################################################
# FUNCTION: McClellan Summation Index
#
# SUMMARY:
# This function computes the McClellan Summation Index used for technical analysis
#
# PARAMETERS:
# TO DO
#
# RETURNS:
#  Vector of results
#######################################################################################################################	
mcsi = function(matr, nr, nc, lag1, lag2, plot=FALSE, ...){

	# declare output
	mc_si = rep(NA, l);
	
	# compute McClellan oscillator
	mc_osc = mcosc(matr,nr,nc,lag1,lag2);

	l = length(mcosc);
	
	# initial value 	
	mcsi[1] = mcosc[1];

	i=2
	while(i<l){
		mc_si[i] = sum(1000,mc_si[i-1],mc_osc[i],na.rm=TRUE)
		i=i+1
	};
		
	mc_si;
}
	

#MG - McGinley Dynamic Indicator:
mcgind = function(X, lag=12, plot=FALSE, ...){  

	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	};

	# series name
	name = deparse(substitute(X))
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	}

	# lagged Close
	lc = Lag(X, 1)
	
	# average lagged Close
	em = ema(lc, lag);
	  
	# compute indicator
	res = em + ((lc - em) / (lc / em) * 125);
	
	class(res) = "oscil";
	attr(res, "type") = "MCGDY";

	if(plot){
		main = paste("McGinley_Dynamic_index: ", name, " - ", "Lags_", lag, sep="")
		plot.oscil(Osc = res, X = X, main = main, ...);
	}
	
	res;	
}


#######################################################################################################################
# FUNCTION: Money Flow
#
# SUMMARY:
# This function computes the Money Flow used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - Low: vector of low prices
# - High: vector of high prices
# - Volume: vector of volumes
#
# RETURNS:
#  Vector of results
#######################################################################################################################	
Mflow = function(Close, High, Low, Volume, plot=FALSE, ...){

	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		Volume = X[, "Volume", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	};

	res = Volume * tyP(High, Low, Close);

	class(res) = "oscil";
	attr(res, "type") = "MFLOW";

	if(plot){
		name = deparse(substitute(Close))
		main = paste("Money_Flow: ", name, sep="")
		plot.oscil(Osc = res, X = Close, main=main, ...);
	}
	
	res;
}


#######################################################################################################################
# FUNCTION: Money Ratio
#
# SUMMARY:
# This function computes the Money Ratio used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - Low: vector of low prices
# - High: vector of high prices
# - Volume: vector of volumes
#
# RETURNS:
#  Vector of results
#######################################################################################################################	
Mflow.ratio = function(Close, High, Low, Volume, plot=FALSE, ...){

	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		Volume = X[, "Volume", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}

	# compute money flow
	mf = Mflow(High, Low, Close, Volume);

	# first difference of typical price
	d_typ = Diff( tyP(High, Low, Close), 1, na.rm=TRUE );

	# calculate increase of typical price
	cond_up = ifelse(d_typ>0, d_typ, 0);
	# calculate decrease of typical price
	cond_do = ifelse(d_typ<0, -d_typ, 0);  

	up = down = rep(0, length(d_typ));

	i = 1
	while(i <= length(d_typ))
	{
		up[i] = sum(cond_up[i:(i+1)], na.rm=TRUE)
		down[i] = sum(cond_do[i:(i+1)], na.rm=TRUE)
		i = i + 1
	};

	res = (up+1) / (down+1);

	class(res) = "oscil";
	attr(res, "type") = "MFLRAT";

	if(plot){
		name = deparse(substitute(Close))
		main = paste("Money_Flow_ratio: ", name, sep="")
		plot.oscil(Osc = res, X = Close, main=main, ...);
	}
	
	res;

}


#######################################################################################################################
# FUNCTION: Money Flow Index
#
# SUMMARY:
# This function computes the Money Flow Index used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - Low: vector of low prices
# - High: vector of high prices
# - Volume: vector of volumes
#
# RETURNS:
#  Vector of results
#######################################################################################################################	
Mflow.ind = function(Close, High, Low, Volume, plot=FALSE, ...){
	
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		Volume = X[, "Volume", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}

	# get money ratio
	mr = mf_ratio(High,Low,Close,Volume);
	
	# get money flow index
	res = 100 - (100 / (1+mr));
	
	class(res) = "oscil";
	attr(res, "type") = "MFLIND";

	if(plot){
		main = paste("Money_Flow_index: ", name, sep="")
		name = deparse(substitute(Close))
		plot.oscil(Osc = res, X = Close, main=main, ...);
	}
	
	res;
		
}	
	
#######################################################################################################################
# FUNCTION: Kairi Relative Index
#
# SUMMARY:
# This function computes the Kairi Relative Index used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - lag1: lags period
# - lag2: lags period
#
# RETURNS:
#  Vector of results
#######################################################################################################################	
kri = function(X, lag1=10, lag2=20, plot=FALSE, ...){
	
	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	};

	# series name
	name = deparse(substitute(X))
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	}

	# calculate index
	res = ( (X - sma(X, min(lag1, lag2))) / sma(X, max(lag1, lag2)) )[-(1:min(lag1, lag2))] * 100;

	class(res) = "oscil";
	attr(res, "type") = "KRI";

	if(plot){
		main = paste("Kairi_Relative_Index: ", name, " - ", "Lags_", lag1, "/", lag2, sep="")
		plot.oscil(Osc = res, X = X, main = main, ...);
	}
	
	res;
}	

#######################################################################################################################
# FUNCTION: Swing Index
#
# SUMMARY:
# This function computes the Swing Index used for technical analysis
#
# PARAMETERS:
# - Close: vector of close prices
# - High: vector of high prices
# - Low: vector of low prices
# - Open: vector of open prices
#
# RETURNS:
#  Vector of results
#######################################################################################################################	

Swing = function(Close, High, Low, Open, ret_cum=FALSE, plot=FALSE, ...){

	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		Open = X[, "Open", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	};

	# calculate differences needed to calculate the index
	c_lc = Diff(Close);
	h_c = High - Close;
	h_lc = High - Lag(Close);
	l_lc = Low - Lag(Close);
	lc_lo = Lag(Close) + Lag(Open);
	
	# numerator
	num =  (abs(cbind(c_lc, Close - Open, lc_lo ))) %*% c(1, 0.5, 0.25);
	
	# get different maximum values per row
	mm2 = abs(cbind(h_lc, l_lc, h_c));
	rr = cbind(h_lc + 0.5*l_lc + 0.25*lc_lo, 
				l_lc + 0.5*h_lc + 0.25*lc_lo, 
				(High - Low) + 0.25*lc_lo);
	dis = t(apply(mm2, 1, function(x) (x == max(x))));

	# denominator
	den = apply(rr * dis, 1, sum, na.rm=TRUE);
	
	fact = rowMax(cbind(h_lc, l_lc)) / vl;
	
	# calculate swing index
	res = 50 * (num / den) * fact;
	
	# return index and cumulative index
	if(ret_cum)
		res = cumSum(res, na.rm=TRUE)

	class(res) = "oscil";
	attr(res, "type") = "SWING";

	if(plot){
		name = deparse(substitute(Close))
		main = paste("Swing_index: ", name, sep="")
		plot(res, X=Close, ...)
	}
	
	res;

}



#############################################################
### START - ADVANCE DECLINE INDEXES ###
#############################################################

# Advance / Decline issues
AdvDec = function(X, lag=5, ret.idx=TRUE, plot=FALSE, ...){
	
	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	}

	# series name
	name = deparse(substitute(X))
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	}

	# number of assets	
	N = NCOL(X);
	R = NROW(X) - lag;

	# lagged difference
	dp = Diff(X, lag, na.rm=TRUE);

	# matrix of resulst
	res = matrix(0, R, 3);
	colnames(res) = c("Advance", "Decline", "Difference");

	# number of advancing issues
	res[, 1] = apply(dp, 1, function(x) length(which(x>0)));
	# number of declining issues
	res[, 2] = N - res[,1];
	# difference
	res[, 3] = res[,1] - res[,2]

	if(ret.idx){
		res2 = cumsum(res[,3])
	
		class(res2) = "oscil";
		attr(res2, "type") = "ADI";

		if(plot){
			main = paste("Advance/Decline_index: ", name, " Lag_", lag, sep="")
			plot(res2, X=X, ...)
		}

		# return index
		return(res2);
	}

	# return table
	res;

}


# Absolute breadth index
Abi = function(X, lag=5, plot=FALSE, ...){
	
	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	}

	# series name
	name = deparse(substitute(X))
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	}

	# calculate advance / decline
	ad = AdvDec(X, lag, FALSE, FALSE);
	
	# calculate abi index
	res = abs(ad[,1] - ad[,2]);

	class(res) = "oscil";
	attr(res, "type") = "ABI";

	if(plot){
		main = paste("Absolute_Breadth_index: ", name, " Lag_", lag, sep="")
		plot(res, X=X, ...)
	}
	
	res;
}

# Breadth Thrust
Breadth = function(X, lag=5, plot=FALSE, ...){

	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	}

	# series name
	name = deparse(substitute(X))
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	}

	# calculate advance / decline
	ad = AdvDec(X, lag, FALSE, FALSE);

	# calculate index
	res = sma(ad[,1]/(ad[,1]+ad[,2]), lag);

	class(res) = "oscil";
	attr(res, "type") = "BRETR";
	
	if(plot){
		main = paste("Breadth_Thrust_index: ", name, " Lag_", lag, sep="")
		plot(res, X=X, ...)
	}
	
	res;
	
}

# Advance Decline Ratio
ADratio = function(X, lag=5, plot=FALSE, ...){

	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	}

	# series name
	name = deparse(substitute(X))
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	}

	# calculate advance / decline
	ad = AdvDec(X, lag, FALSE, FALS);

	res = ad[,1] / ad[,2];

	class(res) = "oscil";
	attr(res, "type") = "ADR";

	if(plot){
		main = paste("Advance/Decline_ratio: ", name, " Lag_", lag, sep="")
		plot(res, X=X, ...)
	}
	
	res;
	
}

# ARMS Index - TRIN
Arms = function(X, Volume, lag, plot=FALSE, ...){

	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		Volume = X[, "Volume", drop = FALSE];
		colnames(X) 	 = attr(Y, "SName");
	}

	# series name
	name = deparse(substitute(X))
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	}

	# calculate trin index 
	res = ADratio(X, lag) / ADratio(XV, lag);
	
	class(res) = "oscil";
	attr(res, "type") = "TRIN";

	if(plot){
		main = paste("Arms_TRIN_index: ", name, " Lag_", lag, sep="")
		plot(res, X=X, ...)
	}
	
	res;

}

#############################################################
### END - ADVANCE DECLINE INDEXES ###
#############################################################


#HHV - Max nelle ultime x barre:
hhv = function(X, lag, na.rm=TRUE){   
	
	maax = rep(NA,length(X))
	
	i = length(X)
	
	# get rolling max
	while (i >= (lag-1)){
		maax[i] = max(X[i:(i - lag+1)])
		i = i - 1
	}
	
	# return results with or without NAs
	if(na.rm)
		maax[-(1:lag)]
	else
		maax	
}



#LLV - Min nelle ultime x barre:
llv = function(X, lag, na.rm=TRUE){   
	
	miin = rep(NA,length(X))
	
	i = length(X)
	
	# get rolling min
	while (i >= (lag-1)){
		miin[i] = min(X[i:(i - lag+1)])
		i = i - 1
	}
	
	# return results with or without NAs
	if(na.rm)
		miin[-(1:lag)]
	else
		miin
}

## PRICE CHANNELS
Pchan = function(CLose, High, Low, lag=20, na.rm=TRUE, plot=FALSE, ...){
	
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}

	# highest high
	up = hhv(High, lag, na.rm)
	# lowest low
	down = llv(Low, lag, na.rm)
	# mid trend
	mid = (up + down) / 2
	
	# matrix of results
	res = cbind(up,down,mid)
	
	# return results with or without NAs
	if(na.rm) 
		res[-(1:lag),]
	else
		res
	
	class(res) = "oscil";
	attr(res, "type") = "PCHAN";
	
	if(plot){
		name = deparse(substitute(Close))
		main = paste("Price_Channel: ", name, " - ", "Lag_", lag, sep="")
		cplot(cbind(res,Close[-(1:lag)]), main = main, ...)
	}
		
}

#####################
## Ichimoku Kinko Hyo 
#####################
Ichkh = function(Close, High, Low, plot=FALSE, ...){
	
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}

	# Tenkan-Sen
	ts = (hhv(High, lag=9) + llv(Low, lag = 9)) / 2

	# Kijun-Sen
	ks = (hhv(High, lag = 26) + llv(Low, lag = 26)) / 2

	# Chikou Span
	cs = Lag(Close,  lag=-26)

	# Senkou-Span A
	ssa = Lag((ts + ks) / 2, lag=26)

	# Senkou-Span A
	ssb = Lag((hhv(High, lag = 52) + llv(Low, lag = 52)) / 2   , lag = 26)

	# matrix of results
	res = cbind(Tenkan_sen = ts,
			Kijun_sen = ks,
			Chikou_span = cs,
			Senkou_spanA = ssa,
			Senkou_spanB = ssb)
	
	class(res) = "oscil";
	attr(res, "type") = "ICHS";
	
	if(plot){
		name = deparse(substitute(X))
		main = paste("Ichimoku_Kinko_Hyo: ", name, sep="")
		cplot(res, main = main, ...)
	}
	
	res;
}

#########################
## FORCE INDEX
#########################

forcidx = function(X, Volume, lag=5, sth=TRUE, sth.lag=13, mov=sma, plot=FALSE, ...){
	
	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		Volume = Y[, "Volume", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	}
	
	# series name
	name = deparse(substitute(X))
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	}

	# price difference
	idx = Diff(X, lag, ...) * Volume
	
	# get moving average type
	#mov =((match.arg(mov)))
	
	# smooth results with mov
	if(sth)
		res = mov(idx, lag, ...)
	else
		res = idx
		
	class(res) = "oscil";
	attr(res, "type") = "FORIDX";

	if(plot){
		main = paste("Force_index: ", name, " Lag_", lag, sep="")
		plot(res, X=X, ...)
	}
		
	# return resutls
	res

}

## ULCER INDEX ##
ulcer = function(X, lag, plot=FALSE, ...){
	
	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	}
	
	# series name
	name = deparse(substitute(X));
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	};
	
	# highest high
	hh = hhv(X, lag);
	
	# calculate index
	res = sqrt(sma(((1 - (X / hh)) * 100)^2));
	
	class(res) = "oscil";
	attr(res, "type") = "ULC";
	
	if(plot){
		main = paste("Ulcer_index: ", name, " Lag_", lag, sep="");
		plot(res, X=X, main = main, ...);
	}
	
	res
	
}

## STOLLER STARC BANDS - Stoller Average Range Channels
starc = function(Close, High=NULL, Low=NULL, atr.mult=2, lag=5, atr.lag=14, mov=c("sma","ema","wma"), plot=FALSE, ...){
	
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}

	# calculate average true range
	atr = trf(Close, High, Low, lag, avg.lag = atr.lag)[-(1:lag)];
	# smoothing type
	mov = match.arg(mov)
	# calculate moving average
	ma = do.call(mov, list(Close,lag))[-(1:lag)];

	# upper band
	upp = ma + (atr.mult * atr);
	# lower band
	low = ma - (atr.mult * atr);

	# matris of results
	res = cbind(Lower_band = low, Middle_band = ma, Upper_band = upp);
	
	class(res) = "oscil";
	attr(res, "type") = "STARC";
	
	if(plot){
		name = deparse(substitute(Close))
		main = paste("Stoller_Starc_Bands: ", name, " - ", "Lag_", lag, sep="")
		cplot(cbind(res, Close), main = main, ...)
	}

	res;

}


########################################
# Performance indicators
Perf = function(X, ini.per=1 ,cut=TRUE, plot=FALSE, ...){

	if(class(X) == "fs") {
		Y = X	
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	}
	
	# series name
	name = deparse(substitute(X))
	
	# convert to matrix and apply names
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	}
	
	# starting reference point 
	start = mean(X[1:ini.per,, drop=FALSE]);

	# calculate index
	if(cut)
		res = (100 * (X/start - 1))[-(1:ini.per),,drop=FALSE]
	else
		res = 100 * (X/start - 1)
	
	class(res) = "oscil";
	attr(res, "type") = "PERF";
	
	if(plot){	
		main = paste("Performance Indicator:", name)
		cplot(res, main=main, ...)
	}
	
	res;
	
}
