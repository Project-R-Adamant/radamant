#######################################################################################################################
# FUNCTION: .getLogCounter
#
# SUMMARY:
# Internal function. Returns the number of log records currently stored into the log buffer.
#
# PARAMETERS:
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
#
# RETURNS:
# Returns the number of log records currently stored into the log buffer.
#
#
#######################################################################################################################
".getLogCounter" =  function(env = getOption("RAdamant")) {
      get("LogCounter", env);
}
#######################################################################################################################
# FUNCTION: .updateLogCounter
#
# SUMMARY:
# Internal function. Update the counter for the log buffer.
#
# PARAMETERS:
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
# - ...: Additional parameters passed to the flushLogBuffer function. 
#
# RETURNS:
# Void
#
#
#######################################################################################################################
".updateLogCounter" = function(env = getOption("RAdamant"), ...) {
	if(.getLogCounter(env) == getLogBufferSize(env)) {
		# Flush log buffer
		flushLogBuffer(env = env,...);
		# Reset log counter
		evalq({LogCounter = 1;}, env)
	} else {
		evalq({LogCounter = LogCounter + 1;}, env)
	}
}
######################################################
".BS.price.std" <- function(under, strike, rfr, sigma, maty, yield, opt.type=c("call","put"))  
{ 
	# differential factor 1
	d1 = (log(under/strike) + (rfr-yield + 0.5*sigma^2)*maty) / (sigma*sqrt(maty));
	# differential factor 2
	d2 = d1 - sigma*sqrt(maty);
	# different calculation for "call" and "put" options
	if (match.arg(opt.type)=="call"){
		# call option
		res = under*exp(-yield*maty)*pnorm(d1) - strike*exp(-rfr*maty)*pnorm(d2);
	} else {
		# put option
		res = under*exp(-yield*maty)*pnorm(-d1) - strike*exp(-rfr*maty)*pnorm(-d2) - under*exp(-maty)+strike*exp(-rfr*maty);
	}
	# return results
	Results = cbind(Price = res, Diff_1 = d1, Diff_2 = d2);
	Results;
}
######################################################
#### Black & Scholes - LogNormal ####
".BS.price.lgn" <- function(under, strike, rfr, sigma, maty, yield, opt.type=c("call","put"))  
{ 
	MOMS = BS.moments(NULL, under, rfr, sigma, yield, maty)
	# differential factor 2
	d2 = (MOMS[3,1] - log(strike)) / sqrt(MOMS[4,1])
	# differential factor 1
	d1 = d2 + sqrt(MOMS[4,1])
	# different calculation for "call" and "put" options
	if (match.arg(opt.type)=="call"){
		# call option
		res = exp(-rfr*maty) * (MOMS[1,1] * pnorm(d1)-strike*pnorm(d2))
	} else {
		# put option
		res = exp(-rfr*maty) * (-1) * (MOMS[1,1] * pnorm(-d1)-strike*pnorm(-d2))
	}
	# return results
	Results = cbind(Price_LGN = res, Diff_1 = d1, Diff_2 = d2);
	Results;
}
######################################################
#### Black & Scholes - Gamma reciprocal ####
".BS.price.gamr" <- function(under, strike, rfr, sigma, maty, yield, opt.type=c("call","put"))  
{ 
	# calculate moments
	MOMS = BS.moments(NULL, under, rfr, sigma, yield, maty)
	# gamma parameters
	alpha = (2*MOMS[2,1] - MOMS[1,1]^2) / (MOMS[2,1] - MOMS[1,1]^2)
	beta = (MOMS[2,1] - MOMS[1,1]^2) / (MOMS[2,1]*MOMS[1,1])
	# different calculation for "call" and "put" options
	if (match.arg(opt.type)=="call"){
		# differential factor 1
		d1 = pgamma(1/strike, shape=alpha-1, rate=beta, scale=beta, lower.tail=TRUE)
		# differential factor 2
		d2 = pgamma(1/strike, shape=alpha, rate=beta, scale=beta, lower.tail=TRUE)
		# call option
		res = exp(-rfr*maty) * (MOMS[1,1] * (d1) - strike * (d2))
	} else {
		# differential factor 1
		d1 = 1-pgamma(1/strike, shape=alpha-1, rate=beta, scale=beta, lower.tail=TRUE)
		# differential factor 2
		d2 = 1-pgamma(1/strike, shape=alpha, rate=beta, scale=beta, lower.tail=TRUE)		
		# put option
		res = exp(-rfr*maty) * (-1) * (MOMS[1,1] *(-d1) - strike * (-d2))
	}
	# return results
	Results = cbind(Price_GR = res, Diff_1 = d1, Diff_2 = d2)
		Results;
}
# BS formula
".BS.formula" <- function(type=c("call","put")){ 
	# different calculation for "call" and "put" options
	if (match.arg(type) == "call") 
		res = expression(under*exp(-yield*maty)*pnorm(log(under/strike) + (((rfr-yield) + 0.5*sigma^2)*maty) / (sigma*sqrt(maty))) - strike*exp(-rfr*maty)*pnorm(log(under/strike) + (((rfr-yield) + 0.5*sigma^2)*maty) / (sigma*sqrt(maty))- sigma*sqrt(maty)))
	else
		res = expression(under*exp(-yield*maty)*pnorm(-(log(under/strike) + (((rfr-yield) + 0.5*sigma^2)*maty) / (sigma*sqrt(maty)))) - strike*exp(-rfr*maty)*pnorm(-(log(under/strike) + (((rfr-yield) + 0.5*sigma^2)*maty) / (sigma*sqrt(maty))- sigma*sqrt(maty))) - under*exp(-yield*maty)+strike*exp(-rfr*maty))

	res
}
# Edgeworth factor
".EdgeFact" <- function(x, s, k){
	# calculate Edgeworth factor
	1 + s * (x^3-3*x) + (k-3) * (x^4-6*x^2+3)/24 + s^2*(x^6-15*x^4+45*x^2-15)/72
}
# Peizer-Pratt Inversion formula
".InvPP" <- function(z, n){
	# Peizer-Pratt Inversion formula
	0.5 + sign(z) * sqrt(0.25 - 0.25*exp( -(z/(n+(1/3)+(0.1/(n+1))) )^2 * (n+(1/6))))
}
# Binomial coefficient
".BinCoef" <- function(N, n){
	# calcualte binomial coefficient
	factorial(N) / (factorial(N-n)*factorial(n))
}
#### Simulate GARCH process ###
".Garch.proc" <- function(pars, order, res, type=c("garch", "mgarch", "tgarch", "egarch"), r, prob){
	#browser()
	# check names of order vector 	
	if(is.null(names(order)))
		names(order) = c("alpha","beta")
	sig0 = mean(res[-(1:max(order))]^2, na.rm=TRUE) 
	res[1:max(order)] = sqrt(sig0)
	# get parameters
	omega = pars[1]
	alpha = pars[2:(order["alpha"]+1)]
	beta = pars[(order["alpha"]+2) : (order["alpha"]+order["beta"]+1)]
	if((match.arg(type) != "garch"))
		phi = pars[length(pars)-2]
	# initial matrix with initial sample variance
	fitt = matrix(NA, length(res)+max(order), 2)
	fitt[1:max(order),] = (sig0)
	if(match.arg(type) == "garch" || match.arg(type) == "mgarch"){
		i = max(order) + 1
		while(i <= length(res) + max(order)){
			# epsilon calculation - arch component
			fitt[i, 1] = omega + alpha %*% res[(i-1):(i - order["alpha"])]^2
			# sigma calculation - garch componentn
			fitt[i, 2] = fitt[i, 1] + beta %*% (fitt[(i-1):(i - order["beta"]), 2])
			i = i + 1;
		};
	} else if(match.arg(type) == "tgarch") {
		## TGARCH
		gb = as.numeric(res<0);
		# asymetry factor
		asym = abs(filter(res^2*gb, c(0,phi), method="c", sides=1))
		asym[1] = sig0
		# arch component
		fit = omega + filter(res^2, c(0,alpha), method="c", sides=1) + asym;
		fit[1] = sig0
		fitt[1:length(fit), 1] = fit
		# garch component
		fit = filter(fit, beta, method="r", init=rep(0, order["beta"]))
		fitt[1:length(fit) ,2] = fit
	} else {
		## EGARCH
		exval = switch(prob,
			"norm" = (sqrt(2/pi)) ,
			"ged" = (gamma(2/r)/sqrt(gamma(1/r)*gamma(3/r))),
			"t" = (sqrt((r-2))*gamma((r-1)*0.5)/gamma(0.5)*gamma(r*0.5)) 
			);   
		# run egarch equation and calculate epsilon and sigma2 separately
		i = 1 + max(order)
		while(i <= length(res) + max(order)){
			# epsilon calculation - arch component
			fitt[i, 1] = omega + alpha %*% ((abs(res[(i-1):(i-max(order))])/sqrt(exp(fitt[(i-1):(i-max(order)), 1]))) - exval) + phi * (res[i-1] / sqrt(exp(fitt[i-1, 1])))
			# sigma calculation - garch component
			fitt[i, 2] = fitt[i, 1] + beta %*% fitt[(i-1):(i-max(order)), 2]
			i = i + 1;
		};
	}
	# return new fitted series Epsilon and Sigma
	fitt[-(1:max(order)),]
}
#######################################################################################################################
# FUNCTION: .genmovav
#
# SUMMARY:
# Generic Moving Average (MA filter). Computes a FIR filtering on each column of the input data
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - weights: Vector of FIR coefficients. equivalent to the 'filter' parameter of the filter function (DEFAULT = 1). 
# - padding: Padding value to fill transient of result (output data rows from 1 to win.size-1). (DEFAULT = 0)
# - rm.transient: LOGICAL. If TRUE, transient is removed, otherwise func is applied to the transient. (DEFAULT = FALSE)
# - normalize.weights: LOGICAL. If TRUE, FIR coefficients are normalized to get sum(weights) = 1 (DEFAULT = FALSE) 
# - type: Charachter attribute attached to the result (DEFAULT: "MA")
#
# RETURNS:
#  A object of class 'ma' with attributes 'type' and 'weights' as given by the corresponding input parameters:
#  - matrix of size NROW(X) by NCOL(X) where each column is the moving average of the corresponding column of X.
#   
#######################################################################################################################
".genmovav" <- function(X, weights = 1, padding = 0, rm.transient = FALSE, normalize.weights = FALSE, type = "MA", desc = "Moving Average", plot = FALSE, ...) {
	# Window size
	win.size = length(weights);
	# Check if input is an instance of the Financial Series class
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Y = X;
		# Process Close data
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		colnames(X) = attr(Y, "SName");
	}
	# Data dimensions
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Declare output
	res = matrix(padding, nrow = N, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = get.row.names(X);
	# Normalize weights to have unitary norm
	if(normalize.weights)
		weights = weights / sum(weights);
	# Apply moving average (filter)
	res[, ] = filter(X
					, filter = weights
					, sides = 1
					, method = "convolution"
					);
	# Filter function does not calculate transient (data points 1:win.size)
	transient.idx = seq(ifelse(rm.transient[1], win.size, 1), win.size-1, len = ifelse(rm.transient[1], 0, win.size-1));
	Tlen = length(transient.idx);
	# Apply moving window filter to each single serie separately
	v = 0;
	while(v < V) {
		v = v + 1;
		n = 0;
		while(n < Tlen) {
			n = n + 1;
			i = transient.idx[n];
			res[i, v] = sum(weights[1:i] * X[i:max(i-win.size+1, 1), v], na.rm = TRUE);
		}
	}
	class(res) = "movav";
	attr(res, "weights") = weights;
	attr(res, "type") = type;
	attr(res, "desc") = desc;
	# Plot Results if required
	if(plot)
		plot(res
			, X = if(fs.flag) Y else X
			, ...
			);
	# Cleanup memory
	cleanup(keep = "res");
	# Return result
	res
}
###################
# different likelihood formulas for Normal, GED and STD-t
".garch.like" <- function(X, S, prob=c("norm","ged","t"), r){
	prob = match.arg(prob)
	# Lambda parameter for GED
	if(prob == "ged")
		lambda = sqrt( (2^(-2/r)) * gamma(1/r)/gamma(3/r));
	# switch likelihood type according to prob: norm, ged, stdT
	l = switch(prob,
			# normal 	
			"norm" = 
			 0.5 * (log(S) + (log(2*pi) + (X^2)/(S)))
			,
			# ged
			"ged" = 
			-( -0.5*log(S) +
			log(r/lambda) - 0.5*abs((X)/(sqrt(S)*lambda))^r - (1+(1/r))*log(2) - log(gamma(1/r)) )
			,
			#stdT
			"t" = 
			- (-0.5*log(S) + 
			log(gamma((r+1)/2) / (sqrt(pi*(r-2))*gamma(r/2))) - ((r+1)/2) * log(1+(((X^2)/S)/(r-2))) ) 
			);
	# return likelihood value
	sum(l, na.rm=TRUE)	
}
######### transform factorise data to weight of evidence #######
".factor2woe" <- function(segm, woe){
	RES = matrix(0, NROW(segm), NCOL(segm))
	colnames(RES) = colnames(segm)
	var = 1
	while(var <= NCOL(segm)){
		idx = which(!is.na(match(woe[ ,1], colnames(segm)[var])))
		for(i in (idx))
			RES[segm[ ,var] == woe[i ,2], var] = as.numeric(woe[i, 10])
		var = var + 1
	}
	invisible(RES);
}
