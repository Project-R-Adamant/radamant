# GARCH CLASS AND FUNCTIONS
Garch = function(x,...) UseMethod("Garch")
newsimp = function(x,...) UseMethod("newsimp")
# extract coefficients
coef.Garch = function(object, names=TRUE, ...){
	cc = as.matrix(object$Results$Estimates)	
	if(names)
		rownames(cc) = rownames(object$Results)
	cc
}
# extraact loglikelihood
logLik.Garch = function(object, ...){
	object$LogLik
}
# extract coefficients covariance matrix
vcov.Garch = function(object, ...){
	object$Vcov
}
print.Garch = function(x, digits=5, ...){
	cat(rep("=", 40), "\n",sep="")
	cat("Model output:", "\n")
	cat(" ", paste(paste(paste("(",round(x[[4]]$Estimates,5), ")", sep=""), rownames(x[[4]]), sep = "*"), collapse = " + "), "\n\n")
	cat(rep("=", 40), "\n",sep="")
	cat("Model Info:", "\n")
	cat("Type:", (x[[1]]), "\n")
	cat("Order:", (x[[2]]), "\n")
	cat("LogLik:", (x$LogLik), "\n")
	cat("Vol_Pers:", (x$Volatility_Persistence), "\n")
	cat("AIC:", (x$AIC), "\n")
	cat(rep("=", 40), "\n",sep="")
	cat("Mean equation:", "\n")
	print(x[[3]], digits=5)
	cat("\n")
	xx = apply(x$Results, 2, round, digits)
	cat(rep("=", 40), "\n",sep="")
	cat("Estimation results:", "\n")
	print(xx)
}
#######################################################################################################################
# FUNCTION: Garch,
#
# SUMMARY:
# This function estimates a Garch(p ,q) model;
#
# PARAMETERS:
# - x: input time series (usally returns);
# - Y: exogenous variable for the mean equation
# - ee: vector of innovations;
# - order: vector of p and q order;
# - type: type of garch to be estimated;
# - prob: probability function to be used;
#
# RETURNS:
#  An object of class "Garch"
#######################################################################################################################	
Garch.default = function(x, Y=NULL, order=c(alpha=1,beta=1), type=c("garch","mgarch","tgarch","egarch"), prob=c("norm","ged","t"), ...){ 
	# coerce input data to matrix and check for NA values
	if(!is.matrix(x))
		x = as.matrix(x)
	if(any(is.na(x)))
			x = x[-is.na(x),]
	# check names of order vector 	
	if(is.null(names(order)))
		names(order) = c("alpha","beta")
	# check for consistency of order parameters
	if(any(order<0) | order["alpha"]==0 | length(order)<2){
		message("Arch order must be at least = 1 and both Arch and Garch order must be positive \n No computation performed")
		return(NULL)
		}
	# get probability function and garch type
	prob = match.arg(prob)
	type = match.arg(type)
	n = NROW(x) 
	# matrix of mean regressors
	if(is.null(Y))
		Y = matrix(1, n, 1)
	# number of regressors 
	k = NCOL(Y)
	# vector of initial innovations
	ee = vector("numeric", n+max(order))
	# initial regression coefficient and residual series 
	reg = lm(x ~ 0 + Y)
	ee[-(1:max(order))] = (x - Y%*%as.matrix(reg$coef))^2
	# initial residual standard deviation
	sig0 = mean(ee[-(1:max(order))], na.rm=TRUE) 
	ee[1:max(order)] = (sig0)
	# parameters initialization 
	theta = double(sum(order)+5)
	# alpha parameter(s)
	theta[(k+2):(order["alpha"] + 2)] = in_a = rep(0.15/order["alpha"], order["alpha"])
	# beta parameter(s)
	theta[(order["alpha"]+k+2) : (order["alpha"]+order["beta"]+2)] = in_b = rep(0.45/order["beta"], order["beta"])
	# constant term
	theta[k+1] = (sig0*(1-sum(in_a)-sum(in_b)))
	# mean equation coefficient
	theta[1:k] = reg$coef
	# phi coefficient (only for egarch and tgarch)
	theta[(length(theta)-2)] = 0.1
	#
	theta[(length(theta))] = 0.01
	# coef position in theta
	pos = (2:(sum(order)+1) + k)
	# probability function and shape parameter when prob = "norms"
	if(prob == "norm"){
		r = 1
		theta[(length(theta))-1] = r
		upper = c(Inf, Inf, rep(1,sum(order)),Inf ,Inf ,Inf)
		lower = c(-Inf, rep(1e-6, sum(order+1)), -Inf, double(0), -Inf)
	} else  if(prob == "t"){
	# probability function and shape parameter when prob = "t"	
		r = 3
		theta[(length(theta))-1] = r
		upper = c(Inf, Inf, as.double(rep(1,sum(order))), Inf, Inf, Inf)
		lower = c(-Inf, rep(1e-6, sum(order+1)), -Inf, double(0), -Inf)
	} else {
	# probability function and shape parameter when prob = "ged"	
		r = 1
		theta[(length(theta))-1] = r
		upper = c(Inf, Inf, as.double(rep(1,sum(order))), -Inf, Inf, Inf)
		lower = c(-Inf, rep(1e-6, sum(order+1)), -Inf, double(0), -Inf)
	}	
	if(type == "egarch"){
		# optimise egarch parameter (hopefully)
		opt = optim(par=theta
					, fn=like.egarch
					, gr=NULL
					, ee=ee
					, x
					, Y
					, order=order
					, prob=prob
					, hessian=TRUE)
	} else if(type == "tgarch") {
		# optimise tgarch parameter (hopefully)		
		opt = optim(par = theta
					, fn=like.tgarch
					, gr=NULL
					, ee=ee
					, x
					, Y
					, order=order
					, prob=prob
					, hessian=TRUE
					, lower=lower
					, upper=upper
					, method="L-BFGS-B"
					, control = list(parscale = c(rep(1,k), signif(theta[k+1],1), rep(0.1, length(pos)), 0.1, 1, 1)
						, ndeps = c(1e-3, 1e-05, rep(1e-03, length(pos)), 1e-03, 1e-03, 1e-03))
					)
	} else if(type == "garch") {
		# optimise garch parameter (hopefully)
		opt = optim(par=theta
					, fn=like.garch
					, gr=NULL
					, ee=ee
					, x
					, Y
					, order=order
					, prob=prob
					, hessian=TRUE
					, lower=lower
					, upper=upper
					, method="L-BFGS-B"
					, control = list(parscale = c(rep(1,k), signif(theta[k+1],1), rep(0.1, length(pos)), 0.1, 1, 1)
						, ndeps = c(1e-3, 1e-06, rep(1e-03, length(pos)), 1e-03, 1e-03, 1e-03) )
					)
	} else {
		# optimise mgarch parameter (hopefully)
		opt = optim(par=theta
					, fn=like.mgarch
					, gr=NULL
					, x
					, Y
					, order=order
					, prob=prob
					, hessian=TRUE
					, lower=lower
					, upper=upper
					, method="L-BFGS-B"
					, control = list(parscale = c(rep(1,k), signif(theta[k+1],1), rep(0.1, length(pos)), 0.1, 1, 1))
					)
	}
	if(type == "garch")
		dum = 0
	else 
		dum = 1
	# extract estimated coefficient
	coef = opt$par[1:(sum(order)+1+dum) + k] 
	# extract hessian matrix and get vcov matrix
	vcov = solve(opt$hessian[(1:(sum(order)+1+dum) + k), (1:(sum(order)+1+dum) + k)])
	# get standard error (finger crossed!)
	parSD = sqrt(diag(vcov))
	# compute t-stat
	tstat = coef / parSD
	# mean equation results
	mc = opt$par[(1:k)] 
	mse = sqrt(diag(solve(opt$hessian[(k),(1:k)])))
	mts = mc / mse
	mpt = 2*(1-pnorm(abs(mts)))
	mean_coef = data.frame(Mean_Coef=mc, Mean_Se=mse, Mean_TStat=mts, Mean_PVal=mpt)
	# volatility persistence 
	if(type == "garch"){
		pers = sum(coef[-1])
	} else {
		pers = sum(coef[-1]) + coef[length(coef)]/2
	}
	# store epsilon and sigma in matrix
	fitted = matrix(0, n, 5)
	colnames(fitted) = c("Returns", "Fitted_ME", "Residuals", "Eps", "Sigma") 
	# calculate new innovations
	ee = x - Y %*% as.matrix(opt$par[1:k])
	# store residuals
	fitted[, 3] = ee
	# get fitted series with Epsilon and Sigma
	fitted[, 4:5] = .Garch.proc(pars=opt$par[-(1:k)], order=order, res=ee, type=type, r=r, prob=prob)
	# store original series
	fitted[,1] = x
	# fitted values of the mean equation
	fitted[,2] = Y %*% as.matrix(opt$par[1:k])
	# table of results
	coef.tab = data.frame(Estimates = coef, Std.Errors = parSD, T_Stat = tstat, P_Value = 2*(1-pnorm(abs(tstat))))
	rownames(coef.tab)[1] = "Omega"
	rownames(coef.tab)[2:(order["alpha"]+1)] = paste("Alpha_",1:order["alpha"],sep="")
	rownames(coef.tab)[(order["alpha"]+2):(order["beta"]+order["alpha"]+1)] = paste("Beta_",1:order["beta"],sep="")
	if(type != "garch")
		rownames(coef.tab)[NROW(coef.tab)] = "Phi"
	# create list of results (output)
	Results=list(
		Type = type,
		Order = order,
		Mean_Equation = mean_coef,
		Results = coef.tab,
		LogLik = opt$value,
		Vcov = vcov,
		Volatility_Persistence = pers,
		AIC = -2*opt$value + 2*(1+sum(order["alpha"]+order["beta"])),
		Fitted = fitted
	)
	class(Results) = "Garch"
	# clean memory
	cleanup("Results")
	Results
}
#######################################################################################################################
# FUNCTION: like.garch
#
# SUMMARY:
# This function computes the log likelihood calculation for a Garch(p ,q) model
#
# PARAMETERS:
# - Theta: vector of parameters
# - ee: vector of innovations
# - order: vector of p and q order
#
# RETURNS:
#  Vector of resutls
#######################################################################################################################	
like.garch = function(theta, ee, x, Y, order, prob=c("norm","ged","t"), r){
	# parameters
	Logger(message = "parameters", from = "like.garch", line = 2, level = 1);
	omega = theta[2] 
	alpha = theta[3:(order["alpha"] + 2)]
	beta = theta[(order["alpha"]+3) : (order["alpha"]+order["beta"]+2)]
	reg = as.matrix(theta[1])
	r = theta[length(theta)-1]
	n = length(x)		 
	# vector of initial innovations
	Logger(message = "vector of initial innovations", from = "like.garch", line = 9, level = 1);
	ee = vector("numeric", n+max(order))
	# initial regression coefficient and residual series 
	Logger(message = "initial regression coefficient and residual series ", from = "like.garch", line = 11, level = 1);
	ee[-(1:max(order))] = (x - Y%*%reg)
	# initial residual standard deviation
	Logger(message = "initial residual standard deviation", from = "like.garch", line = 13, level = 1);
	sig0 = mean(ee^2, na.rm=TRUE) 
	ee[(1:max(order))] = sqrt(sig0)
	# compute ARCH part
	Logger(message = "compute ARCH part", from = "like.garch", line = 16, level = 1);
	fit = omega + filter(ee^2, c(0,alpha), method="c", sides=1)
	fit[1] = sig0
	# compute GARCH part
	Logger(message = "compute GARCH part", from = "like.garch", line = 19, level = 1);
	if(all(!is.na(beta))){
		fit = filter(fit, beta, method="r", init = rep(0, order["beta"]))	}	
	# log-likelihood calculation   
	Logger(message = "log-likelihood calculation   ", from = "like.garch", line = 22, level = 1);
	.garch.like(X=ee[-(1:length(alpha))], S=fit[-(1:max(order))], prob=prob, r=r)  
}
#######################################################################################################################
# FUNCTION: T-Garch likelihood
#
# SUMMARY:
# This function computes the log likelihood calculation for a T-Garch(p ,q) model
#
# PARAMETERS:
# - theta: vector of parameters
# - ee: vector of innovations, residual series
# - prob: probability distribution to consider
# - order: vector of p and q order
# - r: shape parameters (GED) or degrees of freedom (T)
#
# RETURNS:
#  Likelihood value
#######################################################################################################################	
like.tgarch = function(theta, ee, x, Y, order, prob=c("norm","ged","t")){
	# parameters
	Logger(message = "parameters", from = "like.tgarch", line = 2, level = 1);
	omega = theta[2] 
	alpha = theta[3:(order["alpha"] + 2)]
	beta = theta[(order["alpha"]+3) : (order["alpha"]+order["beta"]+2)]
	phi = theta[length(theta)-2]
	reg = as.matrix(theta[1])
	r = theta[length(theta)-1]
	n = length(x)
	# vector of initial innovations
	Logger(message = "vector of initial innovations", from = "like.tgarch", line = 10, level = 1);
	ee = vector("numeric", n+max(order))
	# initial regression coefficient and residual series 
	Logger(message = "initial regression coefficient and residual series ", from = "like.tgarch", line = 12, level = 1);
	ee[-(1:max(order))] = x - Y%*%reg
	# initial residual standard deviation
	Logger(message = "initial residual standard deviation", from = "like.tgarch", line = 14, level = 1);
	sig0 = mean(ee^2, na.rm=TRUE) 
	ee[1:max(order)] = sqrt(sig0)
	# step sign
	Logger(message = "step sign", from = "like.tgarch", line = 17, level = 1);
	gb = as.numeric(ee<0);
	# asymetry factor
	Logger(message = "asymetry factor", from = "like.tgarch", line = 19, level = 1);
	asym = abs(filter(ee^2*gb, c(0,phi), method="c", sides=1))
	asym[1] = sig0
	# arch component
	Logger(message = "arch component", from = "like.tgarch", line = 22, level = 1);
	fit = omega + filter(ee^2, c(0,alpha), method="c", sides=1) + asym;
	fit[1] = sig0
	# garch component
	Logger(message = "garch component", from = "like.tgarch", line = 25, level = 1);
	fit = filter(fit, beta, method="r", init=rep(0, order["beta"]));
	# log-likelihood calculation   
	Logger(message = "log-likelihood calculation   ", from = "like.tgarch", line = 27, level = 1);
	.garch.like(X=ee[-(1:length(alpha))], S=fit[-(1:max(order))], prob=prob, r=r)  
}
#######################################################################################################################
# FUNCTION: E-Garch likelihood
#
# SUMMARY:
# This function computes the log likelihood calculation for an E-Garch(p ,q) model
#
# PARAMETERS:
# - theta: vector of parameters
# - ee: vector of innovations, residual series
# - prob: probability distribution to consider
# - order: vector of p and q order, only 1, 1 for the moment.
# - r: shape parameters (GED) or degrees of freedom (T)
#
# RETURNS:
#  Likelihood value
#######################################################################################################################	
like.egarch = function(theta, ee, x, Y, order=c(alpha=1,beta=1), prob=c("norm","ged","t")){
	# parameters
	Logger(message = "parameters", from = "like.egarch", line = 2, level = 1);
	omega = theta[2] 
	alpha = theta[3:(order["alpha"] + 2)]
	beta = theta[(order["alpha"]+3) : (order["alpha"]+order["beta"]+2)]
	phi = theta[length(theta)-2]
	reg = as.matrix(theta[1])
	r = theta[length(theta)-1]
	n = length(x)		 
	# vector of initial innovations
	Logger(message = "vector of initial innovations", from = "like.egarch", line = 10, level = 1);
	ee = vector("numeric", n+max(order))
	# initial regression coefficient and residual series 
	Logger(message = "initial regression coefficient and residual series ", from = "like.egarch", line = 12, level = 1);
	ee[-(1:max(order))] = x - Y%*%reg
	# initial residual standard deviation
	Logger(message = "initial residual standard deviation", from = "like.egarch", line = 14, level = 1);
	sig0 = mean(ee^2, na.rm=TRUE) 
	ee[1:max(order)] = sqrt(sig0)
    # initial variance
    Logger(message = "initial variance", from = "like.egarch", line = 17, level = 1);
	fit = vector("numeric", n+max(order))
	fit[1:max(order)] = log(sig0)
	exval = switch(prob,
			"norm" = (sqrt(2/pi)) ,
			"ged" = (gamma(2/r)/sqrt(gamma(1/r)*gamma(3/r))),
			"t" = (sqrt((r-2))*gamma((r-1)*0.5)/gamma(0.5)*gamma(r*0.5)) 
			);   
	i = 1 + max(order)
	# calculate likelihood values
	Logger(message = "calculate likelihood values", from = "like.egarch", line = 26, level = 1);
	while(i <= n+max(order)){
			fit[i] = omega + 
				alpha %*% ((abs(ee[(i-1):(i-max(order))])/sqrt(exp(fit[(i-1):(i-max(order))]))) - exval) + 
				beta %*% fit[(i-1):(i-max(order))] + 
				phi * (ee[i-1] / sqrt(exp(fit[i-1])))
		i = i + 1;
	};
	# log-likelihood calculation   
	Logger(message = "log-likelihood calculation   ", from = "like.egarch", line = 34, level = 1);
	.garch.like(X=ee[-(1:length(alpha))], S=exp(fit[-(1:max(order))]), prob=prob, r=r)  
}
#######################################################
like.mgarch = function(theta, x, Y, order, prob=c("norm","ged","t")){
	# parameters
	Logger(message = "parameters", from = "like.mgarch", line = 2, level = 1);
	omega = theta[2] 
	alpha = theta[3:(order["alpha"] + 2)]
	beta = theta[(order["alpha"]+3) : (order["alpha"]+order["beta"]+2)]
	reg = as.matrix(theta[1])
	delta = theta[length(theta)]
	r = theta[length(theta)-1]
	n = length(x)		 
	# vector of initial innovations
	Logger(message = "vector of initial innovations", from = "like.mgarch", line = 10, level = 1);
	ee = vector("numeric", n+max(order))
	# initial regression coefficient and residual series 
	Logger(message = "initial regression coefficient and residual series ", from = "like.mgarch", line = 12, level = 1);
	ee[-(1:max(order))] = (x - Y%*%reg)
	# initial residual standard deviation
	Logger(message = "initial residual standard deviation", from = "like.mgarch", line = 14, level = 1);
	sig0 = mean(ee^2, na.rm=TRUE) 
	ee[(1:max(order))] = sqrt(sig0)
	fit = vector("numeric", n+max(order))
	m = vector("numeric", n)
	fit[1:max(order)] = (sig0)
	i = max(order) + 1
	while(i <= n+max(order)){
		fit[i] = omega + 
			alpha %*% ee[(i-1):(i - order["alpha"])]^2 + 
			beta %*% (fit[(i-1):(i - order["beta"])])
		m[i - max(order)] = Y[i - max(order),] %*% reg + identity(fit[i]) * delta 
		ee[i] = (x[i - max(order)] - m[i])
		i = i + 1;
	 };
	# log-likelihood calculation   
	Logger(message = "log-likelihood calculation   ", from = "like.mgarch", line = 29, level = 1);
	.garch.like(X=ee[-(1:length(alpha))], S=fit[-(1:max(order))], prob=prob, r=r)  
}
#######################################################################################################################
# FUNCTION: News Impact curve (NIC)
#
# SUMMARY:
# This function create the plot of the NIC for Garch-like models
#
# PARAMETERS:
# - X: vector of innvations (x axis of the plot)
# - theta: vector of parameters
# - type: type of garch model
#
# RETURNS:
#  Plot of the News Impact Curve
#######################################################################################################################	
newsimp.default = function(x, theta, order, type=c("garch","mgarch", "egarch","tgarch"), plot=FALSE, ...){
	# coerce input data to matrix
	Logger(message = "coerce input data to matrix", from = "newsimp.default", line = 2, level = 1);
	if(!is.matrix(x) | any(is.na(x)))
		x = as.matrix(x[!is.na(x)])
	# check names of order vector 	
	Logger(message = "check names of order vector 	", from = "newsimp.default", line = 5, level = 1);
	if(is.null(names(order)))
		names(order) = c("alpha","beta")
	mtype = match.arg(type);
	# estimated coefficients 
	Logger(message = "estimated coefficients ", from = "newsimp.default", line = 9, level = 1);
	omega = theta[1]; 
	alpha = theta[ 2:(1+order["alpha"]) ]; 
	beta = theta[(order["alpha"]+2):(order["beta"]+order["alpha"]+1)];
	if(mtype == "garch"){
		phi = NULL
	} else {
		if(length(theta) == (sum(order)+1))
			stop("Parameter PHI is missing!")
		phi = theta[length(theta)-2];
	}
	# NIC garch
	Logger(message = "NIC garch", from = "newsimp.default", line = 20, level = 1);
	if(mtype == "garch"){
		nic = function(x){ 
			sig = omega / (1 - alpha - beta)
			omega + sig * beta + alpha * x^2
		}
	# NIC egarch
	Logger(message = "NIC egarch", from = "newsimp.default", line = 26, level = 1);
	} else if(mtype == "egarch") {
			nic = function(x){
				sig = exp(omega / (1 - beta )) 
				a = (alpha + phi) * x
				b = a - 2*((as.numeric(x < 0) * alpha) * x) 
				sig^(beta) * exp(omega - alpha * sqrt(2/pi)) * exp(b/sqrt(sig))
			}
	# NIC tgarch
	Logger(message = "NIC tgarch", from = "newsimp.default", line = 34, level = 1);
	} else {
		nic = function(x){
			sig = omega / (1 - alpha - beta - phi/2)
			a = (alpha + phi) * (x^2)
			b = a - ((as.numeric(x >= 0) * phi) * (x^2)) 
			omega + sig * beta + b 
		}
	}
	# get values of news impact curve
	Logger(message = "get values of news impact curve", from = "newsimp.default", line = 43, level = 1);
	vv = curve(nic, from=min(x), to=max(x), cex.axis=0.8, type="n")
	dev.off()	
	if(plot)
		cplot( vv[[2]], vv[[1]], main=paste("NIC:",mtype), ... ) 
	# return results
	Logger(message = "return results", from = "newsimp.default", line = 48, level = 1);
	cbind(Sigma = vv[[2]], Innovations = vv[[1]]) 
}
#######################################################################################################################
# FUNCTION: Test ARCH-LM
#
# SUMMARY:
# This function perform the ARCH-LM test for conditional heteroschedasticity of residuals
#
# PARAMETERS:
# - X: vector of residuals series
# - lags: maximum number of lags to use for auxiliary model
#
# RETURNS:
#  List of results with coefficient table and test statistics
#######################################################################################################################
Archlm = function(x, lags, std=FALSE, plot.acf=FALSE){
	if(class(x) == "Garch"){
		# get squared residuals from Garch object
		Logger(message = "get squared residuals from Garch object", from = "Archlm", line = 3, level = 1);
		res = x$Fitted[,3]^2
	} else {
		res = as.numeric(x)^2
	}
	if(std)
		res = scale(res, center = TRUE, scale=FALSE)
	# sample length
	Logger(message = "sample length", from = "Archlm", line = 10, level = 1);
	N = length(res[!is.na(res)])
	# lags frame
	Logger(message = "lags frame", from = "Archlm", line = 12, level = 1);
	mod = MLag(res, lag=0:lags, mode="range");
	# auxiliary regression model
	Logger(message = "auxiliary regression model", from = "Archlm", line = 14, level = 1);
	aux = lm(mod[,1] ~ mod[,-1]);
	# coefficient table
	Logger(message = "coefficient table", from = "Archlm", line = 16, level = 1);
	Coef = summary(aux)$coef;
	# R_Squared from auxiliary regression
	Logger(message = "R_Squared from auxiliary regression", from = "Archlm", line = 18, level = 1);
	r2 = summary(aux)$r.squared;
	# Test based on R-Squared
	Logger(message = "Test based on R-Squared", from = "Archlm", line = 20, level = 1);
	Test1 = N * r2;
	Pval1 = 1 - pchisq(Test1, lags);
	# Test based on F-Statistic
	Logger(message = "Test based on F-Statistic", from = "Archlm", line = 23, level = 1);
	Test2 = (r2 / lags) / ((1 - r2) / (N - lags));
	Pval2 = 1 - pf(Test2, lags, N-lags);
	# Statistics table
	Logger(message = "Statistics table", from = "Archlm", line = 26, level = 1);
	Statistics = data.frame(NxRSq = c(Test1, Pval1), F_Stat = c(Test2, Pval2))
	rownames(Statistics) = c("Statistic","P_Value")
	# list of results
	Logger(message = "list of results", from = "Archlm", line = 29, level = 1);
	Results = list(Coefficients = round(Coef, 6), Results = round(Statistics,5));
	# clean memory
	Logger(message = "clean memory", from = "Archlm", line = 31, level = 1);
	cleanup("Results")
	if(plot.acf)
		mcf(res)
	Results
}
#### Ljung_Box for serial correlations 
LjungBox = function(x, lags, plot.acf=FALSE){
	if(class(x) == "Garch"){
		# get squared residuals from Garch object
		Logger(message = "get squared residuals from Garch object", from = "LjungBox", line = 3, level = 1);
		res = x$Fitted[,3]^2
	} else {
		res = as.numeric(x)^2
	}
	# sample length
	Logger(message = "sample length", from = "LjungBox", line = 8, level = 1);
	T = length(res)
	# calculate squared acf coefficients
	Logger(message = "calculate squared acf coefficients", from = "LjungBox", line = 10, level = 1);
	r = as.numeric(acf(res, lag=lags, plot=FALSE)$acf)[-1]^2
	# Q stat
	Logger(message = "Q stat", from = "LjungBox", line = 12, level = 1);
	Q = T*(T+2) * sum(r / (T - (1:lags)))
	# P-Value
	Logger(message = "P-Value", from = "LjungBox", line = 14, level = 1);
	pv = 1 - pchisq(Q, lags)
	if(plot.acf)
		mcf(res)
	# return results
	Logger(message = "return results", from = "LjungBox", line = 18, level = 1);
	cbind(LB_stat = Q, PVal = round(pv,5))
}
newsimp.Garch = function(x, plot=TRUE, ...){
	## get parameters from object of class "Garch"	
        # residual series
	X = x$Fitted[, 3]
	# estimated parameters
	Logger(message = "estimated parameters", from = "newsimp.Garch", line = 4, level = 1);
	theta = x$Results[, 1]
	# arch - garch order
	Logger(message = "arch - garch order", from = "newsimp.Garch", line = 6, level = 1);
	order = x[[2]]
	# garch type
	Logger(message = "garch type", from = "newsimp.Garch", line = 8, level = 1);
	type = x$Type
	# call generic news impact curve function
	Logger(message = "call generic news impact curve function", from = "newsimp.Garch", line = 10, level = 1);
	newsimp.default(X, theta, order, type, plot, ...)
}
#### Statitc Prediction for Garch models ####
predict.Garch = function(object, plot=TRUE, ...){
	ff = object$Fitted
	cc = object$Results[,1]
	se = sqrt(ff[ ,4])	
	Returns = cbind(Returns_ME = ff[ ,2], Lower_SE = ff[ ,2] - 2*se, Upper_SE = ff[ ,2] + 2*se)
	Var = ff[ ,5]
	# plot predicted values
	if(plot){
		par(mfrow=c(2,1))
		cplot(Returns, main="Predicted values for Return (Mean Equation)") 
		cplot(Var, main="Predicted values for Variance") 
	}
	res = cbind(Returns, Pred_Variance = Var)
	res
}