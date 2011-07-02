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
# - ee: vector of innovations;
# - order: vector of p and q order;
# - type: type of garch to be estimated;
# - prob: probability function to be used;
# - par.init: initial parameters estimation method;
# - r: shape parameters (GED) or degrees of freedom (T);
# - opt.met: optimisation method used for numerical estimation of parameters;
#
# RETURNS:
#  List of results with summary table for estimated parameters and Volatility persistence;
#######################################################################################################################	
Garch.default = function(x, Y=NULL, order=c(alpha=1,beta=1), phi=0, delta=0, type=c("garch","mgarch","tgarch","egarch"), prob=c("norm","ged","t"), ...){ 
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
	theta = rep(0, sum(order)+5)
	theta[3:(order["alpha"] + 2)] = in_a = rep(0.15/order["alpha"], order["alpha"])
	theta[(order["alpha"]+3) : (order["alpha"]+order["beta"]+2)] = in_b = rep(0.45/order["beta"], order["beta"])
	theta[2] = (sig0*(1-sum(in_a)-sum(in_b)))
	theta[1] = reg$coef
	theta[(length(theta)-2)] = phi
	theta[(length(theta))] = delta
		
	
	# probability function and shape parameter
	if(prob == "norm"){
		r = 0
		theta[(length(theta))-1] = r
		upper = c(Inf, Inf, as.double(rep(1,sum(order))),Inf,Inf)
		lower = c(-Inf, rep(1e-5,sum(order)+1), r, delta)
		
	} else  if(prob == "t"){
		r = 3
		theta[(length(theta))-1] = r
		upper = c(Inf, Inf, as.double(rep(1,sum(order))), Inf, Inf, Inf)
		lower = c(-Inf, rep(0L,sum(order)+1), -Inf, r, delta)
	} else {
		r = 1
		theta[(length(theta))-1] = r
		upper = c(Inf, Inf, as.double(rep(1,sum(order))), -Inf, 2, Inf)
		lower = c(-Inf, rep(0L,sum(order)+1), -Inf, r, delta)
	}	
	if(type == "egarch"){
	
		opt = optim(par=theta, fn=like.egarch,  gr=NULL, ee=ee, x, Y, order=order, prob=prob, hessian=TRUE)
	} else if(type == "tgarch") {
	
		opt = optim(par=theta, fn=like.tgarch,  gr=NULL, ee=ee, x, Y, order=order ,prob=prob, hessian=TRUE, lower=lower, upper=upper, method="L-BFGS-B")
	
	} else if(type == "garch") {
	
		opt = optim(par=theta, fn=like.garch,  gr=NULL, ee=ee,  x, Y, order=order, prob=prob, hessian=TRUE, lower=lower, upper=upper, method="L-BFGS-B")
		
	} else {
		
		opt = optim(par=theta, fn=like.mgarch,  gr=NULL, x, Y, order=order, prob=prob, hessian=TRUE, lower=lower, upper=upper, method="L-BFGS-B")
		
	}
	coef = opt$par[1:(sum(order)+1) + k] 
	vcov = solve(opt$hessian[1:(sum(order)+1) + k, 1:(sum(order)+1) + k])
	parSD = sqrt(diag(vcov))
	tstat = coef / parSD
	# mean equation
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
	fitted[,2] = Y %*% as.matrix(opt$par[1])
	
	# table of results
	coef.tab = data.frame(Estimates = coef, Std.Errors = parSD, T_Stat = tstat, P_Value = 2*(1-pnorm(abs(tstat))))
	
	rownames(coef.tab)[1] = "Omega"
	rownames(coef.tab)[2:(order["alpha"]+1)] = paste("Alpha_",1:order["alpha"],sep="")
	rownames(coef.tab)[(order["alpha"]+2):(order["beta"]+order["alpha"]+1)] = paste("Beta_",1:order["beta"],sep="")
	if(type != "garch")
		rownames(coef.tab)[NROW(coef.tab)] = "Phi"
	
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
	omega = theta[2] 
	alpha = theta[3:(order["alpha"] + 2)]
	beta = theta[(order["alpha"]+3) : (order["alpha"]+order["beta"]+2)]
	reg = as.matrix(theta[1])
	r = theta[length(theta)-1]
	n = length(x)		 
	# vector of initial innovations
	ee = vector("numeric", n+max(order))
	
	# initial regression coefficient and residual series 
	ee[-(1:max(order))] = (x - Y%*%reg)
	
	# initial residual standard deviation
	sig0 = mean(ee^2, na.rm=TRUE) 
	ee[(1:max(order))] = sqrt(sig0)
	
	# compute ARCH part
	fit = omega + filter(ee^2, c(0,alpha), method="c", sides=1)
	fit[1] = sig0
	# compute GARCH part
	if(all(!is.na(beta))){
		fit = filter(fit, beta, method="r", init = rep(0, order["beta"]))	}	
	# log-likelihood calculation   
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
	omega = theta[2] 
	alpha = theta[3:(order["alpha"] + 2)]
	beta = theta[(order["alpha"]+3) : (order["alpha"]+order["beta"]+2)]
	phi = theta[length(theta)-2]
	reg = as.matrix(theta[1])
	r = theta[length(theta)-1]
	
	n = length(x)
			 
	# vector of initial innovations
	ee = vector("numeric", n+max(order))
	
	# initial regression coefficient and residual series 
	ee[-(1:max(order))] = x - Y%*%reg
	
	# initial residual standard deviation
	sig0 = mean(ee^2, na.rm=TRUE) 
	ee[1:max(order)] = sqrt(sig0)
	
	# step sign
	gb = as.numeric(ee<0);
   
	# asymetry factor
	asym = abs(filter(ee^2*gb, c(0,phi), method="c", sides=1))
	asym[1] = sig0
	# arch component
	fit = omega + filter(ee^2, c(0,alpha), method="c", sides=1) + asym;
	fit[1] = sig0
	# garch component
	fit = filter(fit, beta, method="r", init=rep(0, order["beta"]));
	# log-likelihood calculation   
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
	omega = theta[2] 
	alpha = theta[3:(order["alpha"] + 2)]
	beta = theta[(order["alpha"]+3) : (order["alpha"]+order["beta"]+2)]
	phi = theta[length(theta)-2]
	reg = as.matrix(theta[1])
	r = theta[length(theta)-1]
	
	n = length(x)		 
	# vector of initial innovations
	ee = vector("numeric", n+max(order))
	
	# initial regression coefficient and residual series 
	ee[-(1:max(order))] = x - Y%*%reg
	
	# initial residual standard deviation
	sig0 = mean(ee^2, na.rm=TRUE) 
	ee[1:max(order)] = sqrt(sig0)
	
    # initial variance
	fit = vector("numeric", n+max(order))
	fit[1:max(order)] = log(sig0)
	
	exval = switch(prob,
			"norm" = (sqrt(2/pi)) ,
			"ged" = (gamma(2/r)/sqrt(gamma(1/r)*gamma(3/r))),
			"t" = (sqrt((r-2))*gamma((r-1)*0.5)/gamma(0.5)*gamma(r*0.5)) 
			);   
	i = 1 + max(order)
	# calculate likelihood values
	while(i <= n+max(order)){
		
			fit[i] = omega + 
				alpha %*% ((abs(ee[(i-1):(i-max(order))])/sqrt(exp(fit[(i-1):(i-max(order))]))) - exval) + 
				beta %*% fit[(i-1):(i-max(order))] + 
				phi * (ee[i-1] / sqrt(exp(fit[i-1])))
			
		i = i + 1;
		
	};
	# log-likelihood calculation   
	.garch.like(X=ee[-(1:length(alpha))], S=exp(fit[-(1:max(order))]), prob=prob, r=r)  
	
}

#######################################################
like.mgarch = function(theta, x, Y, order, prob=c("norm","ged","t")){
	
	# parameters
	omega = theta[2] 
	alpha = theta[3:(order["alpha"] + 2)]
	beta = theta[(order["alpha"]+3) : (order["alpha"]+order["beta"]+2)]
	reg = as.matrix(theta[1])
	delta = theta[length(theta)]
	r = theta[length(theta)-1]

	n = length(x)		 
	# vector of initial innovations
	ee = vector("numeric", n+max(order))
	
	# initial regression coefficient and residual series 
	ee[-(1:max(order))] = (x - Y%*%reg)
	
	# initial residual standard deviation
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
	if(!is.matrix(x) | any(is.na(x)))
		x = as.matrix(x[!is.na(x)])
	# check names of order vector 	
	if(is.null(names(order)))
		names(order) = c("alpha","beta")
	mtype = match.arg(type);
	
	# estimated coefficients 
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
	if(mtype == "garch"){
		
		nic = function(x){ 
			sig = omega / (1 - alpha - beta)
			omega + sig * beta + alpha * x^2
		}
		
	# NIC egarch
	} else if(mtype == "egarch") {
			
			nic = function(x){
				sig = omega / (1 - alpha - beta - phi/2)
				a = (alpha + phi) * x^2
				b = a - 2*((as.numeric(x < 0) * phi) * x^2) 
				exp(omega - alpha * sqrt(2/pi)) * sig^(2*beta) * exp(b/sig)
			}
	# NIC tgarch
	} else {
			
		nic = function(x){
			sig = omega / (1 - alpha - beta - phi/2)
			a = (alpha + phi) * (x^2)
			b = a - ((as.numeric(x >= 0) * phi) * (x^2)) 
			omega + sig * beta + b 
		}
	}
	
	# get values of news impact curve
	vv = curve(nic, from=min(x), to=max(x), cex.axis=0.8, type="n")
	dev.off()	
	if(plot)
		cplot( vv[[2]], vv[[1]], main=paste("NIC:",mtype), ... ) 
	# return results
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
		res = x$Fitted[,3]^2
	} else {
		res = as.numeric(x)^2
	}

	if(std)
		res = scale(res, center = TRUE, scale=FALSE)
	# sample length
	N = length(res[!is.na(res)])
	# lags frame
	mod = MLag(res, lag=0:lags, mode="range");
	# auxiliary regression model
	aux = lm(mod[,1] ~ mod[,-1]);
	# coefficient table
	Coef = summary(aux)$coef;
	# R_Squared from auxiliary regression
	r2 = summary(aux)$r.squared;
	
	# Test based on R-Squared
	Test1 = N * r2;
	Pval1 = 1 - pchisq(Test1, lags);
	# Test based on F-Statistic
	Test2 = (r2 / lags) / ((1 - r2) / (N - lags));
	Pval2 = 1 - pf(Test2, lags, N-lags);
	
	# Statistics table
	Statistics = data.frame(NxRSq = c(Test1, Pval1), F_Stat = c(Test2, Pval2))
	rownames(Statistics) = c("Statistic","P_Value")
	
	# list of results
	Results = list(Coefficients = round(Coef, 6), Results = round(Statistics,5));
	# clean memory
	cleanup("Results")
	
	if(plot.acf)
		mcf(res)

	Results
}


#### Ljung_Box for serial correlations 
LjungBox = function(x, lags, plot.acf=FALSE){
	
	if(class(x) == "Garch"){
		# get squared residuals from Garch object
		res = x$Fitted[,3]^2
	} else {
		res = as.numeric(x)^2
	}

	# sample length
	T = length(res)
	
	# calculate squared acf coefficients
	r = as.numeric(acf(res, lag=lags, plot=FALSE)$acf)[-1]^2
	
	# Q stat
	Q = T*(T+2) * sum(r / (T - (1:lags)))
	
	# P-Value
	pv = 1 - pchisq(Q, lags)
	
	if(plot.acf)
		mcf(res)
	
	# return results
	cbind(LB_stat = Q, PVal = round(pv,5))
	
}


newsimp.Garch = function(x, plot=TRUE, ...){
	
	## get parameters from object of class "Garch"	# original series
	X = x$Fitted[,1]
	# estimated parameters
	theta = x$Results[,1]
	# arch - garch order
	order = x[[2]]
	# garch type
	type = x$Type
	
	# call generic news impact curve function
	newsimp.default(X, theta, order, type, plot, ...)
}
#### Statitc Prediction for Garch models ####
predict.Garch = function(object, plot=TRUE, ...){
	
	ff = object$Fitted
	cc = object$Results[,1]
	se = sqrt(ff[ ,3])	
	Returns = cbind(Returns_ME = ff[ ,2], Lower_SE = ff[ ,2] - 2*se, Upper_SE = ff[ ,2] + 2*se)
	Var = ff[ ,4]
	
	# plot predicted values
	if(plot){
		par(mfrow=c(2,1))
		cplot(Returns, main="Predicted values for Return (Mean Equation)") 
		cplot(Var, main="Predicted values for Variance") 
	}
	
	res = cbind(Returns, Pred_Variance = Var)
	
	res
}
