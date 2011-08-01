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
#############################################
VecAr = function(X, ...) UseMethod("VecAr")
VecAr.default = function(X, ar.lags = 1:2, type=c("const", "trend", "constrend", "none"), regtype = "simple", exog = NULL, ...) { 
    # Get dimensions for X
    Logger(message = "Get dimensions for X", from = "VecAr.default", line = 2, level = 1);
    N = NROW(X);
    V = NCOL(X);
    if(is.null(dim(X)))
        dim(X) = c(N, V);
	# Check for exogenous variables
	Logger(message = "Check for exogenous variables", from = "VecAr.default", line = 7, level = 1);
	exog.names = NULL;
	if(!is.null(exog)) {
		exog.names = get.col.names(exog);
		# Verify dimensions of the exogenous matrix
		Logger(message = "Verify dimensions of the exogenous matrix", from = "VecAr.default", line = 11, level = 1);
		if(is.null(dim(exog))) {
			if(NROW(exog) == N) {
				dim(exog) = c(N, NCOL(exog));
			} else {
				message("Variables have different number of rows! \n No computation performed.");
				return(NULL);
			}
		}
		colnames(exog) = exog.names;
	}
	# Extract VAR type
	Logger(message = "Extract VAR type", from = "VecAr.default", line = 22, level = 1);
	type = match.arg(type);
	Const = NULL;
	Trend = NULL;
	fullnames = c();
	if(type %in% c("const", "constrend")) {
		# Include constant term
		Logger(message = "Include constant term", from = "VecAr.default", line = 28, level = 1);
		Const = rep(1, N);
		fullnames = "Const";
	} 
	if(type %in% c("trend", "constrend")) {
		# Include trend
		Logger(message = "Include trend", from = "VecAr.default", line = 33, level = 1);
		Trend = seq(1,N);
		fullnames = c(fullnames, "Trend");
	}
	# Compute lagged series
	Logger(message = "Compute lagged series", from = "VecAr.default", line = 37, level = 1);
	Xlags = MLag(X, lag = ar.lags, mode = "selected");
	# Define Regression Matrix
	Logger(message = "Define Regression Matrix", from = "VecAr.default", line = 39, level = 1);
	Z = cbind(Const, Trend, Xlags, exog);
	# Get total number of columns of the regression matrix
	Logger(message = "Get total number of columns of the regression matrix", from = "VecAr.default", line = 41, level = 1);
	K = NCOL(Z);
	# Estimate model parameters
	Logger(message = "Estimate model parameters", from = "VecAr.default", line = 43, level = 1);
	mod = mreg(Y = X, X = Z, intercept = FALSE, type = regtype, ...);
	validRows = attr(mod, "validRows");
	Results = vector("list", 2)
	names(Results) = c("Model", "Info_Criteria")
	Results[[1]] = mod;
	# Information statistics
	Logger(message = "Information statistics", from = "VecAr.default", line = 49, level = 1);
	ee = resid(mod)[validRows, drop = FALSE];
	p = max(ar.lags);
	res2 = matrix(NA, 5, 1)
	rownames(res2) = c("N_Obs", "N_Var", "N_Pars", "AIC", "BIC")
	# Number of actual observations used in the regression
	Logger(message = "Number of actual observations used in the regression", from = "VecAr.default", line = 54, level = 1);
	Nobs = sum(validRows);
	res2[1, ] = Nobs
	res2[2, ] = V
	res2[3, ] = K 
	res2[4, ] = log(det(crossprod(ee) / Nobs)) + (2*p/Nobs)
	res2[5, ] = log(det(crossprod(ee) / Nobs)) + (2*p/Nobs) + (p/Nobs) * (log(Nobs) - 2)
	Results[[2]] = res2
	# Assign class 'VecAr'
	Logger(message = "Assign class 'VecAr'", from = "VecAr.default", line = 62, level = 1);
	class(Results) = "VecAr"
	# Store data
	Logger(message = "Store data", from = "VecAr.default", line = 64, level = 1);
	Xlagfull = MLag(X, lag = max(ar.lags), mode = "auto");
	fulldata = cbind(Const, Trend, X, Xlagfull, exog);
	fullnames = c(fullnames, get.col.names(X), colnames(Xlagfull), exog.names);
	colnames(fulldata) = fullnames;
	attr(Results, "Data") = fulldata;
	attr(Results, "Xlag.names") = colnames(Xlagfull);
	# Set attributes
	Logger(message = "Set attributes", from = "VecAr.default", line = 71, level = 1);
	attr(Results, "nser") = V;
	attr(Results, "nobs") = N;
	attr(Results, "npar") = NROW(coef(mod));
	attr(Results, "exog.names") = exog.names;
	attr(Results, "Lag") = max(ar.lags);
	attr(Results, "Type") = type;
	LL = lapply(mod, function(x) round(as.numeric(logLik(x$lm)), 3))
	names(LL) = get.col.names(X);
	attr(Results, "LogLike") = LL;
	cleanup(keep="Results");
	(Results)
}
print.VecAr = function(x, ...) {
	# Print VAR model
	Logger(message = "Print VAR model", from = "print.VecAr", line = 2, level = 1);
	print(x[[1]]);
	# Show Model Matrix
	Logger(message = "Show Model Matrix", from = "print.VecAr", line = 4, level = 1);
	cat("\n===========================================\n");
	cat("===========================================\n");
	cat("Model Matrix:\n")
	show(coef(x[[1]]))
	cat("\n===========================================\n");
	cat("===========================================\n");
}
summary.VecAr = function(object, ...) {
	# Show Summary
	Logger(message = "Show Summary", from = "summary.VecAr", line = 2, level = 1);
	summary(object[[1]])
}
################################################################
## ESTIMATES STRUCTURAL VAR ##
Strvar.VecAr = function(X, A="diag", B=NULL, inter=FALSE, ...){
	if(class(X) != "VecAr"){
		cat(" \'~\' Ops! Argument X must be an object of class \"VecAr\" ", "\n")
		return(NULL)
	}
	# series names
	Logger(message = "series names", from = "Strvar.VecAr", line = 6, level = 1);
	nn = colnames(X[[1]][[1]])
	# get number of series
	Logger(message = "get number of series", from = "Strvar.VecAr", line = 8, level = 1);
	K = attr(X, "nser") 
	# get number of observations
	Logger(message = "get number of observations", from = "Strvar.VecAr", line = 10, level = 1);
	N = attr(X, "nobs")
	# get number of parameters
	Logger(message = "get number of parameters", from = "Strvar.VecAr", line = 12, level = 1);
	P = attr(X, "npar")
	# get denominator df
	Logger(message = "get denominator df", from = "Strvar.VecAr", line = 14, level = 1);
	df = N - P 
	if(inter){
		A = B = matrix(0, K, K)
		A = edit(A)
		B = edit(B)
	}
	# check input matrixes A and B
	Logger(message = "check input matrixes A and B", from = "Strvar.VecAr", line = 21, level = 1);
	if(is.null(A)){
		A = diag(1, K)
	} else if(length(A) == 1 & (A[1] == "diag")){	
		A = diag(NA, K)
	}
	if(is.null(B)){
		B = diag(1, K)
	} else if(length(B) == 1 & (B[1] == "diag")){	
		B = diag(NA, K)
	}
	# control dimensions of A and B
	Logger(message = "control dimensions of A and B", from = "Strvar.VecAr", line = 32, level = 1);
	if(identical(dim(A),dim(B))){
		M = list(A,B)
		names(M) = c("A", "B")
	} else {
		cat("Matrix A and B must have the same dimensions", "\n")
		return(NULL)
	}
	# get residual series
	Logger(message = "get residual series", from = "Strvar.VecAr", line = 40, level = 1);
	ee = residuals(X, TRUE)
	# calculate RSS
	Logger(message = "calculate RSS", from = "Strvar.VecAr", line = 42, level = 1);
	rss =  (t(ee) %*% ee) / (df) 	
	# number of pars to estimate
	Logger(message = "number of pars to estimate", from = "Strvar.VecAr", line = 44, level = 1);
	npars = sapply(M, function(x) length(which(is.na(x))))
	if(sum(npars) == 0){
		cat("\'~\' I have nothing to optimise!", "\n")
		return(NULL)
	}
	pA = npars[1]
	pB = npars[2]
	ctyp = any(npars == 0) 
	# check identification condition
	Logger(message = "check identification condition", from = "Strvar.VecAr", line = 53, level = 1);
	if((ctyp*2+1) * K^2 - sum(npars) < (ctyp*K^2) + K * (K-1) / 2)
		warning("\'~\' Problems with model identification!")
	# set starting values
	Logger(message = "set starting values", from = "Strvar.VecAr", line = 56, level = 1);
	start = rep(0.1, sum(npars))
	# get pars position in A
	Logger(message = "get pars position in A", from = "Strvar.VecAr", line = 58, level = 1);
	idp = lapply(M, function(x) which(is.na(x), TRUE))
	# log likelihood function
	Logger(message = "log likelihood function", from = "Strvar.VecAr", line = 60, level = 1);
	ll = function(start){
		# put starting values in matrix A
		Logger(message = "put starting values in matrix A", from = "Strvar.VecAr", line = 62, level = 1);
		M[[1]][idp$A] = head(start, pA)
		M[[2]][idp$B] = tail(start, pB)
		BI = solve(M[[2]])
		A = M[[1]]
		B = M[[2]]
		(K*N/2) * log(2*pi) - (N/2) * (log(det(A)^2) - log(det(B)^2) - sum(diag(t(A) %*% crossprod(BI) %*% A %*% rss)))
	}
	# optmise parameters
	Logger(message = "optmise parameters", from = "Strvar.VecAr", line = 70, level = 1);
	optpar = optim(start, ll, hessian=TRUE)
	Par = list(head(optpar$par, pA), tail(optpar$par, pB))
	# put estimated parameters in M
	Logger(message = "put estimated parameters in M", from = "Strvar.VecAr", line = 73, level = 1);
	M[[1]][idp[[1]]] = Par[[1]]
	M[[2]][idp[[2]]] = Par[[2]]
	# get SIGMA A and SIGMA B
	Logger(message = "get SIGMA A and SIGMA B", from = "Strvar.VecAr", line = 76, level = 1);
	SigmaA = SigmaB = matrix(0, K, K)
	if (!(is.null(optpar$hessian))){
		ss = sqrt(diag(solve(optpar$hessian)))
		SigmaA[idp$A] = head(ss, pA)
		SigmaB[idp$B] = tail(ss, pB)
	}	
	# get SIGMA U
	Logger(message = "get SIGMA U", from = "Strvar.VecAr", line = 83, level = 1);
	AI = solve(M$A)
	SigmaU = AI %*% crossprod(M$B) %*% t(AI)
	dimnames(A) = dimnames(SigmaA) = dimnames(SigmaU) = list(nn,nn)
	Results = list(EST_Matrix = M, SE = list(SigmaA, SigmaB), SE_U = SigmaU, LogLik = optpar$value)
	Results
}
##################################################
fitted.VecAr = function(object, ...){
	# Extract Model
	Logger(message = "Extract Model", from = "fitted.VecAr", line = 2, level = 1);
	mod = object[[1]];
	# Declare output
	Logger(message = "Declare output", from = "fitted.VecAr", line = 4, level = 1);
	res = vector("list", 2);
	names(res) = c("Fitted", "Residuals");
	# Extract Fitted and Residuals
	Logger(message = "Extract Fitted and Residuals", from = "fitted.VecAr", line = 7, level = 1);
	res[[1]] = predict(mod);
	res[[2]] = resid(mod);
	# Return result
	Logger(message = "Return result", from = "fitted.VecAr", line = 10, level = 1);
	res  
}
coef.VecAr = function(object, ...){
	# get coefficients for object "VecAr"
	if(class(object) == "VecAr")
		coef(object[[1]], ...)
	else
		NULL
}
residuals.VecAr = function(object, na.rm = FALSE, ...){
	# Ectract Residuals for object "VecAr"
	resid(object[[1]], na.rm = na.rm, ...)
}
estVar.VecAr = function(object, ...) {
	# Exctract Regression object
	Logger(message = "Exctract Regression object", from = "estVar.VecAr", line = 2, level = 1);
	mod = object[[1]];
	# Extract Residuals
	Logger(message = "Extract Residuals", from = "estVar.VecAr", line = 4, level = 1);
	err = resid(mod, na.rm = TRUE);
	# Residual degrees of freedom
	Logger(message = "Residual degrees of freedom", from = "estVar.VecAr", line = 6, level = 1);
	resdf = NROW(err) - NROW(coef(mod));
	# Extract Weights from the first model
	Logger(message = "Extract Weights from the first model", from = "estVar.VecAr", line = 8, level = 1);
	w = weights(mod[[1]], na.rm = TRUE);
	# Compute Residuals Covariance matrix
	Logger(message = "Compute Residuals Covariance matrix", from = "estVar.VecAr", line = 10, level = 1);
	if(any(is.null(w))) {
		sigma2 = crossprod(err) / resdf;
	} else {
		sigma2 = t(err) %*% diag(w) %*% err / resdf
	}
	# Return result
	Logger(message = "Return result", from = "estVar.VecAr", line = 16, level = 1);
	sigma2
}
vcov.VecAr = function(object, ...) {
	# Exctract Regression object
	Logger(message = "Exctract Regression object", from = "vcov.VecAr", line = 2, level = 1);
	mod = object[[1]];
	# Compute Residuals Covariance matrix
	Logger(message = "Compute Residuals Covariance matrix", from = "vcov.VecAr", line = 4, level = 1);
	sigma2 = estVar(object);
	# Extract Model Matrix
	Logger(message = "Extract Model Matrix", from = "vcov.VecAr", line = 6, level = 1);
	X.names = attr(mod, "X.names");
	validRows = attr(mod, "validRows");
	modMat = attr(object, "Data")[validRows, X.names, drop = FALSE];
	# Compute QR decomposition of the model matrix
	Logger(message = "Compute QR decomposition of the model matrix", from = "vcov.VecAr", line = 10, level = 1);
	Qr = qr(modMat);
	# Unscaled parameters covariance matrix
	Logger(message = "Unscaled parameters covariance matrix", from = "vcov.VecAr", line = 12, level = 1);
	cov.unscaled = chol2inv(Qr$qr);
	colnames(cov.unscaled) = X.names;
	rownames(cov.unscaled) = X.names;
	# Compute scaled covariane matrix
	Logger(message = "Compute scaled covariane matrix", from = "vcov.VecAr", line = 16, level = 1);
	res = kronecker(sigma2 , cov.unscaled, make.dimnames = TRUE);
	# Return result
	Logger(message = "Return result", from = "vcov.VecAr", line = 18, level = 1);
	res
}
###################################################
## GRANGER CAUSALITY TEST ##
GrangCas.VecAr = function(X, cause = colnames(coef(X)), digits=3, ...){
	if(class(X) != "VecAr"){
		cat(" \'~\' Ops! Argument X must be an object of class \"VecAr\" ", "\n");
		return(NULL);
	}
	# Extract series names
	Logger(message = "Extract series names", from = "GrangCas.VecAr", line = 6, level = 1);
	nn = colnames(coef(X));
	# Check cause
	Logger(message = "Check cause", from = "GrangCas.VecAr", line = 8, level = 1);
	if(is.null(cause)){
		cat(" \'~\' Ops! The Cause variable is not specified!", "\n");
		return(NULL);
	} else if(!all(cause %in% nn)){
		cat(" \'~\' Ops! The Cause variable does not esist", "\n");
		cat("What about chosing one of these? \'_\'", "\n");
		print(nn);
		return(NULL);
	} 
	lc = length(cause);
	# Extract model coefficients
	Logger(message = "Extract model coefficients", from = "GrangCas.VecAr", line = 19, level = 1);
	coefs = coef(X);
	# Extract Variance-Covariance matrix
	Logger(message = "Extract Variance-Covariance matrix", from = "GrangCas.VecAr", line = 21, level = 1);
	cvar = vcov(X);
	# Declare output matrix
	Logger(message = "Declare output matrix", from = "GrangCas.VecAr", line = 23, level = 1);
	res = matrix(NA, lc, 3);
	colnames(res) = c("Wald_Stat", "DF", "P-value");
	rownames(res) = paste(cause, "-> .", sep=" ");
	# Loop through each 'cause' entry
	Logger(message = "Loop through each 'cause' entry", from = "GrangCas.VecAr", line = 27, level = 1);
	i = 0;
	while(i < lc){
		i = i + 1;
		mr = which(sub("(_\\d)$" ,"", rownames(coefs)) %in% cause[i]);
		# Create Selection matrix
		Logger(message = "Create Selection matrix", from = "GrangCas.VecAr", line = 32, level = 2);
		R <- diag(length(coefs))[mr, ];
		# Define quadratic form
		Logger(message = "Define quadratic form", from = "GrangCas.VecAr", line = 34, level = 2);
		qf <- R %*% cvar %*% t(R);
		# Compute Wald statistic
		Logger(message = "Compute Wald statistic", from = "GrangCas.VecAr", line = 36, level = 2);
		wald = t(R %*% as.vector(coefs)) %*% solve(qf) %*% (R %*% as.vector(coefs));
		# Assign results to output
		Logger(message = "Assign results to output", from = "GrangCas.VecAr", line = 38, level = 2);
		res[i ,1] = wald;
		res[i ,2] = NROW(R);
		res[i ,3] = 1 - pchisq(wald, NROW(R));
	}
	# Return result
	Logger(message = "Return result", from = "GrangCas.VecAr", line = 43, level = 1);
	round(res, digits)
}
## MATRIX COEFFICIENTS FOR WOLD DECOMPOSIOTION - MA REPRESENTATION
PHI.VecAr = function(X, steps, ortho=FALSE, ...){
		if(class(X) != "VecAr"){
		cat(" \'~\' Ops! Argument X must be an object of class \"VecAr\" ", "\n")
		return(NULL)
	}
	Lag = attr(X, "Lag")
	type = attr(X, "Type")
	# get dimensions
	Logger(message = "get dimensions", from = "PHI.VecAr", line = 8, level = 1);
	i = steps + 1
	K = attr(X, "nser")
	# index to consider constant and/or trend 
	Logger(message = "index to consider constant and/or trend ", from = "PHI.VecAr", line = 11, level = 1);
	if(type == "const" || type == "trend"){
		# get array of estimated coefficients
		Logger(message = "get array of estimated coefficients", from = "PHI.VecAr", line = 13, level = 1);
		A = array(t(coef(X[[1]])[-1,]), dim=c(K, K, steps))
	}
	else if(type == "constrend"){
		# get array of estimated coefficients
		Logger(message = "get array of estimated coefficients", from = "PHI.VecAr", line = 17, level = 1);
		A = array(t(coef(X[[1]])[-(1:2),]), dim=c(K, K, steps))
	}
	else if(type == "none"){	
		# get array of estimated coefficients
		Logger(message = "get array of estimated coefficients", from = "PHI.VecAr", line = 21, level = 1);
		A = array(t(coef(X[[1]])), dim=c(K, K, steps))
	}
	# adjust array A for number of steps
	Logger(message = "adjust array A for number of steps", from = "PHI.VecAr", line = 24, level = 1);
	if(steps > Lag){
		if(Lag == 1)
			A[, , -1] = 0
		else
			A[, , -(1:Lag)] = 0
	}
	# array of PHI coefficients
	Logger(message = "array of PHI coefficients", from = "PHI.VecAr", line = 31, level = 1);
	PHI = array(0, dim=c(K, K, i+i-1))
	PHI[,,i] = diag(K)
	s = i+1
	# calculate PHI recursively
	Logger(message = "calculate PHI recursively", from = "PHI.VecAr", line = 35, level = 1);
	while(s <= dim(PHI)[3]){
 		temp = matrix(0,K,K)
		for(m in 1:dim(A)[3])
			temp = temp + A[,,m] %*% PHI[,,s-m] 
		PHI[,,s] = temp
		s = s+1
	}
	phi = PHI[,,-(1:(i-1)), drop = FALSE]
	if(ortho){
		# get number of observations
		Logger(message = "get number of observations", from = "PHI.VecAr", line = 45, level = 1);
		N = attr(X, "nobs")
		# get number of parameters
		Logger(message = "get number of parameters", from = "PHI.VecAr", line = 47, level = 1);
		P = attr(X, "npar")
		SE = array(NA, dim=c(K, K, steps))
		se = estVar(X);
		chse = t(chol(se))
		PHI2 =array(NA, dim=dim(phi))
		j = 1
		while(j <=dim(phi)[3]){
			PHI2[,,j] = phi[,,j] %*% chse 	
			j = j + 1
		}
		# return array with calculated PHI orthogonal coefficients
		Logger(message = "return array with calculated PHI orthogonal coefficients", from = "PHI.VecAr", line = 58, level = 1);
		return(PHI2)	
	} else {
		# return array withcalculated PHI coefficients
		Logger(message = "return array withcalculated PHI coefficients", from = "PHI.VecAr", line = 61, level = 1);
		return(phi)
	}	
}
## FORECAST STANDARD ERROR ##
FSE.VecAr = function(X, steps, ...){
	if(class(X) != "VecAr"){
		cat(" \'~\' Ops! Argument X must be an object of class \"VecAr\" ", "\n")
		return(NULL)
	}
	K = attr(X, "nser")
	# get number of observations
	Logger(message = "get number of observations", from = "FSE.VecAr", line = 7, level = 1);
	N = attr(X, "nobs")
	# get number of parameters
	Logger(message = "get number of parameters", from = "FSE.VecAr", line = 9, level = 1);
	Np = attr(X, "npar")
	SE = array(NA, dim=c(K, K, steps))
	P = PHI.VecAr(X, steps)
	# Compute Residuals Covariance matrix
	Logger(message = "Compute Residuals Covariance matrix", from = "FSE.VecAr", line = 13, level = 1);
	se = estVar(X);
	SE[,,1] = se
	s = 2
	# calculate forecast standard error recursively
	Logger(message = "calculate forecast standard error recursively", from = "FSE.VecAr", line = 17, level = 1);
	while(s <= steps){
 		temp = matrix(0,K,K)
		for(m in 1:(s-1))
			temp = temp + P[, , m] %*% se %*% t(P[, , m])
		SE[,,s] = temp + SE[,,1]
		s = s+1
	}
	FSE = matrix(NA, steps, K)
	i = 1 
	while(i <= steps){
		#FSE[i,] = sqrt(diag(SE[,,steps-i+1]))
		Logger(message = "FSE[i,] = sqrt(diag(SE[,,steps-i+1]))", from = "FSE.VecAr", line = 28, level = 2);
		FSE[i,] = sqrt(diag(as.matrix(SE[, , i])))
		i = i + 1
	}
	FSE
}
## PREDICTION FOR VAR ##
predict.VecAr = function(object
						, exog = NULL
						, steps = 5
						, ci = 0.95
						, simulate = FALSE
						, sd.sim = 1
						, aggregate = TRUE
						, scenarios = 1
						, plot = TRUE
						, ...
						) {
	# Extract model coefficients
	Logger(message = "Extract model coefficients", from = "predict.VecAr", line = 2, level = 1);
	beta = coef(object);
	# Extract the data model
	Logger(message = "Extract the data model", from = "predict.VecAr", line = 4, level = 1);
	dm = attr(object, "Data");		
	# Extract number of series
	Logger(message = "Extract number of series", from = "predict.VecAr", line = 6, level = 1);
	V = attr(object, "nser");
	# Extract the number of observations
	Logger(message = "Extract the number of observations", from = "predict.VecAr", line = 8, level = 1);
	N = attr(object, "nobs");
	# Extract the number of parameters
	Logger(message = "Extract the number of parameters", from = "predict.VecAr", line = 10, level = 1);
	P = attr(object, "npar");
	# Extract the number of exogenous variables
	Logger(message = "Extract the number of exogenous variables", from = "predict.VecAr", line = 12, level = 1);
	exog.names = attr(object, "exog.names");
	Nexog = length(exog.names);
	# Extract the order of the model
	Logger(message = "Extract the order of the model", from = "predict.VecAr", line = 15, level = 1);
	Lag = attr(object, "Lag");
	type = attr(object, "Type");
	# Extract regressors names (excluding intercept and trend)
	Logger(message = "Extract regressors names (excluding intercept and trend)", from = "predict.VecAr", line = 18, level = 1);
	X.names = attr(object[[1]], "X.names");
	ct.idx = which(X.names %in% c("Const", "Trend"));
	if(length(ct.idx) > 0)
		X.names = X.names[-ct.idx];
	# Extract names of all state variables
	Logger(message = "Extract names of all state variables", from = "predict.VecAr", line = 23, level = 1);
	state.names = attr(object, "Xlag.names");
	# Compute state selector index
	Logger(message = "Compute state selector index", from = "predict.VecAr", line = 25, level = 1);
	state.idx = which(state.names %in% X.names);
	# Check for exogenous variables
	Logger(message = "Check for exogenous variables", from = "predict.VecAr", line = 27, level = 1);
	if(Nexog > 0 && (is.null(exog) || Nexog != NCOL(exog) || NROW(exog) < steps)) {
		# No exogenous variables were provided. Stopping execution.
		Logger(message = "No exogenous variables were provided. Stopping execution.", from = "predict.VecAr", line = 29, level = 1);
		warning("The input model is based on exogenous variables which were not provided in the input argument 'exog'.");
		return(NULL);
	}
	# Initial state for the prediction (last record of the data matrix)
	Logger(message = "Initial state for the prediction (last record of the data matrix)", from = "predict.VecAr", line = 33, level = 1);
	Const = NULL;
	Trend = NULL;
	Cname = NULL;
	# Extract predictor terms
	Logger(message = "Extract predictor terms", from = "predict.VecAr", line = 37, level = 1);
	state = dm[N, state.names, drop = FALSE];
	if(type == "const" || type == "trend"){
		if(type == "const") {
			# Extract constant term
			Logger(message = "Extract constant term", from = "predict.VecAr", line = 41, level = 1);
			Const = dm[N, 1];
			Cname = "Const";
		} else {
			# Extract trend term
			Logger(message = "Extract trend term", from = "predict.VecAr", line = 45, level = 1);
			Trend = dm[N, 1] + seq(1, steps);
		}
	}
	else if(type == "constrend"){
		# Extract constant and trend terms
		Logger(message = "Extract constant and trend terms", from = "predict.VecAr", line = 50, level = 1);
		Const = dm[N, 1];
		Trend = dm[N, 2] + seq(1, steps);
		Cname = "Const";
	}
	# Set number of regressors
	Logger(message = "Set number of regressors", from = "predict.VecAr", line = 55, level = 1);
	Lstate = length(state);
	# Compute Residuals Covariance matrix
	Logger(message = "Compute Residuals Covariance matrix", from = "predict.VecAr", line = 57, level = 1);
	sigma2 = diag(estVar(object));
	# Declare forecast matrix
	Logger(message = "Declare forecast matrix", from = "predict.VecAr", line = 59, level = 1);
	res = array(NA, dim = c(steps, V, scenarios));
	dimnames(res) = list(paste("Step", seq(1, steps), sep = "_")
						, colnames(beta)
						, paste("Scenario", seq(1, scenarios), sep = "_")
						);
	# Save starting state
	Logger(message = "Save starting state", from = "predict.VecAr", line = 65, level = 1);
	state.base = state;
	# Loop through all scenarios
	Logger(message = "Loop through all scenarios", from = "predict.VecAr", line = 67, level = 1);
	i = 0;
	while(i < scenarios) {
		i = i + 1;
		# Reset the state for the new scenario
		Logger(message = "Reset the state for the new scenario", from = "predict.VecAr", line = 71, level = 2);
		state = state.base;
		# Run forecast scenario
		Logger(message = "Run forecast scenario", from = "predict.VecAr", line = 73, level = 2);
		n = 0;
		while(n < steps) {
			n = n + 1;
			# Determine if a random innovation term should be added to the forecast
			Logger(message = "Determine if a random innovation term should be added to the forecast", from = "predict.VecAr", line = 77, level = 3);
			if(simulate) {
				# Generate random data for the innovation terms
				Logger(message = "Generate random data for the innovation terms", from = "predict.VecAr", line = 79, level = 3);
				eta = sigma2 * sd.sim * rnorm(V);
			} else {
				# No innovation term required
				Logger(message = "No innovation term required", from = "predict.VecAr", line = 82, level = 3);
				eta = rep(0, V);
			}
			# Compute one step forecast
			Logger(message = "Compute one step forecast", from = "predict.VecAr", line = 85, level = 3);
			res[n, , i] = c(Const, Trend[n], state[state.idx], exog[n, exog.names]) %*% beta + eta;
			# Update state for next step
			Logger(message = "Update state for next step", from = "predict.VecAr", line = 87, level = 3);
			state = c(res[n, , i], state[-c(Lstate-seq(0, V-1))]);
			names(state) = state.names;
		}
	}
	# Compute VAR forecast standard errors
	Logger(message = "Compute VAR forecast standard errors", from = "predict.VecAr", line = 92, level = 1);
	fcast.se = -1 * qnorm((1 - ci) / 2) * FSE.VecAr(object, steps);
	colnames(fcast.se) = paste(colnames(beta), "se", sep = ".");
	rownames(fcast.se) = rownames(res);
	if(aggregate) {
		RES = matrix(NA, nrow = steps, ncol = V*3+length(Cname));
		colnames(RES) = c(colnames(res)
						, paste(colnames(res), "lwr", sep = ".")
						, paste(colnames(res), "upr", sep = ".")
						, Cname
						);
		rownames(RES) = rownames(res);
		if(scenarios == 1) {
			RES[, ] = cbind(res[, , 1], res[, , 1]-fcast.se, res[, , 1]+fcast.se, Const);
		} else {
			# Compute Average
			Logger(message = "Compute Average", from = "predict.VecAr", line = 107, level = 1);
			RES[, 1:V] = apply(res, 1:2, mean);
			# Lower quantile
			Logger(message = "Lower quantile", from = "predict.VecAr", line = 109, level = 1);
			RES[, V + 1:V] = apply(res, 1:2, quantile, probs = 1-ci);
			# Upper quantile
			Logger(message = "Upper quantile", from = "predict.VecAr", line = 111, level = 1);
			RES[, 2*V + 1:V] = apply(res, 1:2, quantile, probs = ci);
			if(length(Cname) > 0) {
				# Add constant term
				Logger(message = "Add constant term", from = "predict.VecAr", line = 114, level = 1);
				RES[, 3*V + 1] = Const;
			}
		}
	} else {
		if(length(Cname) > 0) {	
			# Add constant term
			Logger(message = "Add constant term", from = "predict.VecAr", line = 120, level = 1);
			RES = array(NA, dim = dim(res) + c(0, 1, 0));
			RES[, 1:NCOL(res), ] = res;
			RES[, NCOL(res)+1, ] = Const;
			dimnames(RES) = list(rownames(res), c(colnames(res), Cname), dimnames(res)[[3]]);
		} else {
			RES = res;
		}
	}
	# Assign class
	Logger(message = "Assign class", from = "predict.VecAr", line = 129, level = 1);
	attr(RES, "class") = "predVecAr";
	# Assign attributes
	Logger(message = "Assign attributes", from = "predict.VecAr", line = 131, level = 1);
	attr(RES, "snames") = colnames(beta);
	attr(RES, "ci") = ci;
	attr(RES, "aggregate") = aggregate;
	attr(RES, "formula") = formula(object[[1]]);
	attr(RES, "fcast.se") = fcast.se;
	attr(RES, "fitted") = fitted(object);
	if(plot)
		plot(RES, ...);
	# Return result
	Logger(message = "Return result", from = "predict.VecAr", line = 140, level = 1);
	RES
}
print.predVecAr = function(x, ...) {
	aggregate = attr(x, "aggregate");
	# Show forecast
	Logger(message = "Show forecast", from = "print.predVecAr", line = 3, level = 1);
	if(aggregate) {
		show(x[, , drop = FALSE]);
	} else {
		show(x[, , , drop = FALSE]);
	}
	cat("\nForecast Standard Error:\n");
	show(attr(x, "fcast.se"));
}
plot.predVecAr = function(x
						, main = "VAR Forecast"
						, xlabels = NULL
						, legend = NULL
						, theme.params = getCurrentTheme()
						, shaded = FALSE
						, ...
						) {
	# Extract attributes
	Logger(message = "Extract attributes", from = "plot.predVecAr", line = 2, level = 1);
	aggregate = attr(x, "aggregate");
	fit = attr(x, "fitted");
	fcast.se = attr(x, "fcast.se");
	formulas = attr(x, "formula")
	snames = attr(x, "snames");
	ci.pct = sprintf("%.5g%%", attr(x, "ci")*100);
	# Set default plotting parameters
	Logger(message = "Set default plotting parameters", from = "plot.predVecAr", line = 9, level = 1);
	theme.params[["col"]] = theme.params[["col"]][c(1, 2, 2)];
	# Extract dimensions
	Logger(message = "Extract dimensions", from = "plot.predVecAr", line = 11, level = 1);
	steps = NROW(x);
	V = length(snames);
	if(aggregate) {
		scenarios = 1;
	} else {
		scenarios = dim(x)[3];
	}
	# Process legend
	Logger(message = "Process legend", from = "plot.predVecAr", line = 19, level = 1);
	if(is.null(legend))
		legend = formulas;
    # Get plot layout
    Logger(message = "Get plot layout", from = "plot.predVecAr", line = 22, level = 1);
    plot.layout = get.plot.layout(N = V, theme.params = theme.params, overrides = list(...));
    plots.per.window = prod(plot.layout);
	n = 0;
	while(n < scenarios) {
		n = n + 1;
		v = 0;
		while(v < V) {
			v = v + 1;
			# Precess x-labels
			Logger(message = "Precess x-labels", from = "plot.predVecAr", line = 31, level = 3);
			if(is.null(xlabels))
				xlabels = c(rownames(fit[["Fitted"]][[snames[v]]]), rownames(x));			
			if(aggregate) {
				# Collect all data (fitted + aggregate predition)
				Logger(message = "Collect all data (fitted + aggregate predition)", from = "plot.predVecAr", line = 35, level = 3);
				fulldata = rbind(fit[["Fitted"]][[snames[v]]]
								, x[, paste(snames[v], c("", ".lwr", ".upr"), sep = "")]
								);
			} else {
				# Collect all data (fitted + scenario predition)
				Logger(message = "Collect all data (fitted + scenario predition)", from = "plot.predVecAr", line = 40, level = 3);
				fulldata = rbind(fit[["Fitted"]][[snames[v]]]
								, cbind(x[, snames[v], n]
										, x[, snames[v], n] - fcast.se[, paste(snames[v], "se", sep = ".")]
										, x[, snames[v], n] + fcast.se[, paste(snames[v], "se", sep = ".")]
										)
								);			
			}
			if( ((v %% plots.per.window) ==1) || plots.per.window == 1 ) {
				dev.new();
				# Set the number  of plottable areas in the window
				Logger(message = "Set the number  of plottable areas in the window", from = "plot.predVecAr", line = 50, level = 3);
				par(mfrow = plot.layout);
			}
			# Plot Results
			Logger(message = "Plot Results", from = "plot.predVecAr", line = 53, level = 3);
			cplot(fulldata
				, main = main
				, legend = c(legend[v], paste("C.I.", ci.pct))
				, theme.params = theme.params
				, xlabels = xlabels
				, ...
				);
			if(shaded) {
				# Add shaded area
				Logger(message = "Add shaded area", from = "plot.predVecAr", line = 62, level = 3);
				shade.plot(from = fulldata[, 2], to = fulldata[, 3]
							, theme.params = theme.params
							, ...);
				# Replot data points on top of the shade
				Logger(message = "Replot data points on top of the shade", from = "plot.predVecAr", line = 66, level = 3);
				cplot(fulldata
					, append = TRUE
					, main = ""
					, theme.params = theme.params
					, xlabels = xlabels
					, show.legend = FALSE
					, show.xlabels = FALSE
					, ...
					);
			}
		}
	}
}
## IMPULSE RESPONSE FUNCTION ##
IRS.VecAr = function(X, imp, resp=NULL, steps=5, cum=TRUE, ortho=FALSE, ...){
	if(class(X) != "VecAr"){
		cat(" \'~\' Ops! Argument X must be an object of class \"VecAr\" ", "\n")
		return(NULL)
	}
	#get series name
	Logger(message = "get series name", from = "IRS.VecAr", line = 6, level = 1);
	yn = colnames(coef(X[[1]]))	
	# check cause
	Logger(message = "check cause", from = "IRS.VecAr", line = 8, level = 1);
	if(is.null(imp)){
		cat(" \'~\' Ops! The Impulse variable(s) is not specified!", "\n")
		return(NULL)
	} else if(!all(imp %in% yn)) {
		cat(" \'~\' Ops! The Impulse variable does not esist", "\n")
		cat("What about chosing one of these? \'_\'", "\n")
		print(yn)
		return(NULL)
	} 
	# assign default Response variables
	Logger(message = "assign default Response variables", from = "IRS.VecAr", line = 18, level = 1);
	if(is.null(resp)){
		resp = yn[-match(imp,yn)]	
	} 
	# get coefficient PHI
	Logger(message = "get coefficient PHI", from = "IRS.VecAr", line = 22, level = 1);
	P = PHI.VecAr(X, steps, ortho)
	dimnames(P) = list(yn, yn, NULL)
	li = length(imp)
	IRS = vector("list", li)
	i = 1
	while(i <= li){
		if(cum){
			if(length(resp) == 1){
				IRS[[i]] = cumsum(t(P[resp, imp[i], ]))
			} else {
				IRS[[i]] = apply(t(P[resp, imp[i], ]), 2, cumsum)
			}
		} else {
			IRS[[i]] = t(P[resp, imp[i],])
		}	
		i = i + 1		
	}
	IRS
}
