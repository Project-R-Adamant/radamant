#############################################
VecAr = function(X, ...) UseMethod("VecAr")

VecAr.default = function(X, ar.lags = 1:2, type=c("const", "trend",
"constrend", "none"), exog=NULL, ...){
	
	#browser()

	N = NROW(X)
	C = NCOL(X)
	
	# VAR type
	type = match.arg(type)


	# include constant term
	if(type == "const"){
		Z = cbind(Const = rep(1, N), MLag(X, ar.lags, mode="selected"))
	} 

	# include trend
	else if(type == "trend"){
		Z = cbind(Trend = seq(1,N), MLag(X, ar.lags, mode="selected"))
	}

	# include trend AND constant term
	else if(type == "constrend"){
		Z = cbind(Const = rep(1, N), Trend = seq(1,N), MLag(X, ar.lags,
mode="selected"))
	}

	# neither trend nor constant term
	else if(type == "none"){
		Z = MLag(X, ar.lags, mode="selected")
	}	


	# include exogenous variables
	if(!is.null(exog)){
	
		if(!is.matrix(exog))
			exog = as.matrix(exog)
		
		if(NROW(exog) != N){
			message("Variables have different number of rows! \n No computation
performed.")
			return(NULL)	
		} else {
		
			Z = cbind(MLag(X, ar.lags, mode="selected"), exog)
		}
	
	} 

	K = NCOL(Z)
	
	
	# allocate and store data matrix
	
	temp =  matrix(NA, N, K+C)

	if(type == "const" || type == "trend"){
		temp[ ,1] = Z[ ,1]
		temp[ ,2:(C+1)] = X 
		temp[ ,(C+2):NCOL(temp)] = Z[ ,-1] 
	}
	
	else if(type == "constrend"){
		temp[ ,1:2] = Z[ ,1:2]
		temp[ ,3:(C+2)] = X 
		temp[ ,(C+3):NCOL(temp)] = Z[ ,-(1:2)] 
	}
	
	else if(type == "none"){	
		temp[ ,1:C] = X
		temp[ ,(C+1):NCOL(temp)] = Z
	}


	# estimate parameters
	est = lm(X ~ -1 + ., data = as.data.frame(Z))

	Results = vector("list", 2)
	names(Results) = c("Results", "Info_Criteria")

	Results[[1]] = est

	# Information statistics
	ee = resid(est)
	p = max(ar.lags)
	res2 = matrix(NA, 6, 1)
	rownames(res2) = c("N_Obs", "N_Var", "N_Pars", "LogLik", "AIC", "BIC")
	res2[1, ] = N
	res2[2, ] = NCOL(X)
	res2[3, ] = K 
	res2[4, ] = as.numeric(logLik(est))
	res2[5, ] = log(det(crossprod(ee) / N)) + (2*p/N)
	res2[6, ] = log(det(crossprod(ee) / N)) + (2*p/N) + (p/N) * (log(N) - 2)
	#res2[4, ] = det( crossprod(resid(est)) / (N-K) )

	Results[[2]] = res2
	
	class(Results) = "VecAr"
	attr(Results, "Data") = temp
	attr(Results, "nser") = C
	attr(Results, "nobs") = N
	attr(Results, "npar") = NROW(est[[1]])
	attr(Results, "Lag") = max(ar.lags)
	attr(Results, "Type") = type	
	
	cleanup(keep="Results")
	
	(Results)

}

################################################################
################################################################


print.VecAr = function(x, ...){

	cc = coef(x[[1]])

	K = attr(x, "nser")
	
	# Model output in formula
	row = matrix("", NCOL(cc), 1)
	rownames(row) = colnames(cc)
	for(i in 1:NCOL(cc)){
		row[i, ] = paste(paste(paste("(",round(cc[,i],5), ")", sep=""),
rownames(cc), sep = "*"), collapse = " + ")
	}

	## Fitting statistics
	stats = matrix(NA, 3, K)
	rownames(stats) = c("Sigma", "R.Sqr", "Adj.R.Sqr")
	colnames(stats) = colnames(cc)
	for(k in 1:K){
		temp = summary(x[[1]])[[k]]
		stats[ ,k] = c(temp$sigma, temp$r.squared, temp$adj.r.squared)
	}

	# F-Statistics
	ff = matrix(NA, 4, K)
	rownames(ff) = c("F-Stat", "DF_num", "DF_den","p-value")
	colnames(ff) = colnames(cc)
	for(k in 1:K){
		temp = summary(x[[1]])[[k]]
		ff[1:3, k] = c(temp$fstatistic)
		ff[4 ,k] = 1 - pf(ff[1, k], ff[2, k], ff[3, k]) 
		
	}

	cat(rep("=", 20), "\n", sep="")
	cat("Model Information:", "\n", sep="")
	print(x[[2]])
	cat(rep("=", 20), "\n", sep="")

	## print results
	for(i in 1:NCOL(cc)){
	
		cat(rep("=", 40), "\n",sep="")
		cat(paste("Model for: ", rownames(row)[i] ,sep=""), "\n")
		cat(rep("=", 40), "\n\n", sep="")
		cat("  ", paste(rownames(row)[i], row[i,], sep=" = "), "\n\n")
		cat("Fitting statistics: \n")
		cat(" ",paste(rownames(stats), "=", round(stats[,i], 5), collapse="; "),
"\n")
		cat(rep(" ", 20), "\n", sep="")
		cat("F-Statistics: \n")
		cat(paste(rownames(ff), "=", round(ff[,i], 5), collapse="; "), "\n\n\n")
	}

}

################################################################
################################################################
## ESTIMATES STRUCTURAL VAR ##

Strvar.VecAr = function(X, A=NULL, B=NULL, inter=FALSE, ...){

	if(class(X) != "VecAr"){
		cat(" \'~\' Ops! Argument X must be an objekt of class \"VecAr\" ",
"\n")
		return(NULL)
	}

	# series names
	nn = colnames(X[[1]][[1]])

	# get number of series
	K = attr(X, "nreg") 

	# get number of observations
	N = attr(X, "nobs")

	# get number of parameters
	P = attr(X, "npar")

	# get denominator df
	df = N - P 

	if(inter){
		A = B = matrix(0, K, K)
		A = edit(A)
		B = edit(B)
	}

	# check input matrixes A and B
	if(is.null(A))
		A = diag(1, K)

	if(is.null(B))
		B = diag(1, K)

	if(A == "diag")
		A = diag(NA, K)

	if(B == "diag")
		B = diag(NA, K)

	if(identical(dim(A),dim(B))){
		M = list(A,B)
		names(M) = c("A", "B")
	} else {
		cat("Matrix A and B must have the same dimensions", "\n")
		return(NULL)
	}

	# get residual series
	ee = resid(X[[1]])

	# calculate RSS
	rss = crossprod(ee) / (df)

	# number of pars to estimate
	npars = sapply(M, function(x) length(which(is.na(x))))

	if(sum(npars) == 0){
		cat("\'~\' I have nothing to optimise!", "\n")
		return(NULL)
	}

	pA = npars[1]
	pB = npars[2]

	ctyp = any(npars == 0) 

	# check identification condition
	(ctyp*2+1) * K^2 - sum(npars) < (ctyp*K^2) + K * (K-1) / 2

	# set starting values
	start = rep(0.1, sum(npars))
	
	# get pars position in A
	idp = lapply(M, function(x) which(is.na(x), TRUE))

	# log likelihood
	ll = function(start){
	
		# put starting values in matrix A
		M[[1]][idp$A] = head(start, pA)
		M[[2]][idp$B] = tail(start, pB)

		BI = solve(M[[2]])
		A = M[[1]]
		B = M[[2]]

		(K*N/2) * log(2*pi) - (N/2) * (log(det(A)^2) - log(det(B)^2) -
sum(diag(t(A) %*% crossprod(BI) %*% A %*% rss)))
	}

	# optmise parameters
	optpar = optim(start, ll, hessian=TRUE)

	Par = list(head(optpar$par, pA), tail(optpar$par, pB))

	# put estimated parameters in M
	M[[1]][idp[[1]]] = Par[[1]]
	M[[2]][idp[[2]]] = Par[[2]]

	# get SIGMA A
	SigmaA = SigmaB = matrix(0, K, K)
	
	if (!(is.null(optpar$hessian))){

		ss = sqrt(diag(solve(optpar$hessian)))
		SigmaA[idp$A] = head(ss, pA)
		SigmaB[idp$B] = tail(ss, pB)

	}	

	# get SIGMA U
	AI = solve(M$A)
	SigmaU = AI %*% crossprod(M$B) %*% t(AI)

	dimnames(A) = dimnames(SigmaA) = dimnames(SigmaU) = list(nn,nn)

	Results = list(EST_Matrix = M, SE = list(SigmaA, SigmaB), SE_U = SigmaU,
LogLik = optpar$value)

	Results
}


##################################################

fitted.VecAr = function(object, Coefs, ar.lags, ...){
	
	xx = MLag(object, ar.lags, mode="selected", na.rm=FALSE)
	fitt = (cbind(1, xx) %*% Coefs)
	resid = xx - fitt
	
	# list of results
	res = vector("list", 2)
	names(res) = c("Fitted", "Residuals")
	res[[1]] = fitt 
	res[[2]] = resid
	res  
}

###################################################
## GRANGER CAUSALITY TEST ##

GrangCas.VecAr = function(X, cause=NULL, ...){

	if(class(X) != "VecAr"){
		cat(" \'~\' Ops! Argument X must be an objekt of class \"VecAr\" ",
"\n")
		return(NULL)
	}

	# series names
	nn = colnames(coef(X[[1]]))
	
	# check cause
	if(is.null(cause)){
		cat(" \'~\' Ops! The Cause variable is not specified!", "\n")
		return(NULL)
	} else if(!all(cause %in% nn)){
		cat(" \'~\' Ops! The Cause variable does not esist", "\n")
		cat("What abbout chosing one of these? \'_\'", "\n")
		print(nn)
		return(NULL)
	} 
	
	lc = length(cause)
	
	coefs <- coef(X[[1]])
	cvar = vcov((X[[1]]),  use="na.or.complete")
	
	# results matrix
	res = matrix(NA, lc, 3)
	colnames(res) = c("Wald_Stat", "DF", "P-value")
	rownames(res) = paste(cause, "-> .", sep=" ")
	
	i = 1
	while(i <= lc){
	
		mr = which(gsub("_[[:digit:]]", "", rownames(coefs)) %in% cause[i])
		
		# create restrictions matrix
		R <- diag(length(coefs))[mr, ]
		
		# quadratic form
		qf <- R %*% cvar %*% t(R)
		# Wald statistic
		wald = t(R %*% as.vector(coefs)) %*% solve(qf) %*% (R %*%
as.vector(coefs))
		
		# results
		res[i ,1] = wald 
		res[i ,2] = NROW(R)
		res[i ,3] = 1 - pchisq(wald, NROW(R))
		
		i = i + 1
	 
	}
	
	class(res) = "GrangCas"
	
	res

}

print.GrangCas = function(x, ...){
	print.default(round(x, 5))
}


## MATRIX COEFFICIENTS FOR WOLD DECOMPOSIOTION - MA REPRESENTATION

PHI.VecAr = function(X, steps, ortho=FALSE, ...){

	Lag = attr(X, "Lag")
	type = attr(X, "Type")
	
	# get dimensions
	i = steps + 1
	K =
attr(X, "nser")

	# index to consider constant and/or trend 
	if(type == "const" || type == "trend"){
		# get array of estimated coefficients
		A = array(t(coef(X[[1]])[-1,]), dim=c(K, K, steps))
	}
	
	else if(type == "constrend"){
		# get array of estimated coefficients
		A = array(t(coef(X[[1]])[-(1:2),]), dim=c(K, K, steps))
	}
	
	else if(type == "none"){	
		# get array of estimated coefficients
		A = array(t(coef(X[[1]])), dim=c(K, K, steps))
	}

	# adjust array A for number of steps
	if(steps > Lag)
		A[, , steps:(steps-Lag+2)] = 0

	# array of PHI coefficients
	PHI = array(0, dim=c(K, K, i+i-1))
	PHI[,,i] = diag(K)

	s = i+1
	# calculate PHI recursively
	while(s <= dim(PHI)[3]){
 		temp = matrix(0,K,K)
		for(m in 1:dim(A)[3])
			temp = temp + A[,,m] %*% PHI[,,s-m] 
		PHI[,,s] = temp
		s = s+1
	}
	
	phi = PHI[,,-(1:(i-1))]
	
	if(ortho){
			
		# get number of observations
		N = attr(X, "nobs")
		# get number of parameters
		P = attr(X, "npar")
		# get denominator df
		df = N -
P - K + 1

		SE = array(NA, dim=c(K, K, steps))
		se =crossprod(resid(X[[1]])) / df
	
		chse = t(chol(se))
	
		PHI2 =array(NA, dim=dim(phi))

		j = 1
		while(j <=dim(phi)[3]){
			PHI2[,,j] = phi[,,j] %*% chse 	
			j = j +
1
		}
	
		# return array with calculated PHI orthogonal coefficients
		return(PHI2)	
	
	} else {
		# return array withcalculated PHI coefficients
		return(phi)
	}	

}


## FORECAST STANDARD ERROR ##

FSE.VecAr = function(X, steps, ...){

	K = attr(X, "nser")
	# get number of observations
	N = attr(X, "nobs")
	# get number of parameters
	P = attr(X, "npar")
	# get denominator df
	df = N - P - K + 1

	SE = array(NA, dim=c(K, K, steps))
	P = PHI.VecAr(X, steps)
	se = crossprod(resid(X[[1]])) / df
	SE[,,1] = se

	s = 2
	# calculate forecast standard error recursively
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
		FSE[i,] = sqrt(diag(SE[,,steps-i+1]))
		i = i + 1
	}
	FSE

}

## PREDICTION FOR VAR ##
predict.VecAr = function(object, steps=5, CI=0.95, viewby=c("vars","step"), ...){
	
	cc <- coef(object[[1]])
	K = attr(object, "nser")
	dm = attr(object, "Data")
	N = attr(object, "nobs")
	Lag = attr(object, "Lag") 
	type = attr(object, "Type")

	# index to consider constant and/or trend and get the last row of data
matrix
	if(type == "const" || type == "trend"){

		pred = dm[N, -1]
		ct = dm[N, 1]
	}
	
	else if(type == "constrend"){
		pred = dm[N, -(1:2)]
		ct = dm[N, (1:2)]
	}
	
	else if(type == "none"){	
		pred = dm[N, ]
		ct = NULL
	}

	# iteration to calculate steps
	if(steps == 1){
		
		res = c(ct, pred[-( length(pred):(length(pred)-K+1))]) %*% cc

	} else {

		res = c(ct, pred[-( length(pred):(length(pred)-K+1))]) %*% cc
		for(i in 2:steps){
		
			tt = c(ct, res, pred)
			temp = tt[1: (K * Lag + length(ct))] %*% cc
			res = cbind(temp, res)
		}
	}	

	res = matrix(res, steps, K, TRUE)	
	colnames(res) = colnames(cc)
	rownames(res) = paste("STEP_", steps:1, sep="")

	LV = -1 * qnorm((1 - CI) / 2) * FSE.VecAr(object, steps)
	
	viewby = match.arg(viewby)
	
	# display list of results either by variable or step
	if(viewby == "vars"){
		RES = vector("list", K)
		names(RES) = colnames(res)
		for(i in 1:K){
			RES[[i]] = cbind( res[,i] , res[,i] - LV[,i] , res[,i] + LV[,i],
LV[,i] )
			colnames(RES[[i]]) = c("Var_Fst", "Lower", "Upper", "Conf")
		}

	} else {
		RES = vector("list", steps)
		names(RES) = rownames(res)
		for(i in 1:steps){
			RES[[i]] = cbind( res[i,] , res[i,] - LV[i,] , res[i,] + LV[i,],
LV[i,]  )
			colnames(RES[[i]]) = c("Var_Fst", "Lower", "Upper", "Conf")
		}
	}

	RES
}

## IMPULSE RESPONSE FUNCTION ##
IRS.VecAr = function(X, imp, resp=NULL, steps=5, cum=TRUE, ortho=FALSE, ...){
	
	if(class(X) != "VecAr"){
		cat(" \'~\' Ops! Argument
X must be an objekt of class \"VecAr\" ", "\n")
		return(NULL)
	}

	#get series name
	yn = colnames(coef(X[[1]]))	
	
	# check cause
	if(is.null(imp)){
		cat(" \'~\' Ops! The Impulse variable(s) is
not specified!", "\n")
		return(NULL)
	} else if(!all(imp %in% yn)) {
		cat(" \'~\' Ops! The Impulse variable does not esist",
"\n")
		cat("What abbout chosing one of these? \'_\'",
"\n")
		print(yn)
		return(NULL)
	} 

	# assign default Response variables
	if(is.null(resp)){
		resp = yn[-match(imp,yn)]	
	} 

	# get coefficient PHI
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
	
			IRS[[i]] = t(P[resp, imp[i], ])
		}
	
		i = i + 1		
	}
	
	IRS

}





