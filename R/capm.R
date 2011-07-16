########### CAPITAL ASSET PRICING MODEL #########
## PTF = matrix containing one series of returns for each asset
## PTF_M = vector containing returns of the benchmarl ptf
## rf = returns of risk free asset
## rfr = risk free rate
## Capm generic method
Capm = function(PTF, ...) UseMethod("Capm")
Capm.default = function(PTF, PTF_M, rf=NULL, rfr=NULL, ...){
	if(!is.matrix(PTF))
		PTF = as.matrix(PTF)
	if(!is.matrix(PTF_M))
		PTF_M = as.matrix(PTF_M)
	if(is.null(colnames(PTF))){
		nn <- paste("Ser_", 1:NCOL(PTF), sep="")
		colnames(PTF) = nn
	} else {
		nn = colnames(PTF)
	}
	# check risk free assets
	Logger(message = "check risk free assets", from = "Capm.default", line = 12, level = 1);
	if(is.null(rf))
		rf = rep(0, NROW(PTF))
	if(any(is.na(PTF)) | any(is.na(PTF_M)) | any(is.na(rf))){
		cat("'~' Ops! Data contains NAs, \n")
		return(NULL)
	}
	if(!identical(NROW(PTF), NROW(PTF_M), NROW(rf))){
		cat("'~' Ops! Data length differs for some series \n")
		return(NULL)	
	}
	# Returns and Volatility
	Logger(message = "Returns and Volatility", from = "Capm.default", line = 23, level = 1);
	RV = matrix(NA, NCOL(PTF) ,2)
	rownames(RV) = nn
	colnames(RV) = c("Returns", "Volatility")
	RV[,1] = PtfRet(PTF, glob=FALSE, ...)
	RV[,2] = PtfVar(PTF, glob=FALSE, vol=TRUE, ...)
	# Benchmarks
	Logger(message = "Benchmarks", from = "Capm.default", line = 29, level = 1);
	RV_BEN = matrix(NA, NCOL(PTF_M) ,2)
	colnames(RV_BEN) = c("Returns", "Volatility")
	RV_BEN[,1] = PtfRet(PTF_M, glob=FALSE, ...)
	RV_BEN[,2] = PtfVar(PTF_M, glob=FALSE, vol=TRUE, ...)
	# correlation matrix
	Logger(message = "correlation matrix", from = "Capm.default", line = 34, level = 1);
	cm = cor(PTF, use="na.or.complete")
	# excess return for single asset
	Logger(message = "excess return for single asset", from = "Capm.default", line = 36, level = 1);
	ex_ret = PTF - rf
	# excess return for market
	Logger(message = "excess return for market", from = "Capm.default", line = 38, level = 1);
	ex_retm = as.matrix(PTF_M - rf)
	# linear model for capm
	Logger(message = "linear model for capm", from = "Capm.default", line = 40, level = 1);
	mod = lm(ex_ret ~ ex_retm)
	# total risk 
	Logger(message = "total risk ", from = "Capm.default", line = 42, level = 1);
	tot_risk = as.matrix(sd(PTF, na.rm=TRUE))
	# sigma - non-systematic risk
	Logger(message = "sigma - non-systematic risk", from = "Capm.default", line = 44, level = 1);
	sigma = as.matrix(sd(mod$residuals, na.rm=TRUE))
	nsystr =   tot_risk - sigma
	nsystrpc = round(( nsystr / tot_risk ) * 100, 2)
	# beta - systematic risk
	Logger(message = "beta - systematic risk", from = "Capm.default", line = 48, level = 1);
	beta = as.matrix(mod$coefficients)[-1,,drop=TRUE]
	systr =  sqrt(tot_risk^2 - as.matrix( (beta)^2 * var(PTF_M, na.rm=TRUE)))
	systrpc = round((systr / tot_risk) * 100,2)
	# Risk table
	Logger(message = "Risk table", from = "Capm.default", line = 52, level = 1);
	rtable = matrix(0, NCOL(PTF), 5)
	rtable[,c(1,2,4)] = cbind(as.matrix(tot_risk), systr, nsystr)
	rtable[,c(3,5)] = cbind(systrpc, nsystrpc)
	colnames(rtable) = c("Total_Risk",
							"Systematic_Risk",
							"%Syst.",
							"NonSystematic_Risk",
							"%NonSyst.")
	rownames(rtable) = nn
	# alpha - risk free expetced return
	Logger(message = "alpha - risk free expetced return", from = "Capm.default", line = 62, level = 1);
	alpha = as.matrix(mod$coefficients)[1,,drop=TRUE] 
	#exp_ret = apply(mod$fitted.values,2,mean)
	Logger(message = "exp_ret = apply(mod$fitted.values,2,mean)", from = "Capm.default", line = 64, level = 1);
	# get risk free rate
	Logger(message = "get risk free rate", from = "Capm.default", line = 65, level = 1);
	if(is.null(rfr))
		rfr = mean(rf, na.rm=TRUE)
	# expected capm returns  
	Logger(message = "expected capm returns  ", from = "Capm.default", line = 68, level = 1);
	exp_ret = as.matrix(alpha + beta * (mean(ex_retm) - rfr))
	colnames(exp_ret) =  "Expected Return"
	res = list(Correlations = cm,
				Ret_and_Vol = RV,
				Ret_and_Vol_Bench = RV_BEN,
				Beta = beta,
				Alpha = alpha,
				Expected_Return = exp_ret, 
				Risk_Analysis = rtable)
	
	# assign class
	Logger(message = "assign class", from = "Capm.default", line = 78, level = 1);
	class(res) = "Capm"
	# clean memory
	Logger(message = "clean memory", from = "Capm.default", line = 80, level = 1);
	cleanup("res")
	res	
}
#### Portfolio expected returns ####
## ret = vector of returns
## w = vector of weights
## PTF = matrix containing price series
PtfRet = function(PTF, w=NULL, glob=TRUE, calc.ret=FALSE, ...){
	if(!is.matrix(PTF))
		PTF = as.matrix(PTF)
	if(calc.ret){
		ret = colMeans(Ret(PTF, ...), na.rm=TRUE)
	} else {
		ret = colMeans(PTF, na.rm=TRUE)
	}
	# check vector of weights	
	Logger(message = "check vector of weights	", from = "PtfRet", line = 9, level = 1);
	if(is.null(w))
		w = rep(1, length(ret))
	if(length(ret) != NROW(w)){
		message("Vector of weights of wrong length");
		return(NULL);
	};
	if(glob)
		res = t(w) %*% ret
	else
		res = ret
	res;
}
#### Portfolio standard deviation ####
## w = vector of weights
## SIGMA = sample covariance matrix
## PTF = matrix containing price series
## vol = return volatility (standard deviation)
PtfVar = function(PTF, w=NULL, glob=TRUE, vol=FALSE, calc.ret=FALSE, ...){
	if(!is.matrix(PTF))
		PTF = as.matrix(PTF)
	# compute ptf covariance matrix
	Logger(message = "compute ptf covariance matrix", from = "PtfVar", line = 4, level = 1);
	if(calc.ret){
		SIGMA = cov(Ret(PTF, na.rm=TRUE, ...))
	} else {
		SIGMA = cov(PTF, use="na.or.complete")
	}
	# check vector of weights
	Logger(message = "check vector of weights", from = "PtfVar", line = 10, level = 1);
	if(is.null(w))
		w = rep(1, NCOL(SIGMA))
	if(NCOL(SIGMA) != NROW(w)){
		message("Vector of weights of wrong length")
		return(NULL)
	};
	if(glob)
		res = ( t(w) %*% SIGMA %*% w )
	else
		res = diag(SIGMA)
	if(vol) sqrt(res) else res;
}
#### 	Portfolio Beta ####
## beta = vector of Beta coefficients
## w = vector of weights
PtfBeta = function(beta, w=NULL, glob=TRUE){
	if(class(beta) == "Capm"){
		beta = beta$Beta
	}
	# check vector of weights
	Logger(message = "check vector of weights", from = "PtfBeta", line = 5, level = 1);
	if(is.null(w))
		w = rep(1, length(beta))
	if(length(beta) != length(w)){
		warnings("Vector of weights of wrong length")
		return(NULL)
		}
	res = w %*% beta;
	res;
}
#### 4Measures ####
Sharpe = function(PTF, ...) UseMethod("Sharpe")
Jensen = function(PTF, ...) UseMethod("Jensen")
Treynor = function(PTF, ...) UseMethod("Treynor")
Appraisal = function(PTF, ...) UseMethod("Appraisal")
FourMeasures = function(PTF, ...) UseMethod("FourMeasures")
#### Sharpe ratio ####
## PTF = matrix containing price series
## rfr = risk free rate
Sharpe.default = function(PTF, rfr=0, ...){
	# ptf average values
	Logger(message = "ptf average values", from = "Sharpe.default", line = 2, level = 1);
	ptf_mi = PtfRet(PTF, ...) 
	# ptf variances	
	Logger(message = "ptf variances	", from = "Sharpe.default", line = 4, level = 1);
	ptf_sigma = PtfVar(PTF, vol=TRUE, ...)
	# calculate index	
	Logger(message = "calculate index	", from = "Sharpe.default", line = 6, level = 1);
	res = (ptf_mi - rfr) / ptf_sigma ;
	res;
}
Sharpe.Capm = function(PTF, rfr=0, ...){
	ptf_mi = PTF$Ret_and_Vol[,1]
	ptf_sigma = PTF$Ret_and_Vol[,2]
	# calculate index	
	Logger(message = "calculate index	", from = "Sharpe.Capm", line = 4, level = 1);
	res = (ptf_mi - rfr) / ptf_sigma ;
	res;
}
#### Treynor ratio ####
## PTF = matrix containing price series
## PTF_M = market (or benchmark) portfolio
## rfr = risk free rate
## rf = risk free asset
Treynor.default = function(PTF, PTF_M, rfr=0, rf=NULL, ...){
	# check risk free assets
	Logger(message = "check risk free assets", from = "Treynor.default", line = 2, level = 1);
	if(is.null(rf))
		rf = rep(0, NROW(PTF))
	# ptf average values
	Logger(message = "ptf average values", from = "Treynor.default", line = 5, level = 1);
	ptf_mi = PtfRet(PTF, ...)
	# get beta coefficients
	Logger(message = "get beta coefficients", from = "Treynor.default", line = 7, level = 1);
	beta = Capm(PTF, PTF_M=PTF_M, rf=rf, ...)$Beta ;	
	# compute index	
	Logger(message = "compute index	", from = "Treynor.default", line = 9, level = 1);
	res = (ptf_mi - rfr) / beta ;
	res;
}
Treynor.Capm = function(PTF, rfr=0, ...){
	ptf_mi = PTF$Ret_and_Vol[,1]
	beta = PTF$Beta ;	
	res = (ptf_mi - rfr) / beta ;
	res;
}
#### Jensen's Alpha ####
## PTF = matrix containing price series
## PTF_M = market (or benchmark) portfolio
## rfr = risk free rate
## rf = risk free asset
Jensen.default = function(PTF, PTF_M, rf=NULL, rfr=0, ...){
	# check risk free assets
	Logger(message = "check risk free assets", from = "Jensen.default", line = 2, level = 1);
	if(is.null(rf))
		rf = rep(0, NROW(PTF)) ;
	# ptf average values
	Logger(message = "ptf average values", from = "Jensen.default", line = 5, level = 1);
	ptf_mi = PtfRet(PTF, ...) ;
	mkt_mi = PtfRet(PTF_M, ...) ; 
	# get beta coefficients	
	Logger(message = "get beta coefficients	", from = "Jensen.default", line = 8, level = 1);
	beta = Capm(PTF, PTF_M=PTF_M, rf=rf, ...)$Beta ;	
	#compute index
	Logger(message = "compute index", from = "Jensen.default", line = 10, level = 1);
	res = ptf_mi - (rfr + beta * (mkt_mi - rfr));
	res;
}
Jensen.Capm = function(PTF, rfr=0, ...){
	ptf_mi = PTF$Ret_and_Vol[,1];
	mkt_mi = PTF$Ret_and_Vol_Bench[,1] ; 
	beta = PTF$Beta ;
	#compute index
	Logger(message = "compute index", from = "Jensen.Capm", line = 5, level = 1);
	res = ptf_mi - (rfr + beta * (mkt_mi - rfr));
	res;
}
#### Appraisal ratio ####
## PTF = matrix containing price series
## PTF_M = market (or benchmark) portfolio
## rfr = risk free rate
## rf = risk free asset
Appraisal.default = function(PTF, PTF_M, rf=NULL, rfr=0, ...){
	# check risk free assets
	Logger(message = "check risk free assets", from = "Appraisal.default", line = 2, level = 1);
	if(is.null(rf))
		rf = rep(0, NROW(PTF)) ;
	# get alpha and specific risk from CAPM
	Logger(message = "get alpha and specific risk from CAPM", from = "Appraisal.default", line = 5, level = 1);
	capm = Capm(PTF, PTF_M=PTF_M, rf=rf, ...) ;
	alpha = capm$Alpha ;
	sprisk = capm$Risk_Analysis[,2] ;
	#compute index
	Logger(message = "compute index", from = "Appraisal.default", line = 9, level = 1);
	res = alpha / sprisk
	res;
}
Appraisal.Capm = function(PTF, rfr=0, ...){
	alpha = PTF$Alpha ;
	sprisk = PTF$Risk_Analysis[,2] ;	
	#compute index
	Logger(message = "compute index", from = "Appraisal.Capm", line = 4, level = 1);
	res = alpha / sprisk
	res;
}
FourMeasures.default = function(PTF, PTF_M, rf=NULL, rfr=0, ...){
	RES = matrix(NA, 4, NCOL(PTF))
	rownames(RES) = c("Sharpe","Treynor","Jensen","Appraisal")
	if(is.null(colnames(PTF))){
		nn = paste("Ser_", 1:NCOL(PTF), sep="")
		colnames(PTF) = nn
	}
	colnames(RES) = nn
	RES[1, ] = Sharpe(PTF, rfr=rfr, ...)
	RES[2, ] = Treynor(PTF, PTF_M, rfr=rfr, rf=rf, ...)
	RES[3, ] = Jensen(PTF, PTF_M, rf=rf, rfr=rfr, ...)
	RES[4, ] = Appraisal(PTF, PTF_M, rf=rf, rfr=rfr, ...)
	RES;
}
FourMeasures.Capm = function(PTF, rfr=0, ...){
	RES = matrix(NA, 4, NCOL(PTF))
	rownames(RES) = c("Sharpe","Treynor","Jensen","Appraisal")
	colnames(RES) = rownames(PTF$Ret_and_Vol)
	RES[1, ] = Sharpe(PTF, rfr=rfr, ...)
	RES[2, ] = Treynor(PTF, rfr=rfr, ...)
	RES[3, ] = Jensen(PTF, rfr=rfr, ...)
	RES[4, ] = Appraisal(PTF, rfr=rfr, ...)
	RES;
}
## Style ananalysis ##
Styles = function(FUND, IND, W, lower=NULL, upper=NULL, ...){
	if(!identical(NROW(FUND), NROW(IND))){
		cat("'~' Ops! Data length differs for some series \n")
		return(NULL)	
	}
	if(!identical(NCOL(IND), length(W))){
		cat("'~' Ops! Vector of weights is of wrong length \n")
		return(NULL)	
	}
	if(!is.matrix(FUND))
		FUND = as.matrix(FUND)
	if(!is.matrix(IND))
		IND = as.matrix(IND)
	if(is.null(colnames(IND))){
		nn <- paste("Ser_", 1:NCOL(IND), sep="")
		colnames(IND) = nn
	} else {
		nn = colnames(IND)	
	}
	# Returns and Volatility
	Logger(message = "Returns and Volatility", from = "Styles", line = 20, level = 1);
	RV = matrix(NA, NCOL(IND) ,2)
	rownames(RV) = nn
	colnames(RV) = c("Returns", "Volatility")
	RV[,1] = PtfRet(IND, glob=FALSE, ...)
	RV[,2] = PtfVar(IND, glob=FALSE, vol=TRUE, ...)
	# Benchmarks
	Logger(message = "Benchmarks", from = "Styles", line = 26, level = 1);
	RV_F = matrix(NA, NCOL(FUND) ,2)
	colnames(RV_F) = c("Returns", "Volatility")
	RV_F[,1] = PtfRet(FUND, glob=FALSE, ...)
	RV_F[,2] = PtfVar(FUND, glob=FALSE, vol=TRUE, ...)
	# combine base fund and set of indexes
	Logger(message = "combine base fund and set of indexes", from = "Styles", line = 31, level = 1);
	X = cbind(FUND, IND)
	EV <- function(W, X){
 		# Number of assets
 		Logger(message = "Number of assets", from = "Styles", line = 34, level = 1);
		N = NCOL(X)-1;
		# Weights 
		Logger(message = "Weights ", from = "Styles", line = 36, level = 1);
 		w = as.matrix(W[1:N]);
 		err <- (X[,1] - X[,-1] %*% w)
 		track <- var(err)/var(X[, 1]) + ((sum(w) - 1)/sum(w))^2
 		track 
	}
	if(is.null(lower)) 
		lower = c(rep(0, NCOL(IND)))
	if(is.null(upper)) 
		upper = c(rep(1, NCOL(IND)))
	opt = optim(W, EV, X=X, method="L-BFGS-B", lower=lower, upper=upper, ...)
	# get optimal weights
	Logger(message = "get optimal weights", from = "Styles", line = 47, level = 1);
	optw = as.matrix(opt$par);
	rownames(optw) = nn;
	# optimum track
	Logger(message = "optimum track", from = "Styles", line = 50, level = 1);
	opt_track = (X[,1] - X[,-1] %*% optw)
	style_info = matrix(NA, NCOL(FUND), 3)
	colnames(style_info) = c("Error_mean", "Active_Var", "Style_RSQ")
	style_info[ ,1] = mean(opt_track);
	style_info[ ,2] = var(opt_track);
	style_info[ ,3] = 1 - (style_info[,2] / RV_F[,2]^2);	
	# list of results
	Logger(message = "list of results", from = "Styles", line = 57, level = 1);
	Results = list(
		Ret_and_Vol = RV,
		Ret_and_Vol_Bench = RV_F,
		Optim_weights = optw,
		Style_analysis = style_info);
		
	Results;
}
Multi.Styles = function(FUND, IND, W, n_clust=5, lower=NULL, upper=NULL, ...){
	if(!identical(NROW(FUND), NROW(IND))){
		cat("'~' Ops! Data length differs for some series \n")
		return(NULL)	
	}
	# get dimensions
	Logger(message = "get dimensions", from = "Multi.Styles", line = 6, level = 1);
	nci = NCOL(IND)
	ncf = NCOL(FUND) 
	if(!identical(nci, length(W))){
		cat("'~' Ops! Vector of weights is of wrong length \n")
		return(NULL)	
	}
	# check and convert FUND and IND to matrix
	Logger(message = "check and convert FUND and IND to matrix", from = "Multi.Styles", line = 13, level = 1);
	if(!is.matrix(FUND))
		FUND = as.matrix(FUND)
	if(!is.matrix(IND))
		IND = as.matrix(IND)
	# check and assign names to IND
	Logger(message = "check and assign names to IND", from = "Multi.Styles", line = 18, level = 1);
	if(is.null(colnames(IND))){
		nn <- paste("Ser_", 1:nci, sep="")
		colnames(IND) = nn
	} else {
		nn = colnames(IND)	
	}
	# check and assign names to FUND	
	Logger(message = "check and assign names to FUND	", from = "Multi.Styles", line = 25, level = 1);
	if(is.null(colnames(FUND))){
		nnf <- paste("Fund_", 1:ncf, sep="")
		colnames(FUND) = nnf
	} else {
		nnf = colnames(FUND)	
	}
	# combine base fund and set of indexes
	Logger(message = "combine base fund and set of indexes", from = "Multi.Styles", line = 32, level = 1);
	X = cbind(FUND, IND)
	XX = TSClust(x=X, n_clust=n_clust)
	a = summary(XX, PtfRet, glob=FALSE)
	b = summary(XX, PtfVar, glob=FALSE, vol=TRUE)
	# returns and volatility
	Logger(message = "returns and volatility", from = "Multi.Styles", line = 37, level = 1);
	RV <- array(NA, dim=c(nci, 2, n_clust))
	for(i in 1:n_clust){
		RV[1:nci, 1, i] = a[[i]][-ncf]
		RV[1:nci, 2, i] = b[[i]][-ncf]
	}
	dimnames(RV) = list(nn, c("Returns", "Volatility"), names(XX))
	# returns and volatility
	Logger(message = "returns and volatility", from = "Multi.Styles", line = 44, level = 1);
	RV_F = matrix(NA, n_clust, 2)
	colnames(RV_F) = c("Returns", "Volatility")
	#rownames(RV_F) = nnf
	Logger(message = "rownames(RV_F) = nnf", from = "Multi.Styles", line = 47, level = 1);
	for(i in 1:n_clust){
		RV_F[i, 1] = a[[i]][ncf]
		RV_F[i, 2] = b[[i]][ncf]
	}
	EV <- function(W, X){
 		# Number of assets
 		Logger(message = "Number of assets", from = "Multi.Styles", line = 53, level = 1);
		N = NCOL(X)-1;
		# Weights 
		Logger(message = "Weights ", from = "Multi.Styles", line = 55, level = 1);
 		w = as.matrix(W[1:N]);
 		err <- (X[,1] - X[,-1] %*% w)
 		#track <- var(err) + (sum(w) - 1)^2
 		Logger(message = "track <- var(err) + (sum(w) - 1)^2", from = "Multi.Styles", line = 58, level = 1);
 		track <- var(err)/var(X[, 1]) + ((sum(w) - 1)/sum(w))^2
 		track 
	}
	if(is.null(lower)) 
		lower = c(rep(0, NCOL(IND)))
	if(is.null(upper)) 
		upper = c(rep(1, NCOL(IND)))
	WW = matrix(NA, NCOL(IND), length(XX))
	i = 1
	while(i <= length(XX)){
		WW[,i ] = optim(W, EV, X=XX[[i]], method="L-BFGS-B", lower=lower, upper=upper)$par
		cat("Finished___", i, "\n")
		i = i + 1
	}
	rownames(WW) = nn;
	# optimum track
	Logger(message = "optimum track", from = "Multi.Styles", line = 74, level = 1);
	opt_track = (X[,1] - X[,-1] %*% WW)
	colnames(WW) = colnames(opt_track) = 1:n_clust
	style_info = matrix(NA, n_clust, 3)
	colnames(style_info) = c("Error_mean", "Active_Var", "Style_RSQ")
		style_info[ ,1] = apply(opt_track, 2, mean);
		style_info[ ,2] = apply(opt_track, 2, var);
		style_info[ ,3] = 1 - (style_info[,2] / RV_F[,2]^2);	
	# list of results
	Logger(message = "list of results", from = "Multi.Styles", line = 82, level = 1);
	Results = list(
		Ret_and_Vol = RV,
		Ret_and_Vol_Bench = RV_F,
		Optim_weights = WW,
		Style_analysis = style_info);
	Results;
}
