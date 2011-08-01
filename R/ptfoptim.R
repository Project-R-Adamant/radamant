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
###### OPTIMUM MEAN-VARIANCE PORTFOLIO #########
## 1. ret = vector containing averge return for each asset
## 2. ptf = matrix containing one or more series of prices, one time series for each asset
## 3. mi = target return for the portfolio
## 4. SIGMA = sample covariance matrix
## 5. volatility = LOGICAL. If TRUE volatility is returned, else the variance is computed.
PtfOpt = function(ret = NULL, ptf = NULL, mi = NULL, SIGMA = NULL, volatility = TRUE, ...) {
	# get returns
	Logger(message = "get returns", from = "PtfOpt", line = 2, level = 1);
	# if target returns are not specified, calculate them from the portfolio
	Logger(message = "if target returns are not specified, calculate them from the portfolio", from = "PtfOpt", line = 3, level = 1);
	if(is.null(ret)){
		if(is.null(ptf)){
			cat("'~' Specify vector or matrix of returns!");
			return(NULL);
		}
		# number of assests
		Logger(message = "number of assests", from = "PtfOpt", line = 9, level = 1);
		N = NCOL(ptf);
		# Asset Names
		Logger(message = "Asset Names", from = "PtfOpt", line = 11, level = 1);
		asset.names = get.col.names(ptf);
		# Compute Asset Returns
		Logger(message = "Compute Asset Returns", from = "PtfOpt", line = 13, level = 1);
		ret = Ret(ptf, na.rm = TRUE);
		# Compute Average Returns
		Logger(message = "Compute Average Returns", from = "PtfOpt", line = 15, level = 1);
		M = matrix(colMeans(ret, TRUE), nrow = N, ncol = 1);
	} else {
		# number of assests
		Logger(message = "number of assests", from = "PtfOpt", line = 18, level = 1);
		N = length(ret);
		M = ret;
		# Asset Names
		Logger(message = "Asset Names", from = "PtfOpt", line = 21, level = 1);
		asset.names = colnames(ret);
	}
	# if the sample is not specified, it will be calculated from the portfolio
	Logger(message = "if the sample is not specified, it will be calculated from the portfolio", from = "PtfOpt", line = 24, level = 1);
	if(is.null(SIGMA)){
		if(!is.null(ptf)){
			SIGMA = cov(ret)
		} else {
			cat("'~' Specify either portfolio or covariance matrix!");
			return(NULL);
		}
	}
	# Number of target returns for which optimal portfolio has to be computed
	Logger(message = "Number of target returns for which optimal portfolio has to be computed", from = "PtfOpt", line = 33, level = 1);
	Nm = max(1, length(mi));
	## Get sample covariance from ptf series
	Logger(message = "Get sample covariance from ptf series", from = "PtfOpt", line = 35, level = 1);
	# inverse covariance 
	Logger(message = "inverse covariance ", from = "PtfOpt", line = 36, level = 1);
	IS = solve(SIGMA);
	# unit vector
	Logger(message = "unit vector", from = "PtfOpt", line = 38, level = 1);
	U = matrix(1, nrow = N, ncol = 1);
	# precalculate objects needed to get portfolio weights
	Logger(message = "precalculate objects needed to get portfolio weights", from = "PtfOpt", line = 40, level = 1);
	SU = IS %*% U;
	USU = as.numeric(t(U) %*% SU);
	SM = IS %*% M;
	MSM = as.numeric(t(M)%*%SM);
	USM = as.numeric(t(U)%*%SM);
	# calculate portfolio weights
	Logger(message = "calculate portfolio weights", from = "PtfOpt", line = 46, level = 1);
	if (!is.null(mi)) {
		method = "Non Minimum variance";
		# portfolio with expected returns
		Logger(message = "portfolio with expected returns", from = "PtfOpt", line = 49, level = 1);
		#w =((msm*su-usm*sm)/(msm*usu-usm*usm))+mi*((usu*sm-usm*su)/(msm*usu-usm*usm))
		Logger(message = "w =((msm*su-usm*sm)/(msm*usu-usm*usm))+mi*((usu*sm-usm*su)/(msm*usu-usm*usm))", from = "PtfOpt", line = 50, level = 1);
		w = matrix(((MSM*SU-USM*SM)/(MSM*USU-USM*USM)), nrow = N, ncol = Nm) + kronecker(matrix(mi, nrow=1),((USU*SM-USM*SU)/(MSM*USU-USM*USM)));
	} else {
		method = "Minimum variance";
		# weights for portfolio mean-variance efficient
		Logger(message = "weights for portfolio mean-variance efficient", from = "PtfOpt", line = 54, level = 1);
		w = SU / USU;
	}
	# portfolio expected overall return
	Logger(message = "portfolio expected overall return", from = "PtfOpt", line = 57, level = 1);
	eor = t(w) %*% M;
	# portfolio expected overall variance/volatility
	Logger(message = "portfolio expected overall variance/volatility", from = "PtfOpt", line = 59, level = 1);
	eov = rep(0, NCOL(w))
	for(i in 1:length(eov))
		eov[i] = ifelse(volatility, sqrt(t(w[,i]) %*% SIGMA %*% w[,i]), t(w[,i]) %*% SIGMA %*% w[,i]);
	rownames(w) = asset.names;
	# return results in list format
	Logger(message = "return results in list format", from = "PtfOpt", line = 64, level = 1);
	res = list(method = method
				, factors = list(SU = SU, SM = SM, MSM = MSM, USM = USM, USU = USU)
				, weights = w
				, eor = eor
				, eov = eov
				, type = ifelse(volatility, "Volatility", "Variance")
				);
	class(res) = "PtfOpt";
	res
}
print.PtfOpt = function(x, ...) {
	cat("\nMinimum Variance Portfolio (weigths):\n");
	print(x$weights);
	cat("\n")
	y = cbind(x$eor, x$eov);
	colnames(y) = c("Ptf Return", paste("Ptf", x$type));
	print(y);
}
## PORTFOLIO FRONTIER SIMULATION ##
## ret = vector of average prices
## SIGMA = sample covariance matrixs
## n_sim = number of simulations
## perm = vector of weights to use for the frontier
## plot = print plot
PtfFront = function(PTF=NULL
						, ret=NULL
						, SIGMA=NULL
						, mi = NULL
						, n_sim = 10
						, volatility = TRUE
						, plot = TRUE
						, main = paste("Frontier Simulation:", ifelse(is.null(mi), n_sim, length(mi)), "points")
						, xtitle = ifelse(volatility, expression(sigma), expression(sigma^2))
						, ytitle = expression(mu)
						, xlab.srt = 0
						, ytitle.srt = 0
						, type = "o"
						, legend = "Mean-Variance Frontier"
						, ...
						){
	if(is.null(ret) | is.null(SIGMA)){
		if(!is.null(PTF)){
			ret = colMeans(Ret(PTF,...), na.rm=TRUE)
			SIGMA = cov(Ret(PTF, na.rm=TRUE, ...))
		} else {
			cat("'~' ret and SIGMA are required")
			return(NULL)
		}
	}
	# if perm vector is not provided, a set of simulations is set by default with length of n_sim
	Logger(message = "if perm vector is not provided, a set of simulations is set by default with length of n_sim", from = "PtfFront", line = 11, level = 1);
	if(is.null(mi)){
		mi = seq(min(ret, na.rm = TRUE), max(ret, na.rm = TRUE), len = n_sim);
	}
	# number of simulations
	Logger(message = "number of simulations", from = "PtfFront", line = 15, level = 1);
	ns = length(mi);
	# Table of simulations
	Logger(message = "Table of simulations", from = "PtfFront", line = 17, level = 1);
	ST = matrix(0, ns, 2);
	colnames(ST) = c(ifelse(volatility, "Volatility", "Variance"), "Return");
	rownames(ST) = 1:ns;
	# compute different return and volatility for each set of weights
	Logger(message = "compute different return and volatility for each set of weights", from = "PtfFront", line = 21, level = 1);
	opt = PtfOpt(ret = ret, ptf = PTF, mi = mi, SIGMA = SIGMA, volatility = volatility);
	# return 
	Logger(message = "return ", from = "PtfFront", line = 23, level = 1);
	ST[, 2] = opt$eor;
	# volatility/variance
	Logger(message = "volatility/variance", from = "PtfFront", line = 25, level = 1);
	ST[, 1] = opt$eov;
	# plot frontier
	Logger(message = "plot frontier", from = "PtfFront", line = 27, level = 1);
	if(plot){
		cplot(100*ST[,2, drop = FALSE]
				, base = 100*(ST[,1, drop = FALSE])
				, main = main
				, xtitle = xtitle
				, ytitle = ytitle
				, xlab.srt = xlab.srt
				, ytitle.srt = ytitle.srt
				, legend = legend
				, type = type
				, xlab.suffix = "%"
				, ylab.suffix = "%"
				, ...
				)
	}
	ST
}
##### PORTFOLIO UTILITY OPTIMISATION #####
## PTF = matrix of assets
## W = vector of weights
## R = vector of Ptf returns
## SIGMA = Ptf sample covariances
## af = adversion factor
PtfUtility = function(PTF=NULL, W, R=NULL, SIGMA=NULL, af=3, plot=TRUE, ...){
	if(length(R)>2 | NCOL(PTF)>2){
		cat("'_' Hum... for the moment this function works only for two assets... but a multivariate is coming soon! \n")
		return(NULL)
	}
	if(is.null(PTF)){
		if(!is.null(R) & !is.null(SIGMA)){
			# number of assets
			Logger(message = "number of assets", from = "PtfUtility", line = 8, level = 1);
			N = length(R)
			ret = R
		} else {
			cat("Provide both vector of Ptf returns (R) and volatilities (SIGMA) \n")
			return(NULL)
		}
	} else {
		## calculate returns and volatilities from Ptf (X)
		Logger(message = "calculate returns and volatilities from Ptf (X)", from = "PtfUtility", line = 16, level = 1);
		# number of assets
		Logger(message = "number of assets", from = "PtfUtility", line = 17, level = 1);
		N = NCOL(PTF)
		if(any(is.na(PTF))){
			cat("'~' Ops! PTF contains NAs, \n")
			return(NULL)
		}
		if(N > 1){
			ret = colMeans(PTF)
			SIGMA = cov(PTF)
		} else {
			ret = c(colMeans(PTF), R)
			SIGMA = cov(PTF)
		}
	}
	if(length(W) != length(ret)){
		cat("Number of assets and number of weights does not match \n")
		return(NULL)
	}
	names(W) = names(ret) = paste("asset_", 1:length(W), sep="")
	# Quadratic utility function
	Logger(message = "Quadratic utility function", from = "PtfUtility", line = 36, level = 1);
	FU = function(a,b,af){ a - 0.5 * b * af } 
	# original returns and portfolio volatility
	Logger(message = "original returns and portfolio volatility", from = "PtfUtility", line = 38, level = 1);
	or_ret = t(W) %*% ret  #PtfRet(ret, W)
	or_var = sqrt(t(W) %*% SIGMA %*% W)   #PtfVar(W, SIGMA, vol=TRUE)
	# original utility
	Logger(message = "original utility", from = "PtfUtility", line = 41, level = 1);
	or_uti = FU(or_ret, or_var^2, af)
	# simulate a series of utility value for a given adversion factor 
	Logger(message = "simulate a series of utility value for a given adversion factor ", from = "PtfUtility", line = 43, level = 1);
	perm = seq(0,1,0.05)
	# Table of simulations
	Logger(message = "Table of simulations", from = "PtfUtility", line = 45, level = 1);
	ST = matrix(0, length(perm)+1, 3)
	colnames(ST) = c("Sim_return", "Sim_variance", "Sim_utility")
	rownames(ST) = 1:(length(perm)+1)
	# calculate return, volatility and utility for each generated weights
	Logger(message = "calculate return, volatility and utility for each generated weights", from = "PtfUtility", line = 49, level = 1);
	i = 2
	while(i <= (length(perm)+1) ){	
		ws = c(perm[i-1], 1-perm[i-1] )   
		ST[i,1] = t(ws) %*% ret #PtfRet(ret, ws) 
		ST[i,2] = (t(ws) %*% SIGMA %*% ws) #PtfVar(w=ws, SIGMA, vol=FALSE) 
		ST[i,3] = ST[i,1] - 0.5 * af * ST[i,2] 
		i = i + 1
	}
	ST[1,] = c(or_ret, or_var, or_uti)
	# Maximum utility weights
	Logger(message = "Maximum utility weights", from = "PtfUtility", line = 59, level = 1);
	un = c(1,rep(-1,N-1))
	# optimum risky weights
	Logger(message = "optimum risky weights", from = "PtfUtility", line = 61, level = 1);
	#temp = w1[2] + (ret["asset_2"] - ret["asset_1"]) / (af * t(un) %*% SIGMA %*% un)
	temp = ((un*-1) %*% ret + af*(SIGMA[,1] %*% (un))) / (af * t(un) %*% SIGMA %*% un)
	mu_w = c(1-temp, temp)
	names(mu_w) = c("Asset_1","Asset_2")
	# Maximum utility Returns
	Logger(message = "Maximum utility Returns", from = "PtfUtility", line = 66, level = 1);
	mu_r = t(mu_w) %*% ret #PtfRet(ret, mu_w)
	# Maximum utility Volatility
	Logger(message = "Maximum utility Volatility", from = "PtfUtility", line = 68, level = 1);
	mu_v = (t(mu_w) %*% SIGMA %*% mu_w)   #PtfVar(mu_w, SIGMA, vol=FALSE)
	# Maximum utility Value
	Logger(message = "Maximum utility Value", from = "PtfUtility", line = 70, level = 1);
	mu_va = FU(mu_r, mu_v, af)
	# plot results 
	Logger(message = "plot results ", from = "PtfUtility", line = 72, level = 1);
	if(plot){
		# calculate MV optimum portfolio
		Logger(message = "calculate MV optimum portfolio", from = "PtfUtility", line = 74, level = 1);
		PTF_OPT = PtfOpt(ret=ret, SIGMA=SIGMA, volatility = FALSE)
		# retrieve optimum weigths, return and variance
		Logger(message = "retrieve optimum weigths, return and variance", from = "PtfUtility", line = 76, level = 1);
		w1 = PTF_OPT$weights
		ret1 = PTF_OPT$eor
		var1 = PTF_OPT$eov
		plotdata = rbind(cbind(ST[-1,c(1, 3)], NA, NA)
						, c(NA,NA, ret1, NA)
						, c(NA,NA, NA, mu_r)
						);
		colnames(plotdata) = c("Efficient Frontier", "Utility", "Minimum Variance Portfolio", "Maximum Portfolio Utility");
		cplot(100*plotdata#ST[-1,c(1, 3)]
			  , base = 100*c(ST[-1,2], var1, mu_v)
			  , side = c(1, 2, 1, 1)
			  , cex = c(0.5, 0.5, 0.8, 0.8)
			  , main = "Utility & Returns Vs Risk"
			  , xtitle = expression(sigma^2)
			  , ytitle = expression(mu)
			  , ytitle.srt = 0
			  , ytitle2 = "Utility"
			  , type = "o"
			  , pch = c(16, 16, 11, 13)
			  , ylab.suffix = "%"
			  , ylab2.suffix = "%"
			  , xlab.suffix = "%"
			  , ...
			  )
	}
	# return data frame of results
	Logger(message = "return data frame of results", from = "PtfUtility", line = 102, level = 1);
	res = list(Max_Utility_Weight = mu_w , Results = data.frame(Max_Utility_Returns = mu_r, Max_Utility_Volatility = mu_v, Max_Utility_Value = mu_va));
	res;
}
