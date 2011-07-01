###### OPTIMUM MEAN-VARIANCE PORTFOLIO #########
## 1. ret = vector containing averge return for each asset
## 2. ptf = matrix containing one or more series of returns, one time series for each asset
## 3. mi = target return for the portfolio
## 4. SIGMA = sample covariance matrix
## 5. volatility = LOGICAL. If TRUE volatility is returned, else the variance is computed.

PtfOpt = function(ret = NULL, ptf = NULL, mi = NULL, SIGMA = NULL, volatility = TRUE, ...) {

	# get returns
	# if target returns are not specified, calculate them from the portfolio
	if(is.null(ret)){
		
		if(is.null(ptf)){
			cat("'~' Specify vector or matrix of returns!");
			return(NULL);
		}
		
		# number of assests
		N = NCOL(ptf);
		
		# Asset Names
		asset.names = get.col.names(ptf);
	
		# Compute Asset Returns
		ret = Ret(ptf, na.rm = TRUE);
		
		# Compute Average Returns
		M = matrix(colMeans(ret,TRUE), nrow = N, ncol = 1);
	
	} else {
		
		# number of assests
		N = length(ret);
		M = ret;
		# Asset Names
		asset.names = names(ptf);
		
	}

	# if the sample is not specified, it will be calculated from the portfolio
	if(is.null(SIGMA)){
		if(!is.null(ptf)){
			SIGMA = cov(ret)
		} else {
			cat("'~' Specify either portfolio or covariance matrix!");
			return(NULL);
		}
	}
	
	# Number of target returns for which optimal portfolio has to be computed
	Nm = max(1, length(mi));
	
	# Get sample covariance from ptf series
	
	# inverse covariance 
	IS = solve(SIGMA);
	# unit vector
	U = matrix(1, nrow = N, ncol = 1);
	
	# precalculate objects needed to get portfolio weights
	SU = IS %*% U;
	USU = as.numeric(t(U) %*% SU);
	SM = IS %*% M;
	MSM = as.numeric(t(M)%*%SM);
	USM = as.numeric(t(U)%*%SM);
	
	# calculate portfolio weights
	if (!is.null(mi)) {

		method = "Non Minimum variance";
		
		# portfolio with expected returns
		#w =((msm*su-usm*sm)/(msm*usu-usm*usm))+mi*((usu*sm-usm*su)/(msm*usu-usm*usm))
		w = matrix(((MSM*SU-USM*SM)/(MSM*USU-USM*USM)), nrow = N, ncol = Nm) + kronecker(matrix(mi, nrow=1),((USU*SM-USM*SU)/(MSM*USU-USM*USM)));
	} else {
		method = "Minimum variance";
		# weights for portfolio mean-variance efficient
		w = SU / USU;
	}
	
	# portfolio expected overall return
	eor = PtfRet(M ,w);
	# portfolio expected overall variance/volatility
	eov = PtfVar(w, SIGMA, vol = volatility);

	rownames(w) = asset.names;
	# return results in list format
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
	if(is.null(mi)){
		mi = seq(min(ret, na.rm = TRUE), max(ret, na.rm = TRUE), len = n_sim);
	}
	
	# number of simulations
	ns = length(mi);
	
	# Table of simulations
	ST = matrix(0, ns, 2);
	colnames(ST) = c(ifelse(volatility, "Volatility", "Variance"), "Return");
	rownames(ST) = 1:ns;
	
	# compute different return and volatility for each set of weights
	opt = PtfOpt(ret = ret, ptf = PTF, mi = mi, SIGMA = SIGMA, volatility = volatility);
	# return 
	ST[, 2] = opt$eor;
	# volatility/variance
	ST[, 1] = opt$eov;
	
	# plot frontier
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


PtfUtility = function(PTF=NULL, W, R=NULL, SIGMA=NULL, af=3, plot=TRUE, plot.mv=FALSE, plot.mu=FALSE, ...){
	
	if(length(R)>2 | NCOL(PTF)>2){
		cat("'_' Hum... for the moment this function works only for two assets... but a multivariate is coming soon! \n")
		return(NULL)
	}

	if(is.null(PTF)){
		
		if(!is.null(R) & !is.null(SIGMA)){
	
			# number of assets
			N = length(R)
			ret = R
			
		} else {
		
			cat("Provide both vector of Ptf returns (R) and volatilities (SIGMA) \n")
			return(NULL)
	
		}
	
	} else {
		## calculate returns and volatilities from Ptf (X)
		# number of assets
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
	FU = function(a,b,af){ a - 0.5 * b * af } 
	
	# original returns and portfolio variance
	or_ret = PtfRet(ret, W)
	or_var = PtfVar(W, SIGMA, vol=TRUE)
	# original utility
	or_uti = FU(or_ret, or_var^2, af)
	
	
	# simulate a series of utility value for a given adversion factor 
	perm = seq(0,1,0.05)
	# Table of simulations
	ST = matrix(0, length(perm)+1, 3)
	colnames(ST) = c("Sim_return", "Sim_variance", "Sim_utility")
	rownames(ST) = 1:(length(perm)+1)
	
	# calculate return, volatility and utility for each generated weights
	i = 2
	while(i <= (length(perm)+1) ){	
		
		ws = c(perm[i-1], 1-perm[i-1] )   
		ST[i,1] = PtfRet(ret, ws) 
		ST[i,2] = PtfVar(w=ws, SIGMA, vol=FALSE) 
		ST[i,3] = ST[i,1] - 0.5 * af * ST[i,2] 
		i = i + 1
	
	}

	ST[1,] = c(or_ret, or_var, or_uti)
	
	# Maximum utility weights
	un = c(1,rep(-1,N-1))
	# optimum risky weights
	#temp = w1[2] + (ret["asset_2"] - ret["asset_1"]) / (af * t(un) %*% SIGMA %*% un)
	temp = ((un*-1) %*% ret + af*(SIGMA[,1] %*% (un))) / (af * t(un) %*% SIGMA %*% un)
	mu_w = c(1-temp, temp)
	names(mu_w) = c("Asset_1","Asset_2")
	
	# Maximum utility Returns
	mu_r = PtfRet(ret, mu_w)
	
	# Maximum utility Volatility
	mu_v = PtfVar(mu_w, SIGMA, vol=FALSE)
	
	# Maximum utility Value
	mu_va = FU(mu_r, mu_v, af)
	
	# plot results 
	if(plot){
		
		# calculate MV optimum portfolio
		PTF_OPT = PtfOpt(ret=ret, SIGMA=SIGMA, volatility = FALSE)
		# retrieve optimum weigths, return and variance
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
	res = list(Max_Utility_Weight = mu_w , Results = data.frame(Max_Utility_Returns = mu_r, Max_Utility_Volatility = mu_v, Max_Utility_Value = mu_va));
	
	res;
	
}