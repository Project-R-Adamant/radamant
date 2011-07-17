###################################################################################
# FUNCTION: Black & Scholes price
#
# SUMMARY:
# Perform Black & Scholes calculation to estimate option price (Call or Put) 
#
# PARAMETERS:
# x = option price
# under = underline asset price
# maty = maturity
# strike = strike price
# rfr = risk free interest rate
# sigma = sample volatility
# opttype = 1 call; 2 put.
#
# RETURNS:
#  Vector of results
#   
###################################################################################
BS.price = function(under,...) UseMethod("BS.price")
## print method for class "BS.price"
print.BS.price = function(x, mod=1, ...){
	# print BS.price
	Logger(message = "print BS.price", from = "print.BS.price", line = 2, level = 1);
	if(mod==1){
		cat(paste("Black & Scholes price", " (", attr(x,"BS_Type"), ")", " for ", attr(x,"Opt_Type"), " option:", "\n\n", sep=""));
		print.default(x[,])
		cat("=============================", "\n")
		cat( "Parameters:", "\n",
			paste("Opt_Type:", attr(x,"Opt_Type")), "\n",
			paste("BS_Type", attr(x,"BS_Type")), "\n", paste("Under", attr(x,"Under")), "\n",
			paste("Strike", attr(x,"Strike")), "\n", paste("RiskFree", attr(x,"RiskFree")), "\n",
			paste("Volatility", attr(x,"Volatility")), "\n", paste("Maturity", attr(x,"Maturity")), "\n",
			paste("Yield", attr(x,"Yield")), "\n"
		)
	} else if(mod==2){
	# print BS.ImpVol
	Logger(message = "print BS.ImpVol", from = "print.BS.price", line = 15, level = 1);
		cat(paste("Black & Scholes Implied Volatility", " (", attr(x,"BS_Type"), ")", " for ", attr(x,"Opt_Type"), " option:", "\n\n", sep=""));	
		print.default(x[,1])
		cat("=============================", "\n")
		cat( "Parameters:", "\n",
			paste("Ref_Price:", paste(attr(x,"Ref_Price"), collapse="; ")), "\n",
			paste("Opt_Type:", attr(x,"Opt_Type")), "\n",
			paste("BS_Type", attr(x,"BS_Type")), "\n", paste("Under", attr(x,"Under")), "\n",
			paste("Strike", attr(x,"Strike")), "\n", paste("RiskFree", attr(x,"RiskFree")), "\n",
			paste("Volatility", attr(x,"Volatility")), "\n", paste("Maturity", attr(x,"Maturity")), "\n",
			paste("Yield", attr(x,"Yield")), "\n"
		)
	} else {
	# print BS.moments
	Logger(message = "print BS.moments", from = "print.BS.price", line = 28, level = 1);
		cat(paste("Black & Scholes Moments", "\n\n"))
		print(x)
	}
}
BS.price.default=function(under, strike, rfr, sigma, maty, yield, calc.type=c("standard","lognorm","gammarec"), opt.type=c("call","put"), ...)  { 
	# option type
	Logger(message = "option type", from = "BS.price.default", line = 2, level = 1);
	opt.type = match.arg(opt.type)
	# BS type
	Logger(message = "BS type", from = "BS.price.default", line = 4, level = 1);
	calc.type = match.arg(calc.type)
	# BS calculation
	Logger(message = "BS calculation", from = "BS.price.default", line = 6, level = 1);
	switch(calc.type,
		"standard" = (res = .BS.price.std(under, strike, rfr, sigma, maty, yield, opt.type))
			,
		"lognorm" =	(res = .BS.price.lgn(under, strike, rfr, sigma, maty, yield, opt.type))
			,
		"gammarec" = (res = .BS.price.gamr(under, strike, rfr, sigma, maty, yield, opt.type))
	)
	# return results
	Logger(message = "return results", from = "BS.price.default", line = 14, level = 1);
	Results = cbind(Price = res[1], Diff_1 = res[2], Diff_2 = res[3]);
	# Assign class and attributes
	Logger(message = "Assign class and attributes", from = "BS.price.default", line = 16, level = 1);
	class(Results) = "BS.price"
	attr(Results, "Opt_Type") = opt.type
	attr(Results, "BS_Type") = calc.type
	attr(Results, "Under") = under
	attr(Results, "Strike") = strike
	attr(Results, "RiskFree") = rfr
	attr(Results, "Volatility") = sigma
	attr(Results, "Maturity") = maty
	attr(Results, "Yield") = yield
	Results;
}
#######################################################################################################################
# FUNCTION: Black & Scholes greeks
#
# SUMMARY:
# Perform Black & Scholes calculation and calculates the greeks 
#
# PARAMETERS:
# x = option price
# under = underline asset price
# maty = maturity
# strike = strike price
# rfr = risk free interest rate
# sigma = sample volatility
#
# RETURNS:
#  List of results
#   
#######################################################################################################################
BS.greeks = function(X=NULL, ...){
	# parameter X of class "BS.price"
	if(class(X) == "BS.price"){
		opt.type = attr(X, "Opt_Type") 
		calc.type = attr(X, "BS_Type") 
		under = attr(X, "Under") 
		strike = attr(X, "Strike") 
		rfr = attr(X, "RiskFree") 
		sigma = attr(X, "Volatility") 
		maty = attr(X, "Maturity") 
		yield = attr(X, "Yield") 
		# BS calculation
		Logger(message = "BS calculation", from = "BS.greeks", line = 12, level = 1);
		bs_price = X[1]
		d1 = X[2]
		d2 = X[3]
	} else {
	# BS parameters specified in ... parameter
	Logger(message = "BS parameters specified in ... parameter", from = "BS.greeks", line = 17, level = 1);
	if(is.null(X)){
		bs = BS.price(...)
		bs_price = bs[1]
		d1 = bs[2] 
		}
	}
	# option type
	Logger(message = "option type", from = "BS.greeks", line = 24, level = 1);
	if(is.null(X))
		opt.type = attr(bs, "Opt_Type")
	# Delta - sensitivity to underlying price
	Logger(message = "Delta - sensitivity to underlying price", from = "BS.greeks", line = 27, level = 1);
	res_delta = pnorm(d1) * exp(-yield*maty);
	# Vega - sensitivity to volatility 
	Logger(message = "Vega - sensitivity to volatility ", from = "BS.greeks", line = 29, level = 1);
	res_vega = under * dnorm(d1) * sqrt(maty) * exp(-yield*maty);
	# Gamma - convexity to underlying price
	Logger(message = "Gamma - convexity to underlying price", from = "BS.greeks", line = 31, level = 1);
	res_gamma = exp(- yield * rfr) * (dnorm(d1) / (under * sigma * sqrt(maty)));
	# Theta - sensitivity to time
	Logger(message = "Theta - sensitivity to time", from = "BS.greeks", line = 33, level = 1);
	res_theta = rfr * bs_price -(rfr-yield) *under * res_delta - 0.5*res_gamma*under^2*sigma^2;
	# Rho - sensitivity to interest rate
	Logger(message = "Rho - sensitivity to interest rate", from = "BS.greeks", line = 35, level = 1);
	res_rho = strike * maty * exp(-rfr * maty) * pnorm(d1 - sigma*sqrt(maty));
	# Lamda - elasticity to underlying price 
	Logger(message = "Lamda - elasticity to underlying price ", from = "BS.greeks", line = 37, level = 1);
	res_lambda = res_delta * (under / bs_price); 
	# list of results
	Logger(message = "list of results", from = "BS.greeks", line = 39, level = 1);
	res = rbind(Delta = res_delta,
				Vega = res_vega,
				Theta = res_theta,
				Rho = res_rho,
				Lambda = res_lambda,
				Gamma = res_gamma);
	cat(paste("Black & Scholes Greeks", " for ", opt.type=opt.type, " option:", "\n\n", sep=""));
	res;
}
#######################################################################################################################
# FUNCTION: Black & Scholes volatility
#
# SUMMARY:
# Perform Black & Scholes calculation to estimate B&S Volatility
#
# PARAMETERS:
# x = option price
# under = underline asset price
# maty = maturity
# strike = strike price
# rfr = risk free interest rate
# sigma = sample volatility
# opttype = 1 call; 2 put.
#
# RETURNS:
#  Vector of results
#   
#######################################################################################################################
BS.ImpVol=function(P, under, strike, rfr, sigma, maty, yield, calc.type=c("standard","lognorm","gammarec"), opt.type=c("call","put"),interval=c(-20, 20)){
	opt.type = match.arg(opt.type)
	calc.type = match.arg(calc.type)
	# allocate matrix of results
	Logger(message = "allocate matrix of results", from = "BS.ImpVol", line = 4, level = 1);
	res = matrix(NA, length(P), 1)
	dimnames(res) = list(paste("Obs_P_", P, sep=""), "Imp_Vol")
	i = 1
	while(i <= length(P)){
		# Black & Scholes difference
		Logger(message = "Black & Scholes difference", from = "BS.ImpVol", line = 9, level = 2);
		BSdiff=function(sigma){ 
			# calculate BS price
			Logger(message = "calculate BS price", from = "BS.ImpVol", line = 11, level = 2);
			bs = BS.price(under, strike, rfr, sigma, maty, yield)
			bs[1] - P[i]
		}
		# calculation	
		Logger(message = "calculation	", from = "BS.ImpVol", line = 15, level = 2);
		res[i,] = uniroot(BSdiff, interval, tol=1/10^12)$root;
		i = i + 1
	}
	# set attributes
	Logger(message = "set attributes", from = "BS.ImpVol", line = 19, level = 1);
	attr(res, "Ref_Price") = P
	attr(res, "Opt_Type") = opt.type
	attr(res, "BS_Type") = calc.type
	attr(res, "Under") = under
	attr(res, "Strike") = strike
	attr(res, "RiskFree") = rfr
	attr(res, "Volatility") = sigma
	attr(res, "Maturity") = maty
	attr(res, "Yield") = yield
	# return results
	Logger(message = "return results", from = "BS.ImpVol", line = 29, level = 1);
	print.BS.price(res, mod=2);
	res[,1]
}
# BS sample moments 
BS.moments = function(BS=NULL, under, rfr, sigma, yield, maty){
	if(class(BS) == "BS.price"){
		under = attr(BS, "Under") 
		rfr = attr(BS, "RiskFree") 
		sigma = attr(BS, "Volatility") 
		maty = attr(BS, "Maturity") 
		yield = attr(BS, "Yield") 
	}
	# matrix of results
	Logger(message = "matrix of results", from = "BS.moments", line = 9, level = 1);
	MOMS = matrix(NA, 4, 1)
	dimnames(MOMS) = list(c("Mom_1","Mom_2","Mean","Var"), c("BS"))
	# calculate moments for BS 
	Logger(message = "calculate moments for BS ", from = "BS.moments", line = 12, level = 1);
	MOMS[3,1] = (log(under) + (rfr-yield - 0.5*sigma^2)*maty)
	MOMS[4,1] = sigma^2 * maty
	MOMS[1,1] = exp(MOMS[3,1] + MOMS[4,1]/2)
	MOMS[2,1] = exp(2*MOMS[3,1] + 2*MOMS[4,1])
	#print.BS.price(MOMS, mod=3);
	Logger(message = "print.BS.price(MOMS, mod=3);", from = "BS.moments", line = 17, level = 1);
	MOMS;
}
# Jarrow and Rudd (JR) Tree 
JR.BinTree = function(Nsteps, p=0.5, under, strike, rfr, sigma, maty, yield, life, ret.steps=FALSE){
	# control input arguments
	Logger(message = "control input arguments", from = "JR.BinTree", line = 2, level = 1);
	Nsteps = as.integer(Nsteps)
	under = as.numeric(under)
	strike = as.numeric(strike)
	rowID = (Nsteps):0L
	colID = 0L:(Nsteps)
	dt = life / Nsteps
	# movement paramaters
	Logger(message = "movement paramaters", from = "JR.BinTree", line = 9, level = 1);
	# up 
	Logger(message = "up ", from = "JR.BinTree", line = 10, level = 1);
	u = exp((rfr-yield - 0.5*sigma^2)*dt + (sigma*sqrt(dt)))
	# down
	Logger(message = "down", from = "JR.BinTree", line = 12, level = 1);
	d = exp((rfr-yield - 0.5*sigma^2)*dt - (sigma*sqrt(dt)))
	# calculate share price path
	Logger(message = "calculate share price path", from = "JR.BinTree", line = 14, level = 1);
	path = StepMat(under, Nsteps+1, u, d)
	# calculate nodal probabilities
	Logger(message = "calculate nodal probabilities", from = "JR.BinTree", line = 16, level = 1);
	nodprob = as.matrix(.BinCoef(Nsteps, rowID) * p^Nsteps)
	# option payoff - European
	Logger(message = "option payoff - European", from = "JR.BinTree", line = 18, level = 1);
	pay = rowMax(cbind(path[,Nsteps+1] - strike, 0))
	# option expected value
	Logger(message = "option expected value", from = "JR.BinTree", line = 20, level = 1);
	ev_opt = as.vector((t(nodprob) %*% pay) * exp(-rfr*life))
	# table of moments
	Logger(message = "table of moments", from = "JR.BinTree", line = 22, level = 1);
	moms = matrix(NA, 4, 1)
	dimnames(moms) = list(c("Mom_1","Mom_2","Mean","Var"), c("JR"))
	# calculate moments for JR
	Logger(message = "calculate moments for JR", from = "JR.BinTree", line = 25, level = 1);
	moms[1,1] = path[,Nsteps+1] %*% nodprob
	moms[2,1] = path[,Nsteps+1]^2 %*% nodprob
	# calculate mean and variance of log share prices
	Logger(message = "calculate mean and variance of log share prices", from = "JR.BinTree", line = 28, level = 1);
	moms[3,1] = 2 * log(moms[1,1]) - log(moms[2,1])/2
	moms[4,1] = log(moms[2,1]) - 2 * log(moms[1,1]) 
	# calculate moments for BS
	Logger(message = "calculate moments for BS", from = "JR.BinTree", line = 31, level = 1);
	MOMS = cbind(moms, BS.moments(NULL, under, rfr, sigma, yield, maty))
	# calculate BS price
	Logger(message = "calculate BS price", from = "JR.BinTree", line = 33, level = 1);
	BSp = BS.price(under, strike, rfr, sigma, maty, yield, opt.type=c("call","put"))[ ,1]
	# table of results
	Logger(message = "table of results", from = "JR.BinTree", line = 35, level = 1);
	Price_tab = cbind(rowID, nodprob, path[,Nsteps+1], pay)
	colnames(Price_tab) = c("Step", "Prob", "Share", "Option")
	# list to return
	Logger(message = "list to return", from = "JR.BinTree", line = 38, level = 1);
	if(!ret.steps)
		Results = list(Price_eval = Price_tab, Moments = MOMS, Values = cbind(B_S = BSp, LR = ev_opt))
	else
		Results = list(Price_eval = Price_tab, Moments = MOMS, Values = cbind(B_S = BSp, LR = ev_opt), Price_Path = path)
	Results;
}
# Cox, Ross, Rubinstein (CRR) Tree 
CRR.BinTree = function(Nsteps, under, strike, rfr, sigma, maty, yield, life, ret.steps=FALSE){
	# control input arguments
	Logger(message = "control input arguments", from = "CRR.BinTree", line = 2, level = 1);
	Nsteps = as.integer(Nsteps)
	under = as.numeric(under)
	strike = as.numeric(strike)
	rowID <- (Nsteps):0L
	colID <- 0L:(Nsteps)
	# discount factor
	Logger(message = "discount factor", from = "CRR.BinTree", line = 8, level = 1);
	dt = life / Nsteps
	disc = exp((rfr - yield)*dt)
	## movement paramaters
	Logger(message = "movement paramaters", from = "CRR.BinTree", line = 11, level = 1);
	# up
	Logger(message = "up", from = "CRR.BinTree", line = 12, level = 1);
	u = exp(sigma * sqrt(dt))
	# down
	Logger(message = "down", from = "CRR.BinTree", line = 14, level = 1);
	d = 1 / u
	# proability
	Logger(message = "proability", from = "CRR.BinTree", line = 16, level = 1);
	p = (disc - d) / (u - d)
	# calculate share price path
	Logger(message = "calculate share price path", from = "CRR.BinTree", line = 18, level = 1);
	path = StepMat(under, Nsteps+1, u, d)
	# calculate nodal probabilities
	Logger(message = "calculate nodal probabilities", from = "CRR.BinTree", line = 20, level = 1);
	nodprob = as.matrix(.BinCoef(Nsteps, rowID) * p^rowID * (1-p)^(Nsteps-rowID))
	#round(nodprob, 3)
	Logger(message = "round(nodprob, 3)", from = "CRR.BinTree", line = 22, level = 1);
	# option payoff - European
	Logger(message = "option payoff - European", from = "CRR.BinTree", line = 23, level = 1);
	pay = rowMax(cbind(path[,Nsteps+1] - strike, 0))
	# option expected value
	Logger(message = "option expected value", from = "CRR.BinTree", line = 25, level = 1);
	ev_opt = as.vector((t(nodprob) %*% pay) * exp(-rfr*life))
	# table of moments
	Logger(message = "table of moments", from = "CRR.BinTree", line = 27, level = 1);
	moms = matrix(NA, 4, 1)
	dimnames(moms) = list(c("Mom_1","Mom_2","Mean","Var"), c("JR"))
	# calculate moments for JR
	Logger(message = "calculate moments for JR", from = "CRR.BinTree", line = 30, level = 1);
	moms[1,1] = path[,Nsteps+1] %*% nodprob
	moms[2,1] = path[,Nsteps+1]^2 %*% nodprob
	# calculate mean and variance of log share prices
	Logger(message = "calculate mean and variance of log share prices", from = "CRR.BinTree", line = 33, level = 1);
	moms[3,1] = 2 * log(moms[1,1]) - log(moms[2,1])/2
	moms[4,1] = log(moms[2,1]) - 2 * log(moms[1,1]) 
	# calculate moments for BS
	Logger(message = "calculate moments for BS", from = "CRR.BinTree", line = 36, level = 1);
	MOMS = cbind(moms, BS.moments(NULL, under, rfr, sigma, yield, maty))
	# calculate BS price
	Logger(message = "calculate BS price", from = "CRR.BinTree", line = 38, level = 1);
	BSp = BS.price(under, strike, rfr, sigma, maty, yield, opt.type=c("call","put"))[ ,1]
	# table of results
	Logger(message = "table of results", from = "CRR.BinTree", line = 40, level = 1);
	Price_tab = cbind(rowID, nodprob, path[,Nsteps+1], pay)
	colnames(Price_tab) = c("Step", "Prob", "Share", "Option")
	# list of results
	Logger(message = "list of results", from = "CRR.BinTree", line = 43, level = 1);
	if(!ret.steps)
		Results = list(Price_eval = Price_tab, Moments = MOMS, Values = cbind(B_S = BSp, LR = ev_opt))
	else
		Results = list(Price_eval = Price_tab, Moments = MOMS, Values = cbind(B_S = BSp, LR = ev_opt), Price_Path = path)
	Results
}
# Leisen - Reimer (LR) Tree 
LR.BinTree = function(Nsteps, under, strike, rfr, sigma, maty, yield, life, ret.steps=FALSE){
	# control input arguments
	Logger(message = "control input arguments", from = "LR.BinTree", line = 2, level = 1);
	Nsteps = as.integer(Nsteps)
	under = as.numeric(under)
	strike = as.numeric(strike)
	rowID <- (Nsteps):0
	colID <- 0:(Nsteps)
	# discount factor
	Logger(message = "discount factor", from = "LR.BinTree", line = 8, level = 1);
	dt = life / Nsteps
	disc = exp((rfr - yield)*dt)
	# calculate BS price
	Logger(message = "calculate BS price", from = "LR.BinTree", line = 11, level = 1);
	BSres = BS.price(under, strike, rfr, sigma, maty, yield, opt.type=c("call","put"))
	BSp = BSres[,1]
	# calculate BS differential factor D1 and D2
	Logger(message = "calculate BS differential factor D1 and D2", from = "LR.BinTree", line = 14, level = 1);
	d1 = BSres[,2]
	d2 = BSres[,3]
	# probabilities
	Logger(message = "probabilities", from = "LR.BinTree", line = 17, level = 1);
	p = .InvPP(d2, Nsteps)
	ps = .InvPP(d1, Nsteps)
	# movement parameters
	Logger(message = "movement parameters", from = "LR.BinTree", line = 20, level = 1);
	# up
	Logger(message = "up", from = "LR.BinTree", line = 21, level = 1);
	u = disc * ps / p
	# down
	Logger(message = "down", from = "LR.BinTree", line = 23, level = 1);
	d = disc * (1-ps) / (1-p)
	# calculate share price path
	Logger(message = "calculate share price path", from = "LR.BinTree", line = 25, level = 1);
	path = StepMat(under, Nsteps+1, u, d)
	# calculate nodal probabilities
	Logger(message = "calculate nodal probabilities", from = "LR.BinTree", line = 27, level = 1);
	nodprob = as.matrix(.BinCoef(Nsteps, rowID) * p^rowID * (1-p)^(Nsteps-rowID))
	#round(nodprob, 3)
	Logger(message = "round(nodprob, 3)", from = "LR.BinTree", line = 29, level = 1);
	# option payoff - European
	Logger(message = "option payoff - European", from = "LR.BinTree", line = 30, level = 1);
	pay = rowMax(cbind(path[,Nsteps+1] - strike, 0))
	# option expected value
	Logger(message = "option expected value", from = "LR.BinTree", line = 32, level = 1);
	ev_opt = as.vector((t(nodprob) %*% pay) * exp(-rfr*maty))
	# table of moments
	Logger(message = "table of moments", from = "LR.BinTree", line = 34, level = 1);
	moms = matrix(NA, 4, 1)
	dimnames(moms) = list(c("Mom_1","Mom_2","Mean","Var"), c("JR"))
	# calculate moments for JR
	Logger(message = "calculate moments for JR", from = "LR.BinTree", line = 37, level = 1);
	moms[1,1] = path[,Nsteps+1] %*% nodprob
	moms[2,1] = path[,Nsteps+1]^2 %*% nodprob
	# calculate mean and variance of log share prices
	Logger(message = "calculate mean and variance of log share prices", from = "LR.BinTree", line = 40, level = 1);
	moms[3,1] = 2 * log(moms[1,1]) - log(moms[2,1])/2
	moms[4,1] = log(moms[2,1]) - 2 * log(moms[1,1]) 
	# calculate moments for BS
	Logger(message = "calculate moments for BS", from = "LR.BinTree", line = 43, level = 1);
	MOMS = cbind(moms, BS.moments(NULL, under, rfr, sigma, yield, maty))
	# table of results
	Logger(message = "table of results", from = "LR.BinTree", line = 45, level = 1);
	Price_tab = cbind(rowID, nodprob, path[,Nsteps+1], pay)
	colnames(Price_tab) = c("Step", "Prob", "Share", "Option")
	# list of results
	Logger(message = "list of results", from = "LR.BinTree", line = 48, level = 1);
	if(!ret.steps)
		Results = list(Price_eval = Price_tab, Moments = MOMS, Values = cbind(B_S = BSp, LR = ev_opt))
	else
		Results = list(Price_eval = Price_tab, Moments = MOMS, Values = cbind(B_S = BSp, LR = ev_opt), Price_Path = path)
	Results
}
# Edgeworth option price calculation
Edgeworth.price = function(init, under, strike, rfr, sigma, maty, yield){
	init = as.matrix(init)
	Nsteps = nrow(init)-1 
	rownames(init) = 0:(Nsteps)
	# call edgeworth distribution
	Logger(message = "call edgeworth distribution", from = "Edgeworth.price", line = 5, level = 1);
	EW = EdgeWorthDist(init, Nsteps)
	# risk neutral adjustment for edgeworth prob
	Logger(message = "risk neutral adjustment for edgeworth prob", from = "Edgeworth.price", line = 7, level = 1);
	edge_rn = rfr - yield - log(sum(EW[[1]][,2]*exp(sigma*sqrt(maty)*EW[[1]][,1])))/maty
	# calculate share price
	Logger(message = "calculate share price", from = "Edgeworth.price", line = 9, level = 1);
	share_price = under*exp(edge_rn*maty+EW[[1]][,1]*sigma*sqrt(maty))
	# calculate payoff
	Logger(message = "calculate payoff", from = "Edgeworth.price", line = 11, level = 1);
	edge_pay = rowMax(cbind(0, share_price - strike)) * EW[[1]][,2]
	# option expected value
	Logger(message = "option expected value", from = "Edgeworth.price", line = 13, level = 1);
	ev_opt = sum(edge_pay) * exp(-rfr*maty)
	ev_opt;
}
# Edgeworth distribution
EdgeWorthDist = function(init, Nsteps, p=0.5){
	rowID = as.numeric(rownames(init))
	# calculate nodal probabilities
	Logger(message = "calculate nodal probabilities", from = "EdgeWorthDist", line = 3, level = 1);
	nodprob = as.matrix(.BinCoef(Nsteps, rowID) * p^rowID * (1-p)^(Nsteps-rowID))
	edgf = .EdgeFact(init, 0, 5.4)
	adj_prob = edgf * nodprob
	edg_prob = adj_prob / sum(adj_prob)
	# initial edge moments
	Logger(message = "initial edge moments", from = "EdgeWorthDist", line = 8, level = 1);
	moms_i = SampMom(edg_prob, init, moms=1:4)
	std_init = Zscore(init, moms_i[1], sqrt(moms_i[2]))
	# final edge moments
	Logger(message = "final edge moments", from = "EdgeWorthDist", line = 11, level = 1);
	moms_f = SampMom(edg_prob, std_init, 1:4)
	#results
	Logger(message = "results", from = "EdgeWorthDist", line = 13, level = 1);
	list(cbind(Nodes=std_init, Edge_prob=edg_prob) ,Moments=moms_f);
}
# Browniam motion simulation
BroMot = function(nsim, T, S0=0, mi=0, sigma=1, geom=TRUE, same.rnd=TRUE, plot=FALSE, ...){
	mi = as.matrix(mi)
	sigma = as.matrix(sigma)	
	S0 = as.matrix(S0)
	# time frame
	Logger(message = "time frame", from = "BroMot", line = 5, level = 1);
	if(missing(T))
		T = nsim
	nser = nrow(mi)
	# check parameters length
	Logger(message = "check parameters length", from = "BroMot", line = 9, level = 1);
	if(!identical(nrow(mi), nrow(sigma), nrow(S0))){
		cat("'~' Input parameters length differs \n")
		return(NULL)
	}
	# drift factor
	Logger(message = "drift factor", from = "BroMot", line = 14, level = 1);
	if(geom)
		alpha = mi - 0.5 * sigma^2 else
			alpha = mi
	# time increment	
	Logger(message = "time increment	", from = "BroMot", line = 18, level = 1);
	t = (seq(0, 1, length.out=nsim) / nsim) * T	
	# simulate from std normal
	Logger(message = "simulate from std normal", from = "BroMot", line = 20, level = 1);
	# same random path for each series
	Logger(message = "same random path for each series", from = "BroMot", line = 21, level = 1);
	if(same.rnd)
		Z = cumsum(c(0, rnorm(nsim-1, 0, 1))) / nsim
	# simulate nsim-return following (Geom) BM
	Logger(message = "simulate nsim-return following (Geom) BM", from = "BroMot", line = 24, level = 1);
	S = matrix(NA, nsim, nser)
	colnames(S) = paste("BM_", 1:nser, "-- Drift=", alpha, ", Sigma=", sigma, sep="" )
	i = 1
	if(same.rnd){
		# same random path for each series
		Logger(message = "same random path for each series", from = "BroMot", line = 29, level = 1);
		while(i <= nser){
			S[ ,i] = S0[i,] + alpha[i,]*t + sigma[i,] * sqrt(T) * Z
			i = i + 1
		}
	} else {
		while(i <= nser){
			# different random path for each series
			Logger(message = "different random path for each series", from = "BroMot", line = 36, level = 2);
			Z = cumsum(c(0, rnorm(nsim-1, 0, 1))) / nsim
			S[ ,i] = S0[i,] + alpha[i,]*t + sigma[i,] * sqrt(T) * Z
			i = i + 1
		}
	}		
	# check if geometric BM
	Logger(message = "check if geometric BM", from = "BroMot", line = 42, level = 1);
	if(geom)
		S = exp(S) 
	# plot simulated series
	Logger(message = "plot simulated series", from = "BroMot", line = 45, level = 1);
	if(plot){
		cplot(S, base= mean(S0) + t*T, main="Brownian Motion", ...)
		if(nser == 1){
			if(geom)			
				X = exp(mi*t) else
					X = mi*t
			cplot(X, base = t*T, col="red", append=TRUE, show.legend=FALSE, ...)
		}
	}
	# return simulated series
	Logger(message = "return simulated series", from = "BroMot", line = 55, level = 1);
	invisible(S)
}
# Browniam motion simulation bi-dimensional
BroMot2D = function(nsim, T, S0, mi, sigma, geom=TRUE, same.rnd=FALSE, laydisp=NULL, plot=TRUE, ...){
	mi = as.matrix(mi)
	sigma = as.matrix(sigma)	
	S0 = as.matrix(S0)
	nser = nrow(mi)
	# time frame
	Logger(message = "time frame", from = "BroMot2D", line = 6, level = 1);
	if(missing(T))
		T = nsim
	if(!identical(nrow(mi), nrow(sigma), nrow(S0))){
		cat("'~' Input parameters length differs \n")
		return(NULL)
	}
	BM = BroMot(nsim=nsim, T=T, S0=S0, mi=mi, sigma=sigma, geom=geom, same.rnd, plot=FALSE, ...)
	if(nser > 2){
		pos = .BinCoef(nser, 2)
		elem = combn(nser, 2)
		if(is.null(laydisp))
			laydisp = c(floor(pos/2), ceiling(pos/2)-1)
		par(mfrow=laydisp)
		nn = paste("BM_", 1:nser, "-Drift=", mi, ",Sigma=", sigma, sep="" )
		for(i in 1:pos)
			cplot(BM[,elem[1,i]], BM[,elem[2,i]], main=paste(nn[elem[1,i]], nn[elem[2,i]], sep="<-->"), cex.main=0.5, cex.axis=0.5, ...)
	} else {
		nn = paste("BM_", 1:nser, "-- Drift=", mi, ", Sigma=", sigma, sep="" )
		cplot(BM[,1], (seq(0, 1, length.out=nsim) / nsim), main = nn, ...)
	}
	invisible(BM)
}
# Probability density function for hitting barriers
PDFHit = function(t, B=0, S0=0, mi, sigma, cumul=FALSE, plot=FALSE, ...){
	# check if barrier is above starting point
	Logger(message = "check if barrier is above starting point", from = "PDFHit", line = 2, level = 1);
	if(B <= S0){
		a = S0
		b = a + mi*t
	} else {
		a = B - S0
		b = a - mi*t
	}
	if(cumul){
		# cumulative probability distribution
		Logger(message = "cumulative probability distribution", from = "PDFHit", line = 11, level = 1);
		pd = dnorm((-a + mi*t) / (sigma*sqrt(t))) + exp((2*mi*a) / (sigma^2)) * dnorm((-a-mi*t) / (sigma*sqrt(t)))
	} else {
		# probability distribution
		Logger(message = "probability distribution", from = "PDFHit", line = 14, level = 1);
		pd = (a / (sigma*sqrt(2*pi*t^3))) * exp(-b^2 / (2*sigma^2*t))
	}
	if(plot){
		# main title
		Logger(message = "main title", from = "PDFHit", line = 18, level = 1);
		main = paste( ifelse(cumul,"CDF","PDF"), "-> N=", length(t), "; Mi=",mi, "; Sigma=",sigma, "; S0=", S0, "; B=", B, sep="")
		# plot function
		Logger(message = "plot function", from = "PDFHit", line = 20, level = 1);
		cplot(pd, main=main, xtitle="t", ytitle=ifelse(cumul, expression(G(t,S0,B)),expression(g(t,S0,B))), ...)
	}
	invisible(pd)
}
# Probability of hitting barrier
ProbHit = function(B, S0, mi, sigma){
	# double barrier
	Logger(message = "double barrier", from = "ProbHit", line = 2, level = 1);
	if(length(B) > 1){
		if(mi == 0){
			ph = (B[2] - S0) / (B[2] - B[1])
		} else {
			ph = ( exp(-(2*mi*S0)/sigma^2) - exp(-(2*mi*B[1])/sigma^2) ) / (exp(-(2*mi*B[2])/sigma^2) - exp(-(2*mi*B[2])/sigma^2) )
		}
	} else {
		# barrier below starting point
		Logger(message = "barrier below starting point", from = "ProbHit", line = 10, level = 1);
		if(B < S0){
			a = S0 - B
		# barrier above starting point
		Logger(message = "barrier above starting point", from = "ProbHit", line = 13, level = 1);
		} else if(B > S0){
			a = -(B - S0)
		# barrier equal to starting point	
		Logger(message = "barrier equal to starting point	", from = "ProbHit", line = 16, level = 1);
		} else {
			a = S0
		}
		# calculate probability
		Logger(message = "calculate probability", from = "ProbHit", line = 20, level = 1);
		ph = exp(-(2*mi*a) / sigma^2);
	}
	# return probability
	Logger(message = "return probability", from = "ProbHit", line = 23, level = 1);
	cat("Probability to reach Absorption Barrier (5 digits approx): \n")
	if(ph <= 1e-8)
		0 
	else if(ph > 1)
		1
	else 
		ph
}
# First hit time
FirstHit = function(B, S0, mi, sigma, geom=FALSE, nsim=500, plot=FALSE){
	# barrier below starting point
	Logger(message = "barrier below starting point", from = "FirstHit", line = 2, level = 1);
	if(B < S0){
		if(geom){
			if(mi < 0.5*sigma^2)
				ht = (1/(0.5*sigma^2-mi)) * log(S0/B) else
					ht = Inf
		} else {
			if(mi<0)
				ht = (S0 - B) / abs(mi) else 
					ht = Inf
		}
	# barrier above starting point
	Logger(message = "barrier above starting point", from = "FirstHit", line = 13, level = 1);
	} else if(B > S0){
		if(geom){
			if(mi > 0.5*sigma^2)
				ht = (1/(mi-0.5*sigma^2)) * log(B/S0) else
					ht = Inf
		} else {
			if(mi>0)		
				ht = (B - S0) / abs(mi) else 
					ht = Inf
		}
	# barrier equal to starting point	
	Logger(message = "barrier equal to starting point	", from = "FirstHit", line = 24, level = 1);
	} else {
		if(geom){
			if(mi < 0.5*sigma^2)
				ht = (1/(mi-0.5*sigma^2)) * log(S0) else
					ht = Inf
		} else {
			if(mi<0)
				ht = S0 / abs(mi) else
				ht = Inf
		}
	}
	# plot simulated Brownian Motion
	Logger(message = "plot simulated Brownian Motion", from = "FirstHit", line = 36, level = 1);
	if(plot){
		BroMot(nsim=nsim, nsim*B, S0=S0, mi=mi, sigma=sigma, geom=geom, same.rnd=TRUE, plot=TRUE)
		abline(h=B, col="green", lwd=2, lty=4)
	}
	# return expected time
	Logger(message = "return expected time", from = "FirstHit", line = 41, level = 1);
	cat("Expected First Hitting Time: \n")
	ht;
}
# Create Step probability matrix 
StepMat = function(init, n_step, up, down){
	# assign matrix of steps and assign dimnames corresponding to steps number
	Logger(message = "assign matrix of steps and assign dimnames corresponding to steps number", from = "StepMat", line = 2, level = 1);
	mm = matrix(0, n_step, n_step)
	rowID = (n_step-1):0
	colID = 0:(n_step-1)
	dimnames(mm) = list(rowID , colID)
	mm[n_step, 1] = init
	I = n_step
	j = 2
	# loop through columns
	Logger(message = "loop through columns", from = "StepMat", line = 10, level = 1);
	while(j <= n_step){
		# loop through rows
		Logger(message = "loop through rows", from = "StepMat", line = 12, level = 2);
		row <- I : (I-j+1)
		for(i in row)
			mm[i, j] = ifelse(rowID[i] < colID[j] , mm[i, j-1]*down, mm[i+1, j-1]*up)
		j = j + 1
	}
	# return step matrix
	Logger(message = "return step matrix", from = "StepMat", line = 18, level = 1);
	mm;
}
