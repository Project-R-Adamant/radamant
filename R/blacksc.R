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
		cat(paste("Black & Scholes Moments", "\n\n"))
		print(x)
	}
}
BS.price.default=function(under, strike, rfr, sigma, maty, yield, calc.type=c("standard","lognorm","gammarec"), opt.type=c("call","put"), ...)  { 
	# option type
	opt.type = match.arg(opt.type)
	# BS type
	calc.type = match.arg(calc.type)
	# BS calculation
	switch(calc.type,
		"standard" = (res = .BS.price.std(under, strike, rfr, sigma, maty, yield, opt.type))
			,
		"lognorm" =	(res = .BS.price.lgn(under, strike, rfr, sigma, maty, yield, opt.type))
			,
		"gammarec" = (res = .BS.price.gamr(under, strike, rfr, sigma, maty, yield, opt.type))
	)
	# return results
	Results = cbind(Price = res[1], Diff_1 = res[2], Diff_2 = res[3]);
	# Assign class and attributes
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
		bs_price = X[1]
		d1 = X[2]
		d2 = X[3]
	} else {
	# BS parameters specified in ... parameter
	if(is.null(X)){
		bs = BS.price(...)
		bs_price = bs[1]
		d1 = bs[2] 
		}
	}
	# option type
	if(is.null(X))
		opt.type = attr(bs, "Opt_Type")
	# Delta - sensitivity to underlying price
	res_delta = pnorm(d1) * exp(-yield*maty);
	# Vega - sensitivity to volatility 
	res_vega = under * dnorm(d1) * sqrt(maty) * exp(-yield*maty);
	# Gamma - convexity to underlying price
	res_gamma = exp(- yield * rfr) * (dnorm(d1) / (under * sigma * sqrt(maty)));
	# Theta - sensitivity to time
	res_theta = rfr * bs_price -(rfr-yield) *under * res_delta - 0.5*res_gamma*under^2*sigma^2;
	# Rho - sensitivity to interest rate
	res_rho = strike * maty * exp(-rfr * maty) * pnorm(d1 - sigma*sqrt(maty));
	# Lamda - elasticity to underlying price 
	res_lambda = res_delta * (under / bs_price); 
	# list of results
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
	res = matrix(NA, length(P), 1)
	dimnames(res) = list(paste("Obs_P_", P, sep=""), "Imp_Vol")
	i = 1
	while(i <= length(P)){
		# Black & Scholes difference
		BSdiff=function(sigma){ 
			# calculate BS price
			bs = BS.price(under, strike, rfr, sigma, maty, yield)
			bs[1] - P[i]
		}
		# calculation	
		res[i,] = uniroot(BSdiff, interval, tol=1/10^12)$root;
		i = i + 1
	}
	# set attributes
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
	MOMS = matrix(NA, 4, 1)
	dimnames(MOMS) = list(c("Mom_1","Mom_2","Mean","Var"), c("BS"))
	# calculate moments for BS 
	MOMS[3,1] = (log(under) + (rfr-yield - 0.5*sigma^2)*maty)
	MOMS[4,1] = sigma^2 * maty
	MOMS[1,1] = exp(MOMS[3,1] + MOMS[4,1]/2)
	MOMS[2,1] = exp(2*MOMS[3,1] + 2*MOMS[4,1])
	#print.BS.price(MOMS, mod=3);
	MOMS;
}
# Jarrow and Rudd (JR) Tree 
JR.BinTree = function(Nsteps, p=0.5, under, strike, rfr, sigma, maty, yield, life, ret.steps=FALSE){
	# control input arguments
	Nsteps = as.integer(Nsteps)
	under = as.numeric(under)
	strike = as.numeric(strike)
	rowID = (Nsteps):0L
	colID = 0L:(Nsteps)
	dt = life / Nsteps
	# movement paramaters
	# up 
	u = exp((rfr-yield - 0.5*sigma^2)*dt + (sigma*sqrt(dt)))
	# down
	d = exp((rfr-yield - 0.5*sigma^2)*dt - (sigma*sqrt(dt)))
	# calculate share price path
	path = StepMat(under, Nsteps+1, u, d)
	# calculate nodal probabilities
	nodprob = as.matrix(.BinCoef(Nsteps, rowID) * p^Nsteps)
	# option payoff - European
	pay = rowMax(cbind(path[,Nsteps+1] - strike, 0))
	# option expected value
	ev_opt = as.vector((t(nodprob) %*% pay) * exp(-rfr*life))
	# table of moments
	moms = matrix(NA, 4, 1)
	dimnames(moms) = list(c("Mom_1","Mom_2","Mean","Var"), c("JR"))
	# calculate moments for JR
	moms[1,1] = path[,Nsteps+1] %*% nodprob
	moms[2,1] = path[,Nsteps+1]^2 %*% nodprob
	# calculate mean and variance of log share prices
	moms[3,1] = 2 * log(moms[1,1]) - log(moms[2,1])/2
	moms[4,1] = log(moms[2,1]) - 2 * log(moms[1,1]) 
	# calculate moments for BS
	MOMS = cbind(moms, BS.moments(NULL, under, rfr, sigma, yield, maty))
	# calculate BS price
	BSp = BS.price(under, strike, rfr, sigma, maty, yield, opt.type=c("call","put"))[ ,1]
	# table of results
	Price_tab = cbind(rowID, nodprob, path[,Nsteps+1], pay)
	colnames(Price_tab) = c("Step", "Prob", "Share", "Option")
	# list to return
	if(!ret.steps)
		Results = list(Price_eval = Price_tab, Moments = MOMS, Values = cbind(B_S = BSp, LR = ev_opt))
	else
		Results = list(Price_eval = Price_tab, Moments = MOMS, Values = cbind(B_S = BSp, LR = ev_opt), Price_Path = path)
	Results;
}
# Cox, Ross, Rubinstein (CRR) Tree 
CRR.BinTree = function(Nsteps, under, strike, rfr, sigma, maty, yield, life, ret.steps=FALSE){
	# control input arguments
	Nsteps = as.integer(Nsteps)
	under = as.numeric(under)
	strike = as.numeric(strike)
	rowID <- (Nsteps):0L
	colID <- 0L:(Nsteps)
	# discount factor
	dt = life / Nsteps
	disc = exp((rfr - yield)*dt)
	## movement paramaters
	# up
	u = exp(sigma * sqrt(dt))
	# down
	d = 1 / u
	# proability
	p = (disc - d) / (u - d)
	# calculate share price path
	path = StepMat(under, Nsteps+1, u, d)
	# calculate nodal probabilities
	nodprob = as.matrix(.BinCoef(Nsteps, rowID) * p^rowID * (1-p)^(Nsteps-rowID))
	#round(nodprob, 3)
	# option payoff - European
	pay = rowMax(cbind(path[,Nsteps+1] - strike, 0))
	# option expected value
	ev_opt = as.vector((t(nodprob) %*% pay) * exp(-rfr*life))
	# table of moments
	moms = matrix(NA, 4, 1)
	dimnames(moms) = list(c("Mom_1","Mom_2","Mean","Var"), c("JR"))
	# calculate moments for JR
	moms[1,1] = path[,Nsteps+1] %*% nodprob
	moms[2,1] = path[,Nsteps+1]^2 %*% nodprob
	# calculate mean and variance of log share prices
	moms[3,1] = 2 * log(moms[1,1]) - log(moms[2,1])/2
	moms[4,1] = log(moms[2,1]) - 2 * log(moms[1,1]) 
	# calculate moments for BS
	MOMS = cbind(moms, BS.moments(NULL, under, rfr, sigma, yield, maty))
	# calculate BS price
	BSp = BS.price(under, strike, rfr, sigma, maty, yield, opt.type=c("call","put"))[ ,1]
	# table of results
	Price_tab = cbind(rowID, nodprob, path[,Nsteps+1], pay)
	colnames(Price_tab) = c("Step", "Prob", "Share", "Option")
	# list of results
	if(!ret.steps)
		Results = list(Price_eval = Price_tab, Moments = MOMS, Values = cbind(B_S = BSp, LR = ev_opt))
	else
		Results = list(Price_eval = Price_tab, Moments = MOMS, Values = cbind(B_S = BSp, LR = ev_opt), Price_Path = path)
	Results
}
# Leisen - Reimer (LR) Tree 
LR.BinTree = function(Nsteps, under, strike, rfr, sigma, maty, yield, life, ret.steps=FALSE){
	# control input arguments
	Nsteps = as.integer(Nsteps)
	under = as.numeric(under)
	strike = as.numeric(strike)
	rowID <- (Nsteps):0
	colID <- 0:(Nsteps)
	# discount factor
	dt = life / Nsteps
	disc = exp((rfr - yield)*dt)
	# calculate BS price
	BSres = BS.price(under, strike, rfr, sigma, maty, yield, opt.type=c("call","put"))
	BSp = BSres[,1]
	# calculate BS differential factor D1 and D2
	d1 = BSres[,2]
	d2 = BSres[,3]
	# probabilities
	p = .InvPP(d2, Nsteps)
	ps = .InvPP(d1, Nsteps)
	# movement parameters
	# up
	u = disc * ps / p
	# down
	d = disc * (1-ps) / (1-p)
	# calculate share price path
	path = StepMat(under, Nsteps+1, u, d)
	# calculate nodal probabilities
	nodprob = as.matrix(.BinCoef(Nsteps, rowID) * p^rowID * (1-p)^(Nsteps-rowID))
	#round(nodprob, 3)
	# option payoff - European
	pay = rowMax(cbind(path[,Nsteps+1] - strike, 0))
	# option expected value
	ev_opt = as.vector((t(nodprob) %*% pay) * exp(-rfr*maty))
	# table of moments
	moms = matrix(NA, 4, 1)
	dimnames(moms) = list(c("Mom_1","Mom_2","Mean","Var"), c("JR"))
	# calculate moments for JR
	moms[1,1] = path[,Nsteps+1] %*% nodprob
	moms[2,1] = path[,Nsteps+1]^2 %*% nodprob
	# calculate mean and variance of log share prices
	moms[3,1] = 2 * log(moms[1,1]) - log(moms[2,1])/2
	moms[4,1] = log(moms[2,1]) - 2 * log(moms[1,1]) 
	# calculate moments for BS
	MOMS = cbind(moms, BS.moments(NULL, under, rfr, sigma, yield, maty))
	# table of results
	Price_tab = cbind(rowID, nodprob, path[,Nsteps+1], pay)
	colnames(Price_tab) = c("Step", "Prob", "Share", "Option")
	# list of results
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
	EW = EdgeWorthDist(init, Nsteps)
	# risk neutral adjustment for edgeworth prob
	edge_rn = rfr - yield - log(sum(EW[[1]][,2]*exp(sigma*sqrt(maty)*EW[[1]][,1])))/maty
	# calculate share price
	share_price = under*exp(edge_rn*maty+EW[[1]][,1]*sigma*sqrt(maty))
	# calculate payoff
	edge_pay = rowMax(cbind(0, share_price - strike)) * EW[[1]][,2]
	# option expected value
	ev_opt = sum(edge_pay) * exp(-rfr*maty)
	ev_opt;
}
# Edgeworth distribution
EdgeWorthDist = function(init, Nsteps, p=0.5){
	rowID = as.numeric(rownames(init))
	# calculate nodal probabilities
	nodprob = as.matrix(.BinCoef(Nsteps, rowID) * p^rowID * (1-p)^(Nsteps-rowID))
	edgf = .EdgeFact(init, 0, 5.4)
	adj_prob = edgf * nodprob
	edg_prob = adj_prob / sum(adj_prob)
	# initial edge moments
	moms_i = SampMom(edg_prob, init, moms=1:4)
	std_init = Zscore(init, moms_i[1], sqrt(moms_i[2]))
	# final edge moments
	moms_f = SampMom(edg_prob, std_init, 1:4)
	#results
	list(cbind(Nodes=std_init, Edge_prob=edg_prob) ,Moments=moms_f);
}
# Browniam motion simulation
BroMot = function(nsim, T, S0=0, mi=0, sigma=1, geom=TRUE, same.rnd=TRUE, plot=FALSE, ...){
	mi = as.matrix(mi)
	sigma = as.matrix(sigma)	
	S0 = as.matrix(S0)
	# time frame
	if(missing(T))
		T = nsim
	nser = nrow(mi)
	# check parameters length
	if(!identical(nrow(mi), nrow(sigma), nrow(S0))){
		cat("'~' Input parameters length differs \n")
		return(NULL)
	}
	# drift factor
	if(geom)
		alpha = mi - 0.5 * sigma^2 else
			alpha = mi
	# time increment	
	t = (seq(0, 1, length.out=nsim) / nsim) * T	
	# simulate from std normal
	# same random path for each series
	if(same.rnd)
		Z = cumsum(c(0, rnorm(nsim-1, 0, 1))) / nsim
	# simulate nsim-return following (Geom) BM
	S = matrix(NA, nsim, nser)
	colnames(S) = paste("BM_", 1:nser, "-- Drift=", alpha, ", Sigma=", sigma, sep="" )
	i = 1
	if(same.rnd){
		# same random path for each series
		while(i <= nser){
			S[ ,i] = S0[i,] + alpha[i,]*t + sigma[i,] * sqrt(T) * Z
			i = i + 1
		}
	} else {
		while(i <= nser){
			# different random path for each series
			Z = cumsum(c(0, rnorm(nsim-1, 0, 1))) / nsim
			S[ ,i] = S0[i,] + alpha[i,]*t + sigma[i,] * sqrt(T) * Z
			i = i + 1
		}
	}		
	# check if geometric BM
	if(geom)
		S = exp(S) 
	# plot simulated series
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
	invisible(S)
}
# Browniam motion simulation bi-dimensional
BroMot2D = function(nsim, T, S0, mi, sigma, geom=TRUE, same.rnd=FALSE, laydisp=NULL, plot=TRUE, ...){
	mi = as.matrix(mi)
	sigma = as.matrix(sigma)	
	S0 = as.matrix(S0)
	nser = nrow(mi)
	# time frame
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
	if(B <= S0){
		a = S0
		b = a + mi*t
	} else {
		a = B - S0
		b = a - mi*t
	}
	if(cumul){
		# cumulative probability distribution
		pd = dnorm((-a + mi*t) / (sigma*sqrt(t))) + exp((2*mi*a) / (sigma^2)) * dnorm((-a-mi*t) / (sigma*sqrt(t)))
	} else {
		# probability distribution
		pd = (a / (sigma*sqrt(2*pi*t^3))) * exp(-b^2 / (2*sigma^2*t))
	}
	if(plot){
		# main title
		main = paste( ifelse(cumul,"CDF","PDF"), "-> N=", length(t), "; Mi=",mi, "; Sigma=",sigma, "; S0=", S0, "; B=", B, sep="")
		# plot function
		cplot(pd, main=main, xtitle="t", ytitle=ifelse(cumul, expression(G(t,S0,B)),expression(g(t,S0,B))), ...)
	}
	invisible(pd)
}
# Probability of hitting barrier
ProbHit = function(B, S0, mi, sigma){
	# double barrier
	if(length(B) > 1){
		if(mi == 0){
			ph = (B[2] - S0) / (B[2] - B[1])
		} else {
			ph = ( exp(-(2*mi*S0)/sigma^2) - exp(-(2*mi*B[1])/sigma^2) ) / (exp(-(2*mi*B[2])/sigma^2) - exp(-(2*mi*B[2])/sigma^2) )
		}
	} else {
		# barrier below starting point
		if(B < S0){
			a = S0 - B
		# barrier above starting point
		} else if(B > S0){
			a = -(B - S0)
		# barrier equal to starting point	
		} else {
			a = S0
		}
		# calculate probability
		ph = exp(-(2*mi*a) / sigma^2);
	}
	# return probability
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
	if(plot){
		BroMot(nsim=nsim, nsim*B, S0=S0, mi=mi, sigma=sigma, geom=geom, same.rnd=TRUE, plot=TRUE)
		abline(h=B, col="green", lwd=2, lty=4)
	}
	# return expected time
	cat("Expected First Hitting Time: \n")
	ht;
}
# Create Step probability matrix 
StepMat = function(init, n_step, up, down){
	# assign matrix of steps and assign dimnames corresponding to steps number
	mm = matrix(0, n_step, n_step)
	rowID = (n_step-1):0
	colID = 0:(n_step-1)
	dimnames(mm) = list(rowID , colID)
	mm[n_step, 1] = init
	I = n_step
	j = 2
	# loop through columns
	while(j <= n_step){
		
		# loop through rows
		row <- I : (I-j+1)
		for(i in row)
			mm[i, j] = ifelse(rowID[i] < colID[j] , mm[i, j-1]*down, mm[i+1, j-1]*up)
		j = j + 1
	}
	# return step matrix
	mm;
}