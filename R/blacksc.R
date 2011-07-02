###################################################################################
# FUNCTION: Black & Scholes price
#
# AUTHOR: FM, RCC
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
BS.price.default=function(under, strike, rfr, sigma, maty, yield, calc.type=c("standard","lognorm","gammarec"), opt.type=c("call","put"), ...)  { 
	# option type
	opt.type = match.arg(opt.type)
	
	# BS type
	calc.type = match.arg(calc.type)
	
	# BS calculation
	switch(calc.type)(

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
	
	cat(paste("Black & Scholes price", " (", calc.type, ")", " for ", opt.type, " option:", "\n\n", sep=""));
	Results;
}
#######################################################################################################################
# FUNCTION: Black & Scholes greeks
#
# AUTHOR: FM, RCC
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

	if(is.null(X)){
		bs = BS.price(...)
		bs_price = bs[1]
		d1 = bs[2] 
		}
	}
			
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

	cat(paste("Black & Scholes Greeks", " for ", opt.type, " option:", "\n\n", sep=""));
	res;

}
#######################################################################################################################
# FUNCTION: Black & Scholes volatility
#
# AUTHOR: FM, RCC
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
BS.ImpVol=function(P, under, strike, rfr, sigma, maty, yield, interval=c(-20, 20), calc.type=c("standard","lognorm","gammarec"), opt.type=c("call","put")){
	
	opt.type = match.arg(opt.type)
	calc.type = match.arg(calc.type)
		
	BSdiff=function(sigma){ 
		bs = BS.price(under, strike, rfr, sigma, maty, yield, calc.type, calc.type, opt.type)
		bs[1] - P
		}
	
	# calculation	
	res = uniroot(BSdiff, interval, tol=1/10^12)$root;

	cat(paste("Black & Scholes Implied Volatility", " (", calc.type, ")", " for ", opt.type, " option:", "\n\n", sep=""));
	res;

}
# BS formula
BS.formula=function(type=c("call","put")){ 
	# different calculation for "call" and "put" options
	if (match.arg(type) == "call") 
		res = expression(under*exp(-yield*maty)*pnorm(log(under/strike) + (((rfr-yield) + 0.5*sigma^2)*maty) / (sigma*sqrt(maty))) - strike*exp(-rfr*maty)*pnorm(log(under/strike) + (((rfr-yield) + 0.5*sigma^2)*maty) / (sigma*sqrt(maty))- sigma*sqrt(maty)))
	else
		res = expression(under*exp(-yield*maty)*pnorm(-(log(under/strike) + (((rfr-yield) + 0.5*sigma^2)*maty) / (sigma*sqrt(maty)))) - strike*exp(-rfr*maty)*pnorm(-(log(under/strike) + (((rfr-yield) + 0.5*sigma^2)*maty) / (sigma*sqrt(maty))- sigma*sqrt(maty))) - under*exp(-yield*maty)+strike*exp(-rfr*maty))
	
	cat(paste("Black & Scholes for", match.arg(type), "options:", "\n\n"))
	res
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
	
	cat(paste("Black & Scholes Moments", "\n\n"))
	MOMS;
}
# Jarrow and Rudd (JR) Tree 
JR.BinTree = function(Nsteps, p, under, strike, rfr, sigma, maty, yield, life, ret.steps=FALSE){

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
	nodprob = as.matrix(BinCoef(Nsteps, rowID) * p^Nsteps)

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
	MOMS = cbind(moms, BS.moments(under, rfr, sigma, yield, maty))
		
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
	nodprob = as.matrix(BinCoef(Nsteps, rowID) * p^rowID * (1-p)^(Nsteps-rowID))
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
	MOMS = cbind(moms, BS.moments(under, rfr, sigma, yield, maty))

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

	disc = exp((rfr - yield)*dt)

	# calculate BS price
	BSres = BS.price(under, strike, rfr, sigma, maty, yield, opt.type=c("call","put"))
	BSp = BSres[,1]
	
	# calculate BS differential factor D1 and D2
	d1 = BSres[,2]
	d2 = BSres[,3]
	
	# probabilities
	p = InvPP(d2, Nsteps)
	ps = InvPP(d1, Nsteps)

	# movement parameters
	# up
	u = disc * ps / p
	# down
	d = disc * (1-ps) / (1-p)

	# calculate share price path
	path = StepMat(under, Nsteps+1, u, d)

	# calculate nodal probabilities
	nodprob = as.matrix(BinCoef(Nsteps, rowID) * p^rowID * (1-p)^(Nsteps-rowID))
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
	MOMS = cbind(moms, BS.moments(under, rfr, sigma, yield, maty))

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
	nodprob = as.matrix(BinCoef(Nsteps, rowID) * p^rowID * (1-p)^(Nsteps-rowID))
	
	edgf = EdgeFact(init, 0, 5.4)
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
	
	nser = nrow(mi)
	
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
		cplot(S, base=t, main="Brownian Motion", ...)
		if(nser == 1){
			if(geom)			
				X = exp(mi*t) else
					X = mi*t
			cplot(X, base=t*T, col="red", append=TRUE, show.legend=FALSE, ...)
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
	
	if(!identical(nrow(mi), nrow(sigma), nrow(S0))){
		cat("'~' Input parameters length differs \n")
		return(NULL)
	}
		
	BM = BroMot(nsim=nsim, T=T, S0=S0, mi=mi, sigma=sigma, geom=geom, same.rnd, plot=FALSE, ...)
	
	if(nser > 2){
		pos = BinCoef(nser, 2)
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
# Probability density function for hittin barriers
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
ProbHit = function(B=0, S0=0, mi, sigma){

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
	cat("Probability to reach Absorption Barrier: \n")
	
	ph; 
	
}
# First hit time
FirstHit = function(B, S0, mi, geom=FALSE, sigma=NULL){

	# check sigma parameters
	if(is.null(sigma) & geom){
		cat("If Geometric is selected sigma must be specified! \n")
		return(NULL)
	}

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
	
	# return expected time
	cat("Expected First Hitting Time: \n")
	
	ht;
	
}
# Edgeworth factor
EdgeFact = function(x, s, k){
	
	1 + s * (x^3-3*x) + (k-3) * (x^4-6*x^2+3)/24 + s^2*(x^6-15*x^4+45*x^2-15)/72
	
	}
# Sample moments
SampMom = function(P, X, moms=1:2){
	
	res = matrix(NA, length(moms), 1)
	rownames(res) = paste("Mom_", moms, sep="")
	j = 1
	# calculate pdf moments
	while(j <= length(moms)){
		res[j, ] = t(X)^moms[j] %*% P
		j = j +1	
	}		
	
	res
}
# Binomial coefficient
BinCoef = function(N, n){
	factorial(N) / (factorial(N-n)*factorial(n))
}
# Peizer-Pratt Inversion formula
InvPP = function(z, n){
	0.5 + sign(z) * sqrt(0.25 - 0.25*exp( -(z/(n+(1/3)+(0.1/(n+1))) )^2 * (n+(1/6))))
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

	mm;
}







