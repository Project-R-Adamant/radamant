Arma.Spec = function(X, ar_ord = 1, ma_ord = 1, vfreq=NULL){
	
	N = length(X)
	sigma = var(X,na.rm=TRUE)
	
	# check and assign vector of frequencies
	if(is.null(vfreq)){
		
		lambda = 1:ceiling(N/2)
		
		} else {
			
			if(length(vfreq) < 2){
				
				cat("Frequency vector too small, set default. \n")
				lambda = 1:ceiling(N/2)
			
			} else {
			
			lambda = vfreq
			
			}	
		}
	
	# get arima coefficients
	arma = arima(X,c(ar_ord,0,ma_ord))
	
	# ma part
	theta = arma$coef[ (ar_ord+1):(ar_ord+ma_ord)]
	# ar part
	phi = arma$coef[1:ar_ord]
	
	# MA spectrum
	if(ma_ord >= 1){
		spec_ma = matrix(0, length(lambda), ma_ord) 
		for(i in 1: ma_ord)
			spec_ma[,i] = (1 + theta[i]^2 + 2*theta[i]*cos(1*lambda))
			
		spec_ma = apply(spec_ma,1,prod)
		} else {
		spec_ma = 1
	}
	
	# AR spectrum
	if(ar_ord >= 1){
		spec_ar = matrix(0, length(lambda), ar_ord) 
		for(i in 1:ar_ord)
			spec_ar[,i] = (1 + phi[i]^2 + 2*phi[i]*cos(1*lambda))
		spec_ar = apply(spec_ar,1,prod)
		} else {
		spec_ar = 1	
	}
	
	# spectrum calculation
	specx = (sigma / 2*pi) * (spec_ma / spec_ar)
	
	# plot spectogram
	par(mfrow=c(2,1))
	plot(specx, type="l", col="blue", cex.axis=0.7, main="Spectrum")
	grid()
	# plot original series
	plot(X, type="l", cex.axis=0.7, main="Original series")
	grid()
	
}
