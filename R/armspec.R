Arma.Spec = function(X, ar_ord = 1, ma_ord = 1, vfreq=NULL){
	N = length(X)
	sigma = var(X,na.rm=TRUE)
	# check and assign vector of frequencies
	Logger(message = "check and assign vector of frequencies", from = "Arma.Spec", line = 4, level = 1);
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
	Logger(message = "get arima coefficients", from = "Arma.Spec", line = 15, level = 1);
	arma = arima(X,c(ar_ord,0,ma_ord))
	# ma part
	Logger(message = "ma part", from = "Arma.Spec", line = 17, level = 1);
	theta = arma$coef[ (ar_ord+1):(ar_ord+ma_ord)]
	# ar part
	Logger(message = "ar part", from = "Arma.Spec", line = 19, level = 1);
	phi = arma$coef[1:ar_ord]
	# MA spectrum
	Logger(message = "MA spectrum", from = "Arma.Spec", line = 21, level = 1);
	if(ma_ord >= 1){
		spec_ma = matrix(0, length(lambda), ma_ord) 
		for(i in 1: ma_ord)
			spec_ma[,i] = (1 + theta[i]^2 + 2*theta[i]*cos(1*lambda))
		spec_ma = apply(spec_ma,1,prod)
		} else {
		spec_ma = 1
	}
	# AR spectrum
	Logger(message = "AR spectrum", from = "Arma.Spec", line = 30, level = 1);
	if(ar_ord >= 1){
		spec_ar = matrix(0, length(lambda), ar_ord) 
		for(i in 1:ar_ord)
			spec_ar[,i] = (1 + phi[i]^2 + 2*phi[i]*cos(1*lambda))
		spec_ar = apply(spec_ar,1,prod)
		} else {
		spec_ar = 1	
	}
	# spectrum calculation
	Logger(message = "spectrum calculation", from = "Arma.Spec", line = 39, level = 1);
	specx = (sigma / 2*pi) * (spec_ma / spec_ar)
	# plot spectogram
	Logger(message = "plot spectogram", from = "Arma.Spec", line = 41, level = 1);
	par(mfrow=c(2,1))
	plot(specx, type="l", col="blue", cex.axis=0.7, main="Spectrum")
	grid()
	# plot original series
	Logger(message = "plot original series", from = "Arma.Spec", line = 45, level = 1);
	plot(X, type="l", cex.axis=0.7, main="Original series")
	grid()
}
