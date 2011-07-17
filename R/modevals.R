pred_error=function(target, pred, pc=FALSE){
	# trim NAs values
	if(any(is.na(pred))){
		NAs=is.na(pred)
		pred = pred[!NAs]
		target = target[!NAs]
	}
	# calc error	
	res = target - pred
	# results as percentage
	if (pc) 
		res = (res / target) * 100
	# results
	res
}
av_er=function(target, pred, pc=FALSE){
	# remove NAs 
	if (any(is.na(pred))) {
			NAs=is.na(pred)
			pred = pred[!NAs]
			target = target[!NAs]
		}
	# average discard
	res = mean( target - pred)
	# results as percentages
	if (pc) 
		res = res*100
	# return results
	res
}
abs_avdi = function(target, pred, pc=FALSE){
	# remove NAs 
	if(any(is.na(pred))){
		NAs=is.na(pred)
		pred = pred[!NAs]
		target = target[!NAs]
	}
	# average absolute discard
	res = mean( abs( target - pred) )
	# results as percentages
	if (pc) 
		res = res *100 
	# return results
	res
}
mse = function(target, pred){
	# trim NAs
	if(any(is.na(pred))){
		NAs=is.na(pred)
		pred = pred[!NAs]
		target = target[!NAs]
	}
	# average squared discard
	res = mean( (target - pred)^2 )
	# results
	res
}
sde = function(target, pred){
	# trim NAs
	if (any(is.na(pred))) {
		NAs=is.na(pred)
		pred = pred[!NAs]
		target = target[!NAs]
	}
	# error standard deviation
	res = ( mean( (target - pred)^2 ) ) ^ 0.5
	res
}
track_sign = function(target, pred){
	# trim NAs
	if(any(is.na(pred))) {
		NAs=is.na(pred)
		pred = pred[!NAs]
		target = target[!NAs]
	}
	# error track signal
	e = target - pred
	trsi = sum(e) / sum(abs(e))
	trsi
}
track_sign_exp = function(target, pred){
	# trim NAs
	if(any(is.na(pred))) {
		NAs=is.na(pred)
		pred = pred[!NAs]
		target = target[!NAs]
	}
	e = target - pred
	d = ema(e, 2)
	g = ema(abs(e), 2)
	# error exponential track signal
	etrsi = abs(d/g)
	etrsi
}
