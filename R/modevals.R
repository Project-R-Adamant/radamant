pred_error=function(target, pred, pc=FALSE){
	# trim NAs values
	Logger(message = "trim NAs values", from = "pred_error", line = 2, level = 1);
	if(any(is.na(pred))){
		NAs=is.na(pred)
		pred = pred[!NAs]
		target = target[!NAs]
	}
	# calc error	
	Logger(message = "calc error	", from = "pred_error", line = 8, level = 1);
	res = target - pred
	# results as percentage
	Logger(message = "results as percentage", from = "pred_error", line = 10, level = 1);
	if (pc) 
		res = (res / target) * 100
	# results
	Logger(message = "results", from = "pred_error", line = 13, level = 1);
	res
}
av_er=function(target, pred, pc=FALSE){
	# remove NAs 
	Logger(message = "remove NAs ", from = "av_er", line = 2, level = 1);
	if (any(is.na(pred))) {
			NAs=is.na(pred)
			pred = pred[!NAs]
			target = target[!NAs]
		}
	# average discard
	Logger(message = "average discard", from = "av_er", line = 8, level = 1);
	res = mean( target - pred)
	# results as percentages
	Logger(message = "results as percentages", from = "av_er", line = 10, level = 1);
	if (pc) 
		res = res*100
	# return results
	Logger(message = "return results", from = "av_er", line = 13, level = 1);
	res
}
abs_avdi = function(target, pred, pc=FALSE){
	# remove NAs 
	Logger(message = "remove NAs ", from = "abs_avdi", line = 2, level = 1);
	if(any(is.na(pred))){
		NAs=is.na(pred)
		pred = pred[!NAs]
		target = target[!NAs]
	}
	# average absolute discard
	Logger(message = "average absolute discard", from = "abs_avdi", line = 8, level = 1);
	res = mean( abs( target - pred) )
	# results as percentages
	Logger(message = "results as percentages", from = "abs_avdi", line = 10, level = 1);
	if (pc) 
		res = res *100 
	# return results
	Logger(message = "return results", from = "abs_avdi", line = 13, level = 1);
	res
}
mse = function(target, pred){
	# trim NAs
	Logger(message = "trim NAs", from = "mse", line = 2, level = 1);
	if(any(is.na(pred))){
		NAs=is.na(pred)
		pred = pred[!NAs]
		target = target[!NAs]
	}
	# average squared discard
	Logger(message = "average squared discard", from = "mse", line = 8, level = 1);
	res = mean( (target - pred)^2 )
	# results
	Logger(message = "results", from = "mse", line = 10, level = 1);
	res
}
sde = function(target, pred){
	# trim NAs
	Logger(message = "trim NAs", from = "sde", line = 2, level = 1);
	if (any(is.na(pred))) {
		NAs=is.na(pred)
		pred = pred[!NAs]
		target = target[!NAs]
	}
	# error standard deviation
	Logger(message = "error standard deviation", from = "sde", line = 8, level = 1);
	res = ( mean( (target - pred)^2 ) ) ^ 0.5
	res
}
track_sign = function(target, pred){
	# trim NAs
	Logger(message = "trim NAs", from = "track_sign", line = 2, level = 1);
	if(any(is.na(pred))) {
		NAs=is.na(pred)
		pred = pred[!NAs]
		target = target[!NAs]
	}
	# error track signal
	Logger(message = "error track signal", from = "track_sign", line = 8, level = 1);
	e = target - pred
	trsi = sum(e) / sum(abs(e))
	trsi
}
track_sign_exp = function(target, pred){
	# trim NAs
	Logger(message = "trim NAs", from = "track_sign_exp", line = 2, level = 1);
	if(any(is.na(pred))) {
		NAs=is.na(pred)
		pred = pred[!NAs]
		target = target[!NAs]
	}
	e = target - pred
	d = ema(e, 2)
	g = ema(abs(e), 2)
	# error exponential track signal
	Logger(message = "error exponential track signal", from = "track_sign_exp", line = 11, level = 1);
	etrsi = abs(d/g)
	etrsi
}
