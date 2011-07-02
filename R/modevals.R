pred_error=function(actual, predicted, pc=FALSE)
{
	# trim NAs values
	if(any(is.na(predicted))){
			NAs=is.na(predicted)
			predicted = predicted[!NAs]
			actual = actual[!NAs]
	}
	
	# calc error	
	res = actual - predicted
	
	# results as percentage
	if (pc) res = (res / actual) * 100
	
	# results
	res
}
av_er=function(actual, predicted, pc=FALSE)
{
	# remove NAs 
	if (any(is.na(predicted))) {
		
			NAs=is.na(predicted)
			predicted = predicted[!NAs]
			actual = actual[!NAs]
		}
		
	# average discard
	res = mean( actual - predicted)
	
	# results as percentages
	if (pc) res = res *100
	# return results
	res
}
abs_avdi = function(actual, predicted, pc=FALSE)
{
	# remove NAs 
	if(any(is.na(predicted))){
			NAs=is.na(predicted)
			predicted = predicted[!NAs]
			actual = actual[!NAs]
	}
	
	# average absolute discard
	res = mean( abs( actual - predicted) )
	
	# results as percentages
	if (pc) res = res *100 
	# return results
	res
}
mse = function(actual, predicted)
{
	# trim NAs
	if(any(is.na(predicted))){
		NAs=is.na(predicted)
		predicted = predicted[!NAs]
		actual = actual[!NAs]
	}
	
	# average squared discard
	res = mean( (actual - predicted)^2 )
	# results
	res
}
sde = function(actual, predicted)
{
	# trim NAs
	if (any(is.na(predicted))) 
		{
			NAs=is.na(predicted)
			predicted = predicted[!NAs]
			actual = actual[!NAs]
		}
		
	res =  ( mean( (actual - predicted)^2 ) ) ^ 0.5
	res
}
track_sign = function(actual, predicted)
{
	# trim NAs
	if(any(is.na(predicted))) {
			NAs=is.na(predicted)
			predicted = predicted[!NAs]
			actual = actual[!NAs]
		}
		
	e = actual - predicted
	
	trsi = sum(e) / sum(abs(e))
	trsi
}
track_sign_exp = function(actual, predicted, beta=0.1)
{
	# trim NAs
	if (any(is.na(predicted))) 
		{
			NAs=is.na(predicted)
			predicted = predicted[!NAs]
			actual = actual[!NAs]
		}
		
	e = actual - predicted
	
	d=ema(e, 2, TRUE)
	g=ema(abs(e), 2, TRUE)
	
	etrsi = abs(d/g)
	etrsi
}
