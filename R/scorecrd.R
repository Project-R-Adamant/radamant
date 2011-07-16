######### transform input data according to weight of evidence #######
input2woe = function(data, nseg, woe, ...){
	if(is.matrix(data))
		numid = apply(data, 2, is.numeric)
	else
		numid = sapply(data, is.numeric)
	# number of input variables 
	Logger(message = "number of input variables ", from = "input2woe", line = 6, level = 1);
	nvars = NCOL(data)
	# identify categorical variables
	Logger(message = "identify categorical variables", from = "input2woe", line = 8, level = 1);
	catid = (1:nvars)[!numid]
	ncat = length(catid)
	## factorise numerical variables (if any)
	Logger(message = "factorise numerical variables (if any)", from = "input2woe", line = 11, level = 1);
	if(length(numid) > 0){
		temp = Factorise(data[, numid], nseg, ...)
		# allocate segment matrix with names
		Logger(message = "allocate segment matrix with names", from = "input2woe", line = 14, level = 1);
		SEGM = matrix(NA, NROW(temp), NCOL(temp)+ifelse(ncat>0, ncat, 0) )
		colnames(SEGM) = character(ncol(SEGM))
		# place numerical variables
		Logger(message = "place numerical variables", from = "input2woe", line = 17, level = 1);
		SEGM[ , 1:NCOL(temp)] = temp
		colnames(SEGM)[1:NCOL(temp)] = colnames(temp)
		# place categorical variables
		Logger(message = "place categorical variables", from = "input2woe", line = 20, level = 1);
		SEGM[ , (NCOL(temp)+1):NCOL(SEGM)] = as.matrix(data[ ,catid])
		colnames(SEGM)[(NCOL(temp)+1):NCOL(SEGM)] = colnames(data)[catid]
	} else {
		# only categorical variables 
		Logger(message = "only categorical variables ", from = "input2woe", line = 24, level = 1);
		SEGM = as.matrix(data)
	}
	# convert fator data to weight of evidence
	Logger(message = "convert fator data to weight of evidence", from = "input2woe", line = 27, level = 1);
	res = .factor2woe(SEGM, woe)
	cleanup("res")
	# return results as invisible object
	Logger(message = "return results as invisible object", from = "input2woe", line = 30, level = 1);
	invisible(res);
}
## Calculate weight of evidence for a matrix with target variable
WeightEvid = function(data, target, nseg, missing=FALSE, na.replace=NULL, ...){
	if(is.matrix(data))
		numid = apply(data, 2, is.numeric)
	else
		numid = sapply(data, is.numeric)
	# number of input variables 
	Logger(message = "number of input variables ", from = "WeightEvid", line = 6, level = 1);
	nvars = NCOL(data)
	# identify categorical variables
	Logger(message = "identify categorical variables", from = "WeightEvid", line = 8, level = 1);
	catid = (1:nvars)[!numid]
	ncat = length(catid)
	## factorise numerical variables (if any)
	Logger(message = "factorise numerical variables (if any)", from = "WeightEvid", line = 11, level = 1);
	if(length(numid) > 0){
		temp = Factorise(data[, numid], nseg, na.replace)
		# allocate segment matrix with names
		Logger(message = "allocate segment matrix with names", from = "WeightEvid", line = 14, level = 1);
		SEGM = matrix(NA, NROW(temp), NCOL(temp)+ifelse(ncat>0, ncat, 0) )
		colnames(SEGM) = character(ncol(SEGM))
		# place numerical variables
		Logger(message = "place numerical variables", from = "WeightEvid", line = 17, level = 1);
		SEGM[ , 1:NCOL(temp)] = temp
		colnames(SEGM)[1:NCOL(temp)] = colnames(temp)
		# place categorical variables
		Logger(message = "place categorical variables", from = "WeightEvid", line = 20, level = 1);
		SEGM[ , (NCOL(temp)+1):NCOL(SEGM)] = as.matrix(data[ ,catid])
		colnames(SEGM)[(NCOL(temp)+1):NCOL(SEGM)] = colnames(data)[catid]
	} else {
		# only categorical variables 
		Logger(message = "only categorical variables ", from = "WeightEvid", line = 24, level = 1);
		SEGM = as.matrix(data)
	}	
	# check the presence of NAs for each segmented variable
	Logger(message = "check the presence of NAs for each segmented variable", from = "WeightEvid", line = 27, level = 1);
	if(missing){
		
		if(!is.null(na.replace))
			idna = which(apply(SEGM , 2, function(x) as.character(na.replace) %in% x ))
		else
			idna = which(apply(SEGM , 2, function(x) NA %in% x ))
		nna = length(idna)
	} else {
		if(!is.null(na.replace))
			idna = which(apply(SEGM , 2, function(x) as.character(na.replace) %in% x ))
		else
			idna = which(apply(SEGM , 2, function(x) NA %in% x ))
		nna = 0
	}
	# allocate matrix for weigth of evidence results
	Logger(message = "allocate matrix for weigth of evidence results", from = "WeightEvid", line = 41, level = 1);
	numclass = sum(apply(SEGM, 2, function(x) length(unique(x))))
	WOE = (matrix(NA, numclass + ifelse(nna > 0, nna, 0) , 12))
	colnames(WOE) = c("Variable", "Segment", "Obs", "PC.Obs", "Good", "PC.Good", "Bad", "Pc.Bad", "Rate", "Weight.Evidence" ,"Info.Value.Within", "Info.Value")
	## overall target information from the data
	Logger(message = "overall target information from the data", from = "WeightEvid", line = 45, level = 1);
	N = NROW(SEGM)
	# number of goods
	Logger(message = "number of goods", from = "WeightEvid", line = 47, level = 1);
	overall.one = sum(target)
	# number of bads
	Logger(message = "number of bads", from = "WeightEvid", line = 49, level = 1);
	overall.zero = N - overall.one
	j = i = 1
	while(i <= ncol(SEGM)){

		tot = aggregate(target, list(SEGM [,i]), length) 
		one = aggregate(target, list(SEGM [,i]), sum)[, 2] 
		zero = tot[,2] - one
		
		labels = as.character(tot[,1])
		if(i %in% idna){
		
		if(!(missing)){ 
			if(!is.null(na.replace)){
			
				pcone = one / sum(one[-length(one)]) 
				pczero = zero / sum(zero[-length(zero)])
				obs = sum(tot[-nrow(tot),2])
			} else if(is.null(na.replace)){
		
				pcone = one / sum(one) 
				pczero = zero / sum(zero)
				obs = sum(tot[,2])	
			}
		} else if((missing)){
			if(is.null(na.replace)){
			
				labels= c(labels,"Missing")
			
				tot = c(tot[,2], N-sum(tot[,2]))
				one = c(one, N-sum(one))
				zero = c(zero, N-sum(zero))
				pcone = one / overall.one
				pczero = zero / overall.zero
				obs = N
			} else if(!is.null(na.replace)){
				labels[length(labels)] = "Missing"
			
				obs = N
				pcone = one / overall.one
				pczero = zero / overall.zero
			}
		}
		
		} else {
			obs = N
			pcone = one / overall.one
  			pczero = zero / overall.zero
		}
		
		pos = j : (j+length(labels)-1)
		
		# variable name
		Logger(message = "variable name", from = "WeightEvid", line = 90, level = 2);
		WOE[pos, 1] = I(rep(colnames(SEGM )[i], length(labels)))
		# segment labels
		Logger(message = "segment labels", from = "WeightEvid", line = 92, level = 2);
		WOE[pos, 2] = I(labels)
		WOE[pos, 3] = tot[,2]
		WOE[pos, 4] = tot[,2] / obs
		WOE[pos, 5] = one
		WOE[pos, 6] = round(pcone, 5)
		WOE[pos, 7] = zero
		WOE[pos, 8] = round(pczero, 5)
 
		WOE[pos, 9] = round(one / obs , 5)
		# calculate weight of evidence
		Logger(message = "calculate weight of evidence", from = "WeightEvid", line = 101, level = 2);
		WOE[pos, 10] = round(log( pcone / pczero ), 5)
			
		j = j + length(one)
		
		i = i + 1	
	}
	
	# calculate Information Value within each segment
	Logger(message = "calculate Information Value within each segment", from = "WeightEvid", line = 106, level = 1);
	WOE[ ,11] = round((as.numeric(WOE[,6]) - as.numeric(WOE[,8])) * as.numeric(WOE[,10]), 5)
	iv = aggregate(as.numeric(WOE[ ,11]), list(WOE[ ,1]), sum)
	for(i in 1:nrow(WOE)){
		id = which(WOE[ ,1] %in% iv[i, 1])
		WOE[id ,12] = rep(iv[i, 2], length(id)) 
	}
	attr(WOE, "Breaks") = attr(SEGM , "Breaks")
	# return results
	Logger(message = "return results", from = "WeightEvid", line = 114, level = 1);
	(WOE)
}
## Factorise variable
Factorise = function(X, nseg, seg.type = c("freq_equal", "width_equal"), na.replace=NULL){
	if(!is.matrix(X))
		X = as.matrix(X)
	BIN = matrix(NA, NROW(X), NCOL(X)*length(nseg))
	colnames(BIN) = character(NCOL(BIN))
	BR = vector("list", NCOL(BIN))
	I = 1
	j = 1 
	while(j <= length(nseg)){
		i = 1	
		while(i <= NCOL(X)){
			switch(match.arg(seg.type),
			width_equal = {
				# calculate breaks in order to create equally width segments
				w = (max(X[ ,i], na.rm=TRUE) - min(X[ ,i], na.rm=TRUE)) / nseg[j]
				breaks = min(X[ ,i]) + w * seq(0, nseg[j], 1)
				},
			freq_equal = {
				# calculate breaks in order to create equally frequent segments
				breaks = quantile(X[,i], probs = seq(0, 1, 1/nseg[j]), na.rm = TRUE, type=8)
				breaks = unique(breaks)
				}
			)
			# fix the extremes of the breaks	
			Logger(message = "fix the extremes of the breaks	", from = "Factorise", line = 22, level = 3);
			if(length(breaks) > 2){
				breaks[1] <- -Inf
				breaks[length(breaks)] <- Inf
			} else {
				
				breaks = c(-Inf, median(X, na.rm=TRUE), Inf)
			}
			# store breaks in a list
			Logger(message = "store breaks in a list", from = "Factorise", line = 29, level = 3);
			BR[[I]] = breaks #paste(breaks,collapse=" ")
			# factorise variable	
			Logger(message = "factorise variable	", from = "Factorise", line = 31, level = 3);
			fv = as.vector(cut(X[,i], breaks, include.lowest = TRUE, ordered=FALSE))
			# check NAs value in the segment and replace with "UKNOWN"
			if(!is.null(na.replace)){
				if(any(is.na(X[,i])))
					fv[(is.na(fv))] = as.character(na.replace)
			}
			BIN[ ,I] = fv
			colnames(BIN)[I] = paste(colnames(X)[i], nseg[j], sep="__") 
			i = i + 1
			I = I + 1
		}
		j = j + 1
	}
	names(BR) = colnames(BIN)
	attr(BIN, "Breaks") = BR
	class(BIN) = "Factorise"
	invisible(BIN);
}
## print method for Factorise
print.Factorise = function(x, ...){
	
	br = lapply(attr(x, "Breaks"), function(x) paste(round(x,3), collapse=" ") )
	nn = sapply(strsplit(names(br), "__"), function(x) x[[1]])
	nb = sapply(strsplit(names(br), "__"), function(x) x[[2]])
	res = matrix(NA, length(unique(nn)), length(unique(nb))+1)
	colnames(res) = c("Variable", paste("SEG_",unique(nb),sep=""))
	res[,1] = unique(nn)
	for(i in 1:nrow(res)){
	
		res[i, -1] = unlist(br[which(!is.na(match(nn, res[i,1])))])
	}	
	cat("Info:", "\n")
	cat("Number of variables:", ncol(res), "\n", sep=" ")
	cat("Number of segments:", unique(nb), "\n", sep=" ")
	cat(rep("=",40), "\n \n", sep="")
	print(as.data.frame(res))
}
## extract specificied break from an object of class Factorise
extrBreak = function(var, Factors){
	
	if(!is.character(var))
		var = as.character(var)
	
	br = lapply(attr(Factors, "Breaks"), function(x) paste(round(x,3), collapse=" ") )
	nn = sapply(strsplit(names(br), "__"), function(x) x[[1]])
	nb = sapply(strsplit(names(br), "__"), function(x) x[[2]])
	
	check = which(!is.na(match(nn, var)))
	
	if(length(check) == 0){
		cat("'~' The variable specified is not in the list! \n")
		cat("Chose one of these: \n", paste(unique(nn), collpse="\n"))
		return(NULL)
	} else {
		
		br[check]	
	
	}
	
}
Score.card = function(X, Y, nseg=2, col.classes=NULL){
	# change format of the variables (if requested)
	Logger(message = "change format of the variables (if requested)", from = "Score.card", line = 2, level = 1);
	if(!is.null(col.classes) & length(col.classes) == NCOL(X)){
		X = reformat(X, col.classes)	
	}
	# calculate Weight of Evidence
	Logger(message = "calculate Weight of Evidence", from = "Score.card", line = 6, level = 1);
	woe = WeightEvid(X, Y, nseg, missing = FALSE, na.replace=NULL)
	# replace original data with weight of evidence
	Logger(message = "replace original data with weight of evidence", from = "Score.card", line = 8, level = 1);
	X = input2woe(X, nseg, woe, na.replace=0)
	# assign names to new data
	Logger(message = "assign names to new data", from = "Score.card", line = 10, level = 1);
	colnames(X) = unique(woe[,1])
	# allocate matrix for overall Scorecard results
	Logger(message = "allocate matrix for overall Scorecard results", from = "Score.card", line = 12, level = 1);
	SC = matrix(NA, nrow(woe), 9)
	SC[, 1:2] = woe[, 1:2]
	SC[, 3] = woe[, 10]
	colnames(SC) = c("Variable", "Segment", "WoE", "Est.Coef", "Wald-Z", "P-Val", "Ods_ratio", "Score", "Round.Score")
	##Run Logistic Regression to calculate the parameter estimates
	Logger(message = "Run Logistic Regression to calculate the parameter estimates", from = "Score.card", line = 17, level = 1);
	model = glm(Y ~., data = as.data.frame(X), family=binomial(link="logit"))
	reg.res = summary(model)$coefficients
	# number of coeffs for numerical variables
	Logger(message = "number of coeffs for numerical variables", from = "Score.card", line = 20, level = 1);
	k = nrow(reg.res)
	i = 2
	while(i <= k){
		idx = which(!is.na(match(SC[ ,1], rownames(reg.res)[i])))
		SC[ ,4][idx] = round(reg.res[i ,1], 5)
		SC[ ,5][idx] = round(reg.res[i ,3], 5)
		SC[ ,6][idx] = round(reg.res[i ,4], 5)
		i = i + 1
	}
	# get ods ratios
	Logger(message = "get ods ratios", from = "Score.card", line = 30, level = 1);
	SC[ ,7] = round(exp(as.numeric(SC[ ,4])), 5)
	# get score
	Logger(message = "get score", from = "Score.card", line = 32, level = 1);
	SC[ ,8] = round( exp(as.numeric(SC[ ,4]) * as.numeric(SC[ ,3]))*100 , 5)
	# round score
	Logger(message = "round score", from = "Score.card", line = 34, level = 1);
	SC[ ,9] = round(as.numeric(SC[ ,8]), 0)
	# matrix of results
	Logger(message = "matrix of results", from = "Score.card", line = 36, level = 1);
	Results = list(Scorecard = as.data.frame(SC), 
			Model = model, 
			WeightOfEvidence = woe);
	# assign class
	Logger(message = "assign class", from = "Score.card", line = 40, level = 1);
	class(Results) = "scorecard"
	# clean memory
	Logger(message = "clean memory", from = "Score.card", line = 42, level = 1);
	cleanup(keep="Results");
	# return results
	Logger(message = "return results", from = "Score.card", line = 44, level = 1);
	Results;
}
print.scorecard = function(x, ...){

	cat("============================","\n", sep="")
	cat("!Score Card results! '_' \n")
	cat("============================","\n\n", sep="")

	cat("==========","\n", sep="")
	cat("Input Info: '_' \n")
	cat("N. observations: ", length(x[[2]]$effects), "\n" )
	cat("N. variables: ", length(unique(x[[1]][,1])), "\n\n" )

	cat("==========","\n", sep="")
	cat("Results: '_' \n\n")
	print(x[[1]])
}
summary.scorecard = function(object, plot=FALSE, ...){
	# confusion matrix
	Logger(message = "confusion matrix", from = "summary.scorecard", line = 2, level = 1);
	cmat = round(confusionM(object, th=0.5)[[1]], 2)
	# accuracies
	Logger(message = "accuracies", from = "summary.scorecard", line = 4, level = 1);
	acc = round(accuracy(object, 0.5), 2)
	# get predicted values
	Logger(message = "get predicted values", from = "summary.scorecard", line = 6, level = 1);
	pred = predict(object)
	# get target
	Logger(message = "get target", from = "summary.scorecard", line = 8, level = 1);
	target = object[[2]]$y
	# ROC info
	Logger(message = "ROC info", from = "summary.scorecard", line = 10, level = 1);
	rocinfo = rbind(SomerD = SomerD(target, pred),t(KendallTau(target, pred)), Gamma = GKgamma(target, pred), Gini = Gini(object))
	# plot results
	Logger(message = "plot results", from = "summary.scorecard", line = 12, level = 1);
	if(plot){
		# 4x4 window of plots
		Logger(message = "4x4 window of plots", from = "summary.scorecard", line = 14, level = 1);
		par(mfrow=c(2,2))
		# lift chart
		Logger(message = "lift chart", from = "summary.scorecard", line = 16, level = 1);
		Lift(object)
		# cumulative gain
		Logger(message = "cumulative gain", from = "summary.scorecard", line = 18, level = 1);
		Gain(object)
		# ROC
		Logger(message = "ROC", from = "summary.scorecard", line = 20, level = 1);
		ROCplot(object)
	}
	# return list of results
	Logger(message = "return list of results", from = "summary.scorecard", line = 23, level = 1);
	list(Results = object[[1]], Confusion_Matrix = cmat, Accuracy = acc, ROC = rocinfo)
}
predict.scorecard = function(object, ...){
	(1/(1+exp(-predict(object[[2]]))))
}
confusionM = function(target, ...) UseMethod("confusionM")
confusionM.default = function(target, pred, th=0.5, ...){
	rec = pred <= th
		pred[rec] = 0
		pred[!rec] = 1
	freq = table(target, pred)
	res = matrix(0, 9, 3)
	colnames(res) = c("Bad", "Good", "Tot")
	rownames(res) = c("Bad", "%", "r%", "c%", "Good", "%", "r%", "c%", "Tot")
	stats = matrix(0, 1, 5)
	colnames(stats) = c("Accuracy", "FPR", "FNR", "Sensitivity", "Specificity")
	if(all(rec)){
		res[c(1,5) ,1] = freq
	} else if(!any(rec)){
		res[c(1,5) ,2] = freq
	} else {
		res[c(1,5) ,c(1,2)] = freq
	}
	res[9 ,] = colSums(res)
	res[ ,3] = rowSums(res)
	res[c(2,6) ,] = res[c(1,5),] / res[9,3] *100
	res[c(3,7) , 1:2] = res[c(1,5), 1:2] / res[c(1,5),3] *100
	res[c(4,8), 1] = c(res[1,1], res[5,1]) / res[9,1] *100
	res[c(4,8), 2] = c(res[1,2], res[5,2]) / res[9,2] *100
	res[c(3,4,7,8),3] = NA
	stats[ ,1] = (1 - (res[5,1] + res[1,2])/res[9,3])
	stats[ ,2] = res[1,2] / (res[1,2] + res[1,1])
	stats[ ,3] = res[5,1] / (res[5,1] + res[5,2])
	stats[ ,4] = 1 - stats[ ,2]
	stats[ ,5] = 1 - stats[ ,3]
	#cat("Observed * Predicted \n \n")
	list(Confusion.Matrix = res,
		Accuracy = stats*100)
}
confusionM.scorecard = function(target, th=0.5, ...){
	# get predicted values
	Logger(message = "get predicted values", from = "confusionM.scorecard", line = 2, level = 1);
	pred = predict(target)
	# calculate confusion matrix
	Logger(message = "calculate confusion matrix", from = "confusionM.scorecard", line = 4, level = 1);
	confusionM.default(target=target[[2]]$y, pred=pred , th)
}
accuracy = function(x, ...) UseMethod("accuracy")
accuracy.scorecard = function(x, th=0.5, ...){
	# calculate confusion matrix
	Logger(message = "calculate confusion matrix", from = "accuracy.scorecard", line = 2, level = 1);
	confusionM(x , th)$Accuracy
	
}
CalcPairs = function(target, pred, segm_fact=0.002){
	# define variable segmentation
	Logger(message = "define variable segmentation", from = "CalcPairs", line = 2, level = 1);
	seg_prob = trunc(pred / segm_fact)
	# calculate cross frequencies
	Logger(message = "calculate cross frequencies", from = "CalcPairs", line = 4, level = 1);
	tab = table(target, seg_prob)
	# get freq table dimensions
	Logger(message = "get freq table dimensions", from = "CalcPairs", line = 6, level = 1);
	I = nrow(tab)
	J = ncol(tab)
	# concordant pairs
	Logger(message = "concordant pairs", from = "CalcPairs", line = 9, level = 1);
	conc = rep(0, J-1)
	for(i in 1:(I-1)){
		temp = rep(0, J-1)
		j = 1
		while(j < J){
			temp[j] = tab[i ,j] * sum( tab[-(1:i), -(1:j)] )
			j = j + 1
		}
		
		conc[i] = sum(temp)
		
	}
	# discordant pairs
	Logger(message = "discordant pairs", from = "CalcPairs", line = 20, level = 1);
	disc = rep(0, I-1)
	for(i in 1:(I-1)){
		temp = rep(0, J-1)
		j = J
		while(j > 1){
			temp[J-j+1] = tab[i ,j] * sum( tab[-(1:i), -(j:J)] )
			j = j - 1
		}
		disc[i] = sum(temp) 
	}
	# tied on X
	Logger(message = "tied on X", from = "CalcPairs", line = 31, level = 1);
	Tx = rep(0,J)
	for(j in 1:J){
		temp = rep(0, I)
		i = 1
		while(i < I){
			temp[i] = tab[i, j] * sum(tab[(i+1):I, j]) 	
			i = i + 1
		}
		Tx[j] = sum(temp)
	}
	# tied on Y
	Logger(message = "tied on Y", from = "CalcPairs", line = 42, level = 1);
	Ty = rep(0, I)
	for(i in 1:I){
		temp = rep(0, J-1)
		j = 2
		while(j <= J){
			temp[j-1] = tab[i, j-1] * sum(tab[i, j:J]) 	
			j = j + 1
		}
		Ty[i] = sum(temp)
	}
	# get overall values for the pairs
	Logger(message = "get overall values for the pairs", from = "CalcPairs", line = 53, level = 1);
	P = sum(conc)
	Q = sum(disc)
	Ty = sum(Ty)
	Tx = sum(Tx)
	res = cbind(Conc = P, Disc = Q, TiedX = Tx, TiedY = Ty )
	# return results
	Logger(message = "return results", from = "CalcPairs", line = 59, level = 1);
	attr(res, "info_tab") = c(nrow(tab), ncol(tab), sum(tab))
	res
}
SomerD = function(target, pred, ...){
	# calculate pairs
	Logger(message = "calculate pairs", from = "SomerD", line = 2, level = 1);
	pair = CalcPairs(target, pred, ...)
	# calculate index
	Logger(message = "calculate index", from = "SomerD", line = 4, level = 1);
	res = (pair[1] - pair[2]) / (pair[1] + pair[2] + pair[3])
	res
}
KendallTau = function(target, pred, ...){
	# calculate pairs
	Logger(message = "calculate pairs", from = "KendallTau", line = 2, level = 1);
	pair = CalcPairs(target, pred, ...)
	
	info = attr(pair, "info_tab")
	# calculate index
	Logger(message = "calculate index", from = "KendallTau", line = 5, level = 1);
	tb = (pair[1] - pair[2]) / sqrt((pair[1] + pair[2] + pair[4]) * (pair[1] + pair[2] + pair[3]))
	tc = (pair[1] - pair[2]) / ( 0.5*info[3]^2*(1-(1/min( c(info[1], info[2])))) )
	res = cbind(TauB = tb, TauC = tc)
	res
}
GKgamma = function(target, pred, ...){
	# calculate pairs
	Logger(message = "calculate pairs", from = "GKgamma", line = 2, level = 1);
	pair = CalcPairs(target, pred, ...)
	# calculate index
	Logger(message = "calculate index", from = "GKgamma", line = 4, level = 1);
	res = (pair[1] - pair[2]) / (pair[1] + pair[2])
	res
}
ROCplot = function(x, ...) UseMethod("ROCplot")
ROCplot.scorecard = function(x, ...){
	TH = seq(0, 1, length.out=100)
	tab = matrix(0, 2, length(TH))
	for(i in 1:length(TH)){
		acc = accuracy(x, th=TH[i])
		tab[1,i] = acc[ ,4]
		tab[2,i] = acc[ ,5]
	}
	tab = tab / 100
	if(dev.cur() == 1) 
		par(mfrow=c(2,1)) else dev.set(dev.cur())
	
	cplot(X = tab[1,], base = 1-tab[2,], lwd = 1, type="o", col="red", xtitle="1-Specificity", ytitle="Sensitivity", main="ROC curve", legend=c("ROC","Random"), legend.col=c("red", "green"))
	abline(a=0,b=1, col="green", lty=4)
	cplot(cbind(tab[1,],tab[2,]),type="o", lwd=1, col=c("green", "yellow"),
	xtitle="%Pop", ytitle="Sens / Spec", main="Sensitivity VS Specificity",
	legend=c("Sensitivity","Specificity"))
}
Lift = function(x, ...) UseMethod("Lift")
Gain = function(x, ...) UseMethod("Gain")
Lift.scorecard = function(x, pc=0.1, ...){
	pred = predict(x)
	ordid = order(pred, decreasing=TRUE)
	ordtg = x[[2]]$y[ordid]
	n = length(ordid)
	TG = sum(x[[2]]$y)
	TPG = sum(x[[2]]$y) / n
	probs=seq(0, 1, pc)
	segs = quantile(1:n, probs=probs)
	lift = rep(0, length(segs))
	i = 1
	while(i <= length(segs)){
		lift[i] = (sum(ordtg[1:segs[i]]) / length(ordtg[1:segs[i]])) / TPG
		i = i + 1
	}
	cplot(lift[-1], base=probs[-1]*100, col="red", lwd=2, main="Lift chart", ytitle="Lift", xtitle="Pop_segments", legend="Lift")
	abline(h=1, col="green", lty=4)
	invisible(lift)
}
Gain.scorecard = function(x, pc=0.1, ...){
	pred = predict(x)
	ordid = order(pred, decreasing=TRUE)
	ordtg = x[[2]]$y[ordid]
	n = length(ordid)
	TG = sum(x[[2]]$y)
	TPG = sum(x[[2]]$y) / n
	probs=seq(0, 1, pc)
	segs = quantile(1:n, probs=probs)
	gain = rep(0, length(segs))
	i = 2
	while(i <= length(segs)){
		gain[i] = (sum(ordtg[1:segs[i]] / TG)) * 100
		i = i + 1
	}
	cplot(gain, base=probs*100, col="red", lwd=2, main="Cumulated Gain", ytitle="% Cumul. Good", xtitle="Pop_segments", , legend="Gain")
	abline(a=0, b=1, col="green", lty=4)
	invisible(gain)
}
Gini = function(x, ...) UseMethod("Gini")
Gini.scorecard = function(x, glob=TRUE, ...){
	# extract table with needed information for Gini calculation
	Logger(message = "extract table with needed information for Gini calculation", from = "Gini.scorecard", line = 2, level = 1);
	tab = apply(x[[3]][,-(1:2)], 2, as.numeric)
	# calculate Gini for the entire model
	Logger(message = "calculate Gini for the entire model", from = "Gini.scorecard", line = 4, level = 1);
	if(glob){
		gini = Gini.default(x=tab[ ,c(3,5)])
		return(gini)
	} else {
		## calculate Gini for each variable separately
		Logger(message = "calculate Gini for each variable separately", from = "Gini.scorecard", line = 9, level = 1);
		# get variable group
		Logger(message = "get variable group", from = "Gini.scorecard", line = 10, level = 1);
		group = x[[3]][,1:2] 
		gid = unique(group[, 1])
		# allocate matrix of results
		Logger(message = "allocate matrix of results", from = "Gini.scorecard", line = 13, level = 1);
		gini = matrix(NA, length(gid), 1)
		dimnames(gini) = list(gid,"Gini")
		# loop through each variable and calculate corresponding Gini
		Logger(message = "loop through each variable and calculate corresponding Gini", from = "Gini.scorecard", line = 16, level = 1);
		for(i in 1:length(gid)){
			lab = which(!is.na(match(group[,1], gid[i])))
			gini[i] = Gini.default(x=tab[lab ,c(3,5)])
		}
		cleanup(keep = "gini")	
		return(gini)
	}
}
Gini.default = function(x, ...){
	# calculate "rate" if needed
	if(ncol(x) == 2)
		x = cbind(x, apply(x, 1, function(x) x[1]/sum(x)) ) 
	# calculate Gini for the entire model based on the info table provided
	Logger(message = "calculate Gini for the entire model based on the info table provided", from = "Gini.default", line = 5, level = 1);
	x = x[order(x[,3]),]
	cumgood = cumsum(x[,1]/sum(x[,1]))
	temp =  c(0, cumgood[-length(cumgood)])
	# Gini formula
	Logger(message = "Gini formula", from = "Gini.default", line = 9, level = 1);
	gini = (2*(1 - (sum((x[,2]/sum(x[,2]))*((cumgood + temp)/2))))) - 1
	cleanup(keep = "gini")	
	gini
}


