

######### transform input data according to weight of evidence #######
input2woe = function(data, nseg, woe, na.replace, seg.type=c("freq_equal", "width_equal")){
	
	#if(class(woe) != "WoE"){
	#	cat("'~' Object must be of class 'WoE'")
	#	return(NULL)
	#}
	
	seg.type = match.arg(seg.type)
	# factorise data
	SEGM = Factorise(data, nseg, na.replace, seg.type=seg.type)
	# convert fator data to weight of evidence
	res = factor2woe(SEGM, woe)
	
	# return results as invisible object
	invisible(res)
	
}


######### transform factorise data to weight of evidence #######
factor2woe = function(segm, woe){
	
	RES = matrix(0, NROW(segm), NCOL(segm))
	colnames(RES) = colnames(segm)
	var = 1
	while(var <= NCOL(segm)){

		idx = which(!is.na(match(woe[ ,1], colnames(segm)[var])))

		for(i in (idx))
			RES[segm[ ,var] == woe[i ,2], var] = as.numeric(woe[i, 10])
		
		var = var + 1

	}
	
	invisible(RES);
}
######################################################################


## Calculate weight of evidence for a matrix with target variable
WeightEvid = function(data, target, nseg, missing=FALSE, na.replace=NULL){
	
	if(!is.matrix(data))
		data = as.matrix(data)

	# number of input variables 
	nvars = NCOL(data)
	
	## factorise variable
	SEGM = Factorise(data, nseg, na.replace)
	
	# check the presence of NAs for each segmented variable
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
	WOE = (matrix(NA, sum(nvars * nseg) + ifelse(nna > 0, nna, 0) , 12))
	colnames(WOE) = c("Variable", "Segment", "Obs", "PC.Obs", "Good", "PC.Good", "Zero", "Pc.Zero", "Rate", "Weight.Evidence" ,"Info.Value.Within", "Info.Value")
	
	## overall target information from the data
	N = NROW(data)
	# number of goods
	overall.one = sum(target)
	# number of bads
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
		WOE[pos, 1] = I(rep(colnames(SEGM )[i], length(labels)))
		# segment labels
		WOE[pos, 2] = I(labels)

		WOE[pos, 3] = tot[,2]
		WOE[pos, 4] = tot[,2] / obs

		WOE[pos, 5] = one
		WOE[pos, 6] = round(pcone, 5)

		WOE[pos, 7] = zero
		WOE[pos, 8] = round(pczero, 5)
 
		WOE[pos, 9] = round(one / obs , 5)


		# calculate weight of evidence
		WOE[pos, 10] = round(log( pcone / pczero ), 5)
			
		j = j + length(one)
		
		i = i + 1	
	}
	
	# calculate Information Value within each segment
	WOE[ ,11] = round((as.numeric(WOE[,6]) - as.numeric(WOE[,8])) * as.numeric(WOE[,10]), 5)

	iv = aggregate(as.numeric(WOE[ ,11]), list(WOE[ ,1]), sum)

	for(i in 1:nrow(WOE)){
		id = which(WOE[ ,1] %in% iv[i, 1])
		WOE[id ,12] = rep(iv[i, 2], length(id)) 
	}

	attr(WOE, "Breaks") = attr(SEGM , "Breaks")

	# return results
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
			if(length(breaks) > 2){
				
				breaks[1] <- -Inf
				breaks[length(breaks)] <- Inf
					
			} else {
				
				breaks = c(-Inf, median(X, na.rm=TRUE), Inf)
					
			}
			
			# store breaks in a list
			BR[[I]] = breaks #paste(breaks,collapse=" ")
		
			# factorise variable	
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
	#class(BIN) = "Factorise"
	
	invisible(BIN);

}

#########################################################

## print method for Factorise
print.Factorise = function(x, ...){
	
	br = lapply(attr(x, "Breaks"), function(x) paste(round(x,3), collapse=" ") )
	nn = sapply(strsplit(names(br), "__"), function(x) x[[1]])
	nb = sapply(strsplit(names(br), "__"), function(x) x[[2]])

	res = matrix(NA, length(unique(nn)), length(unique(nb))+1)
	colnames(res) = c("Variable", paste("SEG_",unique(nb),sep=""))

	res[,1] = unique(nn)

	for(i in 1:ncol(res)){
	
		res[i, -1] = unlist(br[which(!is.na(match(nn, res[i,1])))])
	}	

	cat("Info:", "\n")
	cat("Number of variables:", ncol(res), "\n", sep=" ")
	cat("Number of segments:", unique(nb), "\n", sep=" ")
	cat(rep("=",40), "\n \n", sep="")
	print(as.data.frame(res))

}


## extract specificied break from an object of class Factorise
extrBreak = function(var){
	
	if(!is.character(var))
		var = as.character(var)
	
	br = lapply(attr(x, "Breaks"), function(x) paste(round(x,3), collapse=" ") )
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





Score.card = function(Y, X, woe){

	SC = matrix(NA, nrow(woe), 8)
	colnames(SC) = c("Variable", "Segment", "WoE", "Est.Coef", "Wald-Z", "P-Val", "Score", "Round.Score")

	##Run Logistic Regression to calculate the parameter estimates
	model <- glm(Y ~., data = as.data.frame(X), family=binomial(link="logit"))

	SC[ ,1:2] = (woe[,1:2])
	SC[ ,3] = woe[ ,10]

	reg.res = summary(model)$coefficients

	# number of coeffs
	k = nrow(reg.res)

	i = 2
	while(i <= nrow(reg.res)){

		idx = which(!is.na(match(SC[ ,1], rownames(reg.res)[i])))
		SC[idx ,4] = round(reg.res[i ,2], 5)
		SC[idx ,5] = round(reg.res[i ,3], 5)	
		SC[idx ,6] = round(reg.res[i ,4], 5)
		i = i + 1
	
	}

	SC[ ,7] = round( -(as.numeric(SC[ ,3]) * as.numeric(SC[ ,4])+reg.res[1 ,1]/k) * 100, 5)
	
	SC[ ,8] = round(as.numeric(SC[ ,7]), 0)

	Results = list(Scorecard = as.data.frame(SC), 
			Model = model, 
			WeightOfEvidence = woe);
	
	class(Results) = "scorecard"

	Results;

}


summary.scorecard = function(sc, plot=FALSE, ...){

	# confusion matrix
	cmat = round(confusionM(sc, th=0.5)[[1]], 2)

	# accuracies
	acc = round(accuracy(sc, 0.5), 2)
	
	pred = (1/(1+exp(-predict(sc[[2]]))))
	target = sc[[2]]$y

	# ROC info
	rocinfo = rbind(SomerD = SomerD(target, pred),t(KendallTau(target, pred)),Gamma = 	GKgamma(target, pred), Gini = Gini(sc))
	
	if(plot){
		par(mfrow=c(2,2))
		Lift(sc)
		Gain(sc)
		ROCplot(sc)
	}
	
	list(Results = sc[[1]], Confusion_Matrix = cmat, Accuracy = acc, ROC = rocinfo)
}

#############################################################################
confusionM = function(x, ...) UseMethod("confusionM")

confusionM.default = function(orig, pred, th){

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


confusionM.scorecard = function(x, th){

	pred = (1/(1+exp(-predict(x[[2]]))))

	confusionM.default(orig=x[[2]]$y, pred=pred , th)
}

accuracy = function(x, ...) UseMethod("accuracy")

accuracy.scorecard = function(x, th){

	#cat("Accuracy measures: \n")
	confusionM(x , th)$Accuracy
	
}


#############################################################################

CalcPairs = function(target, pred, segm_fact=0.002){
	
	seg_prob = trunc(pred / segm_fact)
	tab = table(target, seg_prob)

	I = nrow(tab)
	J = ncol(tab)

	# concordant pairs
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

	# tied X
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


	# tied Y
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


	P = sum(conc)
	Q = sum(disc)
	Ty = sum(Ty)
	Tx = sum(Tx)

	res = cbind(Conc = P, Disc = Q, TiedX = Tx, TiedY = Ty )
	
	attr(res, "info_tab") = c(nrow(tab), ncol(tab), sum(tab))

	res
}



#############################################################################
SomerD = function(target, pred, ...){

	pair = CalcPairs(target, pred, ...)
	
	res = (pair[1] - pair[2]) / (pair[1] + pair[2] + pair[4])

	res

}

KendallTau = function(target, pred, ...){

	pair = CalcPairs(target, pred, ...)
	
	info = attr(pair, "info_tab")

	tb = (pair[1] - pair[2]) / sqrt((pair[1] + pair[2] + pair[4]) * (pair[1] + pair[2] + pair[3]))
	tc = (pair[1] - pair[2]) / ( 0.5*info[3]^2*(1-(1/min( c(info[1], info[2])))) )

	res = cbind(TauB = tb, TauC = tc)

	res

}

GKgamma = function(target, pred, ...){

	pair = CalcPairs(target, pred, ...)
	
	res = (pair[1] - pair[2]) / (pair[1] + pair[2])

	res
}


##############################################################

ROCplot = function(x, ...) UseMethod("ROCplot")

ROCplot.scorecard = function(sc){

	TH = seq(0, 1, length.out=100)
	tab = matrix(0, 2, length(TH))

	for(i in 1:length(TH)){
		acc = accuracy(sc, th=TH[i])
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



#############################################

Lift = function(x, ...) UseMethod("Lift")
Gain = function(x, ...) UseMethod("Gain")

Lift.scorecard = function(sc, pc=0.1, ...){
	
	pred = (1/(1+exp(-predict(sc[[2]]))))

	ordid = order(pred, decreasing=TRUE)
	ordtg = sc[[2]]$y[ordid]

	n = length(ordid)
	
	TG = sum(sc[[2]]$y)
	TPG = sum(sc[[2]]$y) / n

	probs=seq(0, 1, pc)
	segs = quantile(1:n, probs=probs)
	
	lift = rep(0, length(segs))
	i = 1
	while(i <= length(segs)){

		lift[i] = (sum(ordtg[1:segs[i]]) / length(ordtg[1:segs[i]])) / TPG
		i = i + 1
	}
	
	cplot(lift, base=probs*100, col="red", lwd=2, main="Lift chart", ytitle="Lift", xtitle="Pop_segments", legend="Lift")
	abline(h=1, col="green", lty=4)
	
	invisible(lift)

}


Gain.scorecard = function(sc, pc=0.1, ...){
	
	pred = (1/(1+exp(-predict(sc[[2]]))))

	ordid = order(pred, decreasing=TRUE)
	ordtg = sc[[2]]$y[ordid]
	
	n = length(ordid)
		
	TG = sum(sc[[2]]$y)
	TPG = sum(sc[[2]]$y) / n
	
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


############################

Gini = function(x, ...) UseMethod("Gini")

Gini.scorecard = function(sc, glob=TRUE, ...){
	
	tab = apply(sc[[3]][,-(1:2)], 2, as.numeric)

	if(glob){
	
		tab = tab[order(tab[,7]),]
		cumgood <- cumsum(tab[,3]/sum(tab[,3]))
		temp =  c(0, cumgood[-length(cumgood)])
	
		gini = (2*(1 - (sum((tab[,5]/sum(tab[,5]))*((cumgood + temp)/2))))) - 1
	
		return(gini)
		
	} else {

		group = sc[[3]][,1:2] 
		gid = unique(group[, 1])
	
		gini = matrix(NA, length(gid), 1)
		dimnames(gini) = list(gid, "Gini")

		for(i in 1:length(gid)){

			lab = which(!is.na(match(group[,1], gid[i])))
		
			temptab = tab[lab, ][order(tab[lab,7]),]

			cumgood = cumsum(temptab[ ,3]/sum(temptab[ ,3]))
			temp =  c(0, cumgood[-length(cumgood)])
			gini[i] = (2*(1 - (sum((temptab[ ,5]/sum(temptab[ ,5]))*((cumgood + temp)/2))))) - 1
		
		}
	
		return(gini)
	}
	
}
