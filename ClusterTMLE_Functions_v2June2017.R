#======================================================================================
# R functions to implement the hierarchical TMLEs 
#		described in Balzer et al: "A New Approach to Hierarchical Data Analysis: 
#		Targeted Maximum Likelihood Estimation of Cluster-Based Effects Under Interference
# 
# Programer: Laura Balzer
# Contact: lbbalzer@gmail.com
#
# Last update: June 2, 2017
#
#
# This code is focused on (cluster-level) population average treatment effect (PATE):
#		E[ Yc(1)- Yc(0)]
#		expected difference in the counterfactual cluster-level outcomes Yc(a)
#		
# When there is informative cluster size, this causal parameter will generally not equal
#		the pooled individual-level causal effect E[Y(1)-Y(0)]
#		where Y(a) denotes the individual-level counterfactual outcome
#
# Disclaimer: Without loss of generality, this code assumes both the cluster-level outcome Yc 
#		and the individual-level outcome Y are in [0,1] (e.g. binary or a proportion)
#		If not, the outcome should be rescaled as outlined in Gruber and van der Laan 2010
#		Rescaling will improve robustness and not introduce bias. 
#		(i.e. this is a good thing)
#
# Disclaimer II: The cluster-level outcome is Yc_j = 1/n_j  sum_i^{n_j} Y_{ij}
#		Corresponding to the weights alpha_{ij}= 1/n_j
#		No other weighting scheme is currently implemented
#===========================================================================================


#-----------------------------------------------------#-----------------------------------------------------
# do.estimation.inference: function to implement Hierarchical TMLE
#	 with SuperLearner, using Adaptive Prespecification, or parametric regression
#
# input: 
#		psi (true value of PATE), 
#		observed data (must have a column id indicating cluster membership), 
#		individual-level adjustment variables (ind.adj), cluster-level adjustment variables (clust.adj),
#		indicator if initial estimation of the conditional mean outcome is at the individual-level (Qinit.Indv),
#		adjustment variables for the conditional mean outcome (QAdj), 
#		indicator if the working model is considered to hold (i.e. estimation under smaller sub-model),
#		adjustment variables for the propensity score (gAdj),
#		marginal probability of the treatment (prob.txt),
#		indicator to run the full SuperLearner,
#		indicator to do Adaptive Prespecification (only appropriate for trials),
#		SuperLearner libraries for the conditional mean outcome and pscore (SL.library.Q, SL.library.g),
#		indicator to print updates (verbose)
#
# output: estimation and inference
#-----------------------------------------------------#-----------------------------------------------------

do.estimation.inference<- function(psi=NA, data, ind.adj=NULL, clust.adj=NULL, 
	Qinit.Indv=NA, QAdj=NULL,  work.model=F,  gAdj=NULL, prob.txt=0.5,
	Do.SuperLearner=F, Do.AdaptivePrespec=F, SL.library.Q, SL.library.g,  
	verbose=T){	

	# specify the number of independent units (clusters)
	n.indpt<- length(unique(data$id)) 

	if(Do.SuperLearner){		
		# If running the full SuperLearner algorithm  
		est <-  suppressWarnings( doTMLE(train=data, 
						Qinit.Indv=Qinit.Indv, QAdj=QAdj, work.model=work.model, 
						gAdj=gAdj,  prob.txt=prob.txt,  
						Do.SuperLearner =T, SL.library.Q=SL.library.Q, SL.library.g=SL.library.g, 
						verbose=verbose) )
		QAdj<- gAdj<- 'SL'
		DoCVinf<- F
	
	} else if (Do.AdaptivePrespec){
		 # If instead doing adaptive pre-specification to minimize variance and maximize power in a trial
		 # Reference: Balzer et al. 2016 Stat Med. 
		select<- do.adaptive.prespec(data=data, ind.adj=ind.adj, clust.adj=clust.adj, 
						work.model=work.model,  prob.txt=prob.txt)
		
		# After selecting the candidate TMLE (corresponding to adjustment variables for the outcome regression
		#	and propensity score), run the TMLE algorithm
		Qinit.Indv=select$Qinit.Indv
		QAdj = select$QAdj
		gAdj = select$gAdj
		est <-  doTMLE(train=data, Qinit.Indv=Qinit.Indv, QAdj=QAdj, work.model=work.model, 
						gAdj=gAdj,  prob.txt=prob.txt, verbose=verbose) 
		# Requires a cross-validated variance estimate. See Balzer et al. for further details
		DoCVinf<- T
	
	}else{
		# if estimating via a priori specified regression (i.e not doing data-adaptive estimation)
		est<- doTMLE( train=data, Qinit.Indv= Qinit.Indv, QAdj=QAdj,  work.model=work.model,  
						 gAdj=gAdj, prob.txt=prob.txt, verbose=verbose)
		QAdj <- gAdj<- 'glm'
		DoCVinf<- F	
	}
		
	
	# to get inference we need to respect the independent unit as the cluster
	if(!DoCVinf){ 
		
		# if not obtaining a cross-validated variance estimate	
		if(work.model) {
			# if estimating under the submodel, calculate the individual-level influence-curve
			IC.ind<- est$H.AW*(data$Y-est$QbarAW) + est$Qbar1W - est$Qbar0W - est$ATE
			# average to the cluster-level
			IC<- aggregate(data.frame(IC.ind), by=list(id=data$id), mean)[,2]	
		
		}else{
			# if estimating under the general (larger) model,  calculate the cluster-level influence curve
			Yc <- aggregate(data.frame(data$Y), by=list(data$id), mean)[,2]
			IC<- est$H.AW*(Yc-est$QbarAW) + est$Qbar1W - est$Qbar0W - est$ATE
		}
		
		# sanity check that solving EffIC
		if(verbose) print(paste('proof solve EffIC:', mean(IC)))
	
	} else{		 
		# if obtaining a cross-validated variance estimate
		IC<- getIC.CV(data=data, work.model=work.model, Qinit.Indv=Qinit.Indv, QAdj=QAdj,  
			gAdj=gAdj, prob.txt=prob.txt)$IC
		# already aggregated to the cluster-level
	} 
	
	var.IC<- var(IC)/n.indpt

	# Use the Student's T distribution in place of the Std Normal if J<41
	df<- ifelse(n.indpt<41, (n.indpt - 2), NA)
	
	inference<- get.Inference(psi=psi, psi.hat=est$ATE, var= var.IC,  df=df)
	Estimator<- data.frame(est$TxtSpecificMean1,  est$TxtSpecificMean0,  est$ATE, inference, QAdj,  gAdj, work.model) 
	if(is.na(psi)){
		# then application where don't know the truth
		colnames(Estimator)<- c('Risk1', 'Risk0', 'RiskDiff',  'var', 'tstat', 'CI.lo', 'CI.hi', 'pval', 'QAdj', 'gAdj', 'WorkModel')		
	}else{
		# if simulation where we know the truth	
		colnames(Estimator)<- c('Risk1', 'Risk0', 'RiskDiff',  'var', 'tstat', 'cover', 'reject', 'QAdj', 'gAdj', 'WorkModel')		
	}

	Estimator
}

#-----------------------------------------------------#-----------------------------------------------------
# doTMLE: function to run full hierarchical TMLE under the submodel or nonparametric model 
#
# input: data (train), 
#		indicator if initial estimation of the conditional mean outcome is at the individual-level (Qinit.Indv),
#		adjustment variables for the conditional mean outcome (QAdj),
#		indicator if the working model is considert to hold (i.e. estimation under smaller sub-model),
#		adjustment variables for the propensity score (gAdj),
#		marginal probability of the treatment (prob.txt),
#		indicator to run the full SuperLearner,
#		SuperLearner libraries for the conditional mean outcome (SL.library.Q, SL.library.g),
#		indicator to print updates (verbose)
# output: 
#		clever covariate, 
#		targeted predictions, 
#		risk under treatment E[E(Y|A=1,W)], 
#		risk under control E[E(Y|A=0,W)]
#   	risk difference (point estimate of treatment effect), 
#		initial estimator of the conditional mean outcome,
#		estimator of the propensity score,
#		estimated fluctuation coefficient,
#		indicator to do the targeting
#-----------------------------------------------------#-----------------------------------------------------

doTMLE<- function(train, Qinit.Indv=F, QAdj,  work.model, gAdj=NULL, prob.txt,
	Do.SuperLearner=F, SL.library.Q, SL.library.g, 
	verbose=F){	
	
	# if no weight vector add a dummy one
	if( sum(grep('alpha', colnames(train)))==0){
		train <- cbind(train, alpha=1)
	} 		
		
	#------------------------------------------------------------------------------
	# Step1 - initial estimation of E(Y|A,W)= Qbar(A,W)
	#-----------------------------------------------------------------------------

	# if running entire TMLE algorithm at cluster-level, aggregate data now
	if(!Qinit.Indv){  			
		train<- aggregate(x=train, by=list(id=train$id), mean)[,2:(ncol(train)+1)]
		# weights for aggregated data should be 1
		train$alpha<- 1
	}
	# otherwise, initial estimation of the mean outcome Qbar is at the individual level
	
	if(Do.SuperLearner){ 
		# run full SuperLearner
		Q<- do.Init.Qbar.SL(train=train, QAdj=QAdj, SL.library=SL.library.Q, verbose=verbose)
		Q.out <- Q$outcome.SL
	} else{
		# run glm on the adjustment set
		Q<- do.Init.Qbar(train=train, QAdj=QAdj, verbose=verbose)
		Q.out=Q$glm.out
	}		

	
	# if we believe our working model (i.e. if estimating under the submodel)
	if(work.model){  
		# do not aggregate initial predictions to the cluster-level
		QbarAW.train<- Q$Qinit$QbarAW.train
		Qbar1W.train<- Q$Qinit$Qbar1W.train
		Qbar0W.train<- Q$Qinit$Qbar0W.train
		A.train<- train$A
		Y.train<- train$Y
	
	} else{ 
		# aggregate initial predictions to the cluster-level
		QC<- aggregate(x=Q$Qinit, by=list(id=train$id), mean)			
		QbarAW.train<- QC$QbarAW.train
		Qbar1W.train<- QC$Qbar1W.train
		Qbar0W.train<- QC$Qbar0W.train		
		
		train<- aggregate(x=train, by=list(id=train$id), mean)[,2:(ncol(train)+1)]
		A.train<- train$A
		Y.train<- train$Y
		# weights for aggregated data should be 1
		train$alpha <-1
	}
	

	#---------------------------------------------------------------------
	# Step2: Calculate the clever covariate
	#--------------------------------------------------------------------

	if( is.null(gAdj) ){	
		pscore.train <-  rep(prob.txt, length(A.train) )
		p.out<- NA
	} else if ( length(gAdj)==1 & gAdj[1]=='U'){	
		pscore.train <-  rep(prob.txt,  length(A.train)  )
		p.out<- NA
	} else{  
		# fit pscore on training set 	
		gAdj.train <- train[, gAdj]

		if(Do.SuperLearner){
			# if running SuperLearner
			p.out <- SuperLearner(Y=A.train, X=as.data.frame(gAdj.train),  family="binomial",
				SL.library= SL.library.g, id=train$id,
				cvControl=list(V=5), obsWeights=train$alpha, verbose=verbose)
			pscore.train <- p.out$SL.predict
			
		}else{
			# run glm on the adjustment set
			data.temp<- data.frame(gAdj.train, A=A.train)
			colnames(data.temp)<- c(gAdj, 'A')
			p.out<-   suppressWarnings( glm( A~. , family='binomial', data=data.temp, weights=train$alpha) )
			if(verbose) print(p.out)
			pscore.train <- predict(p.out, newdata=data.temp,  type="response")
		}
		
	}

	# bound g
	pscore.train[pscore.train< 0.025] <- 0.025
	pscore.train[pscore.train> 0.975] <- 0.975
	
	# Clever covariate corresponding to the GComp identifiability result E[ E(Y|A=1,W)- E(Y|A=0,W)]
	#	for the population average treatment effect E[ Y(1) - Y(0)]
	H.AW.train<- A.train/pscore.train - (1-A.train)/(1-pscore.train)

	#---------------------------------------------------------------------------
	# Step3: Targeting
	#---------------------------------------------------------------------------

	# Skip the updating step if outcome=0 in either txt or control group
	Skip.update<-  mean(Y.train[A.train==1])==0 | mean(Y.train[A.train==0])==0 |  mean(Y.train[A.train==1])==1 | mean(Y.train[A.train==0])==1 
	
	if( !Skip.update  ){
		do.update<- T
		# updating should be the same whether its at the indv or cluster
		
		# Note: updating on the logistic scale is strongly recommended. See Gruber & van der Laan 2010
		logitUpdate<- suppressWarnings( glm(Y.train ~ -1 +offset(qlogis(QbarAW.train)) + H.AW.train, family="binomial", weights=train$alpha))
		if(verbose) print(logitUpdate)
		eps<-logitUpdate$coef
	
		# updated QbarAW estimates for training set. 
		QbarAW.train<- plogis( qlogis(QbarAW.train) + eps*H.AW.train)	
		Qbar1W.train<- plogis( qlogis(Qbar1W.train) + eps/pscore.train)
		Qbar0W.train<- plogis( qlogis(Qbar0W.train) - eps/(1-pscore.train))
	} else{
		do.update<- F	
		eps<- 0
	}
	
	#----------------------------------------------------------------------------------
	# Step4: Parameter estimation
	#-----------------------------------------------------------------------------------

	# Point estimation 
	#	if believe the working model (i.e. estimating under the submodel), need to aggregate first within clusters
	# then take the empirical mean across clusters.
	TxtSpecificMean1<- mean( aggregate(data.frame(Qbar1W.train), by=list(train$id), mean)[,2] )
	TxtSpecificMean0<- mean( aggregate(data.frame(Qbar0W.train), by=list(train$id), mean)[,2] )
	ATE<- TxtSpecificMean1 - TxtSpecificMean0
	
    		
	RETURN<- list(H.AW=H.AW.train, QbarAW=QbarAW.train, 
		 	 Qbar1W=Qbar1W.train, Qbar0W=Qbar0W.train, TxtSpecificMean1=TxtSpecificMean1, 
		 	 TxtSpecificMean0=TxtSpecificMean0, ATE=ATE,
			 Q.out=Q.out, p.out=p.out, eps=eps, do.update=do.update)	
	RETURN
}




#-----------------------------------------------------#-----------------------------------------------------
# do.Init.Qbar - function to do initial estimation of the outcome regression
# 	input: training set, adjustment variable(s), verbose
# 	output: glm.fit, data.frame w/ initial predictions: Qbar(A,W), Qbar(1,W) and Qbar(0,W)
#-----------------------------------------------------#-----------------------------------------------------
do.Init.Qbar<- function(train, QAdj, verbose=F){
	
	X1<- X0<- train
	X1$A<-1; X0$A<- 0	
	train.temp<- train[, c(QAdj, 'A', 'Y')]

	# this could be indv or cluster level 
	glm.out<- suppressWarnings( glm( as.formula(Y~.), family='binomial', data=train.temp, weights=train$alpha ) )	
	if(verbose) print(glm.out)
	QbarAW.train <- predict(glm.out, newdata=train, type='response')
	Qbar1W.train<-  predict(glm.out, newdata=X1, type='response')
	Qbar0W.train<-  predict(glm.out, newdata=X0, type='response')
	
	list(glm.out=glm.out, Qinit=data.frame(QbarAW.train, Qbar1W.train, Qbar0W.train))
}	


#-----------------------------------------------------#-----------------------------------------------------
# do.Init.Qbar.SL - function to do initial estimation of the outcome regression with SuperLearner
# input: training set, adjustment variable(s), SuperLearner library, verbose
# output: glm.fit, data.frame w/ initial predictions Qbar(A,W), Qbar(1,W) and Qbar(0,W)
#
# Disclaimer: under the general (np) model, care must be taken to specify both the
#		SuperLearner loss function and library. under this statistical model, candidate estimators 
#		should target the conditional mean of the cluster-level mean Qbar^c
#		(e.g. averages of individual-level parametric or data-adaptive regressions)
# 		and the loss function should be optimized for Qbar^c (Appendix C)
#-----------------------------------------------------#-----------------------------------------------------
do.Init.Qbar.SL<- function(train, QAdj, SL.library, verbose=F){
	
	# now do initial estimation at the cluster level
	train.temp <- train[, c(QAdj, 'A')] 
	X1<- X0<- train.temp
	X1$A<-1; X0$A<- 0
	
	outcome.SL<- SuperLearner(Y=train$Y, X=train.temp, family='binomial', SL.library= SL.library, 
			id=train$id,	cvControl=list(V=5), obsWeights=train$alpha, verbose=verbose)
	QbarAW.train<- predict(outcome.SL, newdata=train.temp)$pred
	Qbar1W.train<- predict(outcome.SL, newdata=X1)$pred
	Qbar0W.train<- predict(outcome.SL, newdata=X0)$pred	
	list(outcome.SL=outcome.SL, Qinit=data.frame(QbarAW.train, Qbar1W.train, Qbar0W.train))
}	


#-----------------------------------------------------#-----------------------------------------------------
# get.Inference: function to calculate 95% confidence intervals and test the null hypothesis ATE=0
#	input: psi (true value of PATE), psi.hat (point estimate), var (variance estimate), 
#		df (degrees of freedom if using a Student's t-dist ) 
# output: 
#		if truth known (i.e. simulation): variance, test statistic, 95% CI coverage, indicator of rejecting the null
#		if truth unknown (i.e. application) :variance, test statistic, confidence intervals, pval,
#-----------------------------------------------------#-----------------------------------------------------	
get.Inference<- function(psi=NA, psi.hat, var, df=NA){
  
  	# assumed signifcance level of 0.05
  	sig.level=0.05
	if(is.na(df)){
   	 	cutoff<- qnorm(sig.level/2, lower.tail=F)
	}else{ 
  		# cutoff based on t-dist for testing and CI	
  		cutoff <- qt(sig.level/2, df=df, lower.tail=F)
  	}
  	
 	# standard error (square root of the variance)
  	se<- sqrt(var)
  	
  	# test statistic and pvalue
  	tstat <- psi.hat/se
  	
  	# 95% confidence interval
  	CI.lo <- psi.hat - cutoff*se
  	CI.hi <-  psi.hat + cutoff*se
  	# p-value 
  	if(is.na(df)){
  		pval <- 2*pnorm(abs(tstat), lower.tail=F) 
  	}else{
  		pval<- 2*pt(abs(tstat), df=df, lower.tail=F) 
    }
  	
  	# if application vs. simulation
  	if(is.na(psi)){
  		# application where the truth is not known
 		out<-   	data.frame(var, tstat, CI.lo, CI.hi, pval)

  	} else{
  		# simulation where truth is known 
  	
  		# confidence interval coverage
  		cover<- ( CI.lo <= psi & psi <= CI.hi )
  		
  		# reject the null
  		reject <- pval < sig.level 
  
  		out<- data.frame(var, tstat, cover, reject)
  	}
  	out 
}



#-----------------------------------------------------#-----------------------------------------------------
# do.adaptive.prespec: function to implement adaptive prespecification as described in 
#		Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching"
#	input: dataset, indv-level and cluster-level adjustment variables
#		indicator working model is considered to hold (i.e. estimation under the sub-model)
#		marginal probability of the exposure
#	output: selection for candidate TMLE
#-----------------------------------------------------#-----------------------------------------------------
do.adaptive.prespec<- function(data, ind.adj, clust.adj, work.model, prob.txt ){

	#do adaptive prespecification with indv and cluster covariates
	select.Q<- suppressWarnings( CV.selector(data=data, ind.adj=ind.adj, clust.adj=clust.adj, 
						work.model=work.model, forQ=T, prob.txt=prob.txt) )
	# print(select.Q)			
	if(select.Q$Adj!='U'){ 
		# if did not select the unadjusted for initial estimation of Qbar(A,W), 
		# then need to remove this variable from the candidate for pscore estimation
		if(select.Q$Ind){
			ind.adj<- ind.adj[-which(ind.adj==select.Q$Adj) ]
		} else{
			clust.adj<- clust.adj[-which(clust.adj==select.Q$Adj)]
		}
	}
	select.G<- suppressWarnings( CV.selector(data=data,  ind.adj=ind.adj, clust.adj=clust.adj, 
			work.model=work.model, forQ=F, 
			Qinit.Indv= select.Q$Ind, QAdj=select.Q$Adj, 
			prob.txt=prob.txt))
	# print(select.G)
	list(Qinit.Indv=select.Q$Ind, QAdj=select.Q$Adj, gAdj=select.G$Adj)
}




#-----------------------------------------------------#-----------------------------------------------------
# CV.selector: function to estimate the cross-validated risk
#		Loss function is the squared-IC; Risk is then the variance of the TMLE
#		See Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching"
#
#	input: dataset, indv-level and cluster-level adjustment variables
#		indicator if for the conditional mean outcome (forQ)
#		indicator if outcome regression is fit at indv-level (Qinit.Indv)
#		selected adjustment variable for the outcome regression (QAdj)
#		indicator working model is considered to hold (i.e. estimation under the sub-model)
#		marginal probability of the exposure
#	output: selection for adjustment variable (corresponding to a TMLE)
#-----------------------------------------------------#-----------------------------------------------------
CV.selector <- function(data, ind.adj, clust.adj,forQ, Qinit.Indv, QAdj, work.model,  prob.txt ){

	# covariates correspond to different TMLEs  
	adj.set<- c(ind.adj, clust.adj)
	Qinit.Indv.temp<- c(rep(T, length(ind.adj)), rep(F, length(clust.adj))   )

	num.tmles <- length(adj.set)
	CV.risk <- rep(NA, num.tmles)
	
	for(k in 1: num.tmles){	

		if(forQ){
			# if selecting the adjustment set for the outcome regression
			IC.temp<- getIC.CV(data=data,  Qinit.Indv=Qinit.Indv.temp[k], QAdj=adj.set[k], work.model=work.model, gAdj=NULL, prob.txt=prob.txt)
		} else{
			# if collaboratively selecting the adjustment set for the pscore
			IC.temp<- getIC.CV(data=data, Qinit.Indv=Qinit.Indv, QAdj=QAdj,  work.model=work.model, gAdj=adj.set[k], prob.txt=prob.txt)
		}
		# estimating the CV risk for each candidate
		CV.risk[k]<- mean(IC.temp$Risk)
	}
	# select the candidate estimator resulting in the smallest CV-risk
	this<- which.min(CV.risk)
 	list(Ind= Qinit.Indv.temp[this], Adj=adj.set[this ] )
}


#-----------------------------------------------------#-----------------------------------------------------
# getIC.CV: function to obtain a cross-validated estimate of the influence curve
#		See Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching"
#
#	input: dataset, 
#		indicator if outcome regression is fit at indv-level (Qinit.Indv),
#		selected adjustment variable for the outcome regression (QAdj),
#		indicator working model is considered to hold (i.e. estimation under the sub-model),
#		selected adjustment variable for the pscore (gAdj),
#		marginal probability of the exposure
#	output: cross-validated estimate of the IC,
#		cross-validated risk estimate (loss=IC^2)
#-----------------------------------------------------#-----------------------------------------------------
getIC.CV<- function(data, Qinit.Indv, QAdj,  work.model, gAdj, prob.txt){
	
	#******** Disclaimer: this is implemented for leave-one-CLUSTER-out cross validation  *********
	#	but could easily be modified for other CV schemes
	#	Need to respect the experimental unit as the unit of idpt
	Folds<- data$id
	units <- unique(data$id)
	nFolds<- length(units)   
	
	IC<- Risk <- rep(NA,  nFolds)

	# doing leave-one-out cross-validation
	for(i in 1: nFolds) {
		
		valid <- data[Folds==units[i] , ]
		train <- data[Folds!=units[i], ]
	
		# run full TMLE algorithm on the training set
		train.out<- doTMLE(train=train, Qinit.Indv=Qinit.Indv, 
			QAdj=QAdj, work.model=work.model,  gAdj=gAdj, prob.txt=prob.txt, Do.SuperLearner=F, verbose=F)		
	
		# get the relevant components of the IC for the validiation set, using fits based on the training set
		valid.out<- do.TMLE.for.valid(valid=valid, Qinit.Indv=Qinit.Indv,
			QAdj=QAdj, work.model=work.model, gAdj=gAdj, prob.txt=prob.txt, 
			Q.out= train.out$Q.out, p.out= train.out$p.out, eps= train.out$eps, do.update= train.out$do.update, ATE= train.out$ATE)
		
		# for cross-validation selector, the loss fuction is the IC-squared	
		IC.sq <- valid.out$IC*valid.out$IC
		# empirical cross-validated risk estimate for fold i	
		# (not really needed bc leave-one-out )
		Risk[i] <- mean(IC.sq)
		
		# cross-validated IC estimate (already aggregated to cluster level)
		# again assumes LOOCV
		IC[i] <-  valid.out$IC
	
	}
	
	list(IC=IC, Risk=Risk)
}



#-----------------------------------------------------#-----------------------------------------------------
# do.TMLE.for.valid: function to obtain a cross-validated estimate of the influence curve
#	for observations in the validation set
#		See Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching"
#
#	input: dataset, 
#		indicator if outcome regression is fit at indv-level (Qinit.Indv),
#		selected adjustment variable for the outcome regression (QAdj),
#		indicator working model is considered to hold (i.e. estimation under the sub-model),
#		selected adjustment variable for the pscore (gAdj),
#		marginal probability of the exposure,
#		initial estimator of the outcome regression - fit on the training set (Q.out)
#		estimator of the pscore - fit on the training set (p.out)
#		estimated fluctuation coefficient - fit on the training set (eps)
#		indicator to the update (do.update)
#		point estimate of the ATE
#	output: cross-validated estimate of the IC,
#		cross-validated risk estimate (loss=IC^2)
#-----------------------------------------------------#-----------------------------------------------------

do.TMLE.for.valid<- function(valid, Qinit.Indv, QAdj, work.model, gAdj, prob.txt, Q.out, p.out, eps, do.update=T, ATE) {
	
	#------------------------------------------------------------------------------
	# Step1 - initial estimation of E(Y|A,W)= Qbar(A,W)
	#-----------------------------------------------------------------------------

	# if running entire TMLE algorithm at cluster-level 
	if(!Qinit.Indv){  	
		valid<- aggregate(x=valid, by=list(id=valid$id), mean)
	}

	V1<- V0<- valid
	V1$A= 1; V0$A=0
	QbarAW.valid <-  predict(Q.out, newdata=valid, type='response')
	Qbar1W.valid <-  predict(Q.out, newdata=V1, type='response')
	Qbar0W.valid <-  predict(Q.out, newdata=V0, type='response')

	if(!work.model){  
		# if estimating under the larger np model (i.e. treating the assumptions as working)
		QbarAW.valid <-  aggregate(QbarAW.valid, by=list(id=valid$id), mean)[,2]
		Qbar1W.valid <-  aggregate(Qbar1W.valid, by=list(id=valid$id), mean)[,2]
		Qbar0W.valid <-  aggregate(Qbar0W.valid, by=list(id=valid$id), mean)[,2]
	
		valid<- aggregate(x=valid, by=list(id=valid$id), mean)[,2:ncol(valid)]

	}
	A.valid <- valid$A
	Y.valid <- valid$Y
	
	#--------------------------
	# Step2: Calculate the clever covariate
	#---------------------------
	if( is.null(gAdj)){
		pscore.valid<- rep( prob.txt, length(A.valid))
	} else if (length(gAdj)==1 & gAdj[1]=='U'){	
		pscore.valid<- rep( prob.txt, length(A.valid))
	} else{
		data.temp<- valid[,c(gAdj, 'A')]
		colnames(data.temp)<- c(gAdj, 'A')
		pscore.valid<- predict(p.out, newdata=data.temp,  type='response')
	}
	# bound g
	pscore.valid[pscore.valid< 0.025] <- 0.025	
	pscore.valid[pscore.valid> 0.975] <- 0.975	

	# clever covariate for observations in the validation set
	H.AW.valid<- A.valid/pscore.valid - (1-A.valid)/(1-pscore.valid)

	#-----------------------------
	# Step3: Targeting
	#------------------------------
		
	if(do.update){
				
		# update
		QbarAW.valid<- plogis( qlogis(QbarAW.valid) + eps*H.AW.valid)	
		Qbar1W.valid<- plogis( qlogis(Qbar1W.valid) + eps/pscore.valid)
		Qbar0W.valid<- plogis( qlogis(Qbar0W.valid) - eps/(1-pscore.valid))	
		
	}	
	# obtain a CV-estimate of the IC
	IC <- H.AW.valid*(valid$Y - QbarAW.valid) + Qbar1W.valid - Qbar0W.valid - ATE
	# aggregate to the cluster-level 
	IC <- aggregate(data.frame(IC), by=list(id=valid$id), mean)
	IC
}


