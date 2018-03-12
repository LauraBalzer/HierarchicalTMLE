#-----------------------------#-----------------------------#-----------------------------#-----------------------------
# R code to generate simulated data and implement the hierarchical TMLEs
#		described in Balzer et al: "A New Approach to Hierarchical Data Analysis: 
#		Targeted Maximum Likelihood Estimation of Cluster-Based Effects Under Interference"
#
# Simulation 1 - Simple Observational Settings
#
# Programer: Laura Balzer (lbalzer@umass.edu)
#
# Last update: Feb 23, 2018
# Explore scenarios +/- covariate interference +/- dependent UY
#-----------------------------#-----------------------------#-----------------------------#-----------------------------


rm(list=ls())
library('MASS')
source('ClusterTMLE_Functions_v2June2017.R')
set.seed(1)

#-----------------------------#-----------------------------#-----------------------------#-----------------------------
# getY: function to generate the outcome
#	note: the coefficients beta are global variables
#-----------------------------#-----------------------------
getY<- function( A, W, Z, UY){
	as.numeric(UY < plogis(.25 + betaA*A + betaWc*mean(W) + betaZc*mean(Z) + betaW*W + betaZ*Z) )
}


#-----------------------------#-----------------------------#-----------------------------#-----------------------------
# generate.cluster.data:
#		input:
#			indicator of dependence in the unmeasured factors
#			indicator of effect (or under null), 
#			number of indv per cluster (n)
#			treatment indicators (A), 
#			number of clusters (j)
#		output: data.frame of observed data and the counterfactual outcomes for one cluster
#-----------------------------#-----------------------------
generate.cluster.data<- function(deptU, effect, n, A, j){

	UE <- runif(1, -1, 1)
	W<- rnorm(n, UE, .5) 
	Z <- rnorm(n, UE, .5)
	if(is.na(A)){
		A <- rbinom(1, 1, prob=plogis( 0.75*mean(W))	)
	}
	
	if(deptU){
		UY <- pnorm( rnorm(n, UE, .3)) 
	}else{
		UY <- runif(n, 0,1)
	}
	Y0 <- getY(0, W, Z, UY)
	
	if(effect){
		Y1 <- getY(1, W, Z, UY)
	} else{ 
		#if under the null, then set the counterfactual outcome Y1 to Y0
		Y1 <- Y0
	}
	
	if(A){
		# if exposed cluster
		Y <- Y1
	}else{
		# if unexposed cluster
		Y <- Y0	
	}
	# calculate the weights as 1/n_j
	alpha <- 1/n
	# return data.frame with id as the cluster indicator and dummy variable U=1 (for unadjusted)
	data.frame( cbind(id=j, n=n, U=1,  W=W, Wc=mean(W), Z=Z, Zc=mean(Z), A=A, Y, alpha, Y1, Y0) )
}

#-----------------------------#-----------------------------#-----------------------------#-----------------------------
# get.full.data:
#		input: number of clusters (j),  number of indv per cluster (n), 
#			indicator of randomized trial (trial), indicator of dependence in unmeasured factors (deptU)
#			 indicator of effect vs. the null (effect)
#		output: data.frame of observed data and the counterfactual outcomes 
#-----------------------------#-----------------------------
get.full.data <- function(J, n, deptU, trial, effect){
	
	if(trial){
		# if randomized trial, then assign treatment with equal allocation (J/2 treated & J/2 control)
		A.1 <- rbinom(J/2, 1, 0.5)
		A.2 <- ifelse(A.1==1, 0, 1)
		A <- sample( c(A.1,A.2))
	}else{
		A<- rep(NA, J)
	}
	
	n.indv<- round(rnorm(J, n, 10))
	# deal with negative values: set to n if less than 1. 
	n.indv[n.indv<1] <- n

	full.data<- NULL
	for(j in 1: J){
		data.comm.j <- generate.cluster.data(deptU=deptU, effect=effect, n=n.indv[j], A=A[j], j) 
		full.data <- rbind(full.data, data.comm.j)
	}
	full.data

}

#--------------------------------------------------#--------------------------------------------------#-----------------------------#-----------------------------
#--------------------------------------------------#--------------------------------------------------#-----------------------------#-----------------------------

# SPECIFY THE PARAMETERS FOR THE SIMULATION

# Indicator if randomized trial setting
trial <- F

# Indicator if there is an effect (vs. the null)
effect<- F

# Indicator that there is covariate interference
interference <-  T

# Indicator that dependent Us
deptU <- T



# Number of clusters in the target population
nPop<- 10000
# Number of clusters in each study
J <- 100
# Average number of individuals per cluster
n <- 50
# Number of repetitions of the data generating experiment
nReps<- 5000

file.name<- paste('Sim1_', 'observational', !trial, 'interference', interference,  'deptU', deptU, '_effect', effect, '_J',J, '_n', n, '_nReps', nReps,'_v', format(Sys.time(), "%d%b%Y"), '.Rdata', sep='')
#--------------------------------------------------#--------------------------------------------------#-----------------------------#-----------------------------
#--------------------------------------------------#--------------------------------------------------#-----------------------------#-----------------------------



if(!interference){
	betaA <<- 0.1 #change (Feb23) to decrease power
	betaWc <<-  0.15 	
	betaW  <<- 1.15 
	betaZc<<-  0
	betaZ  <<-  1
} else{
	betaA <<- 0.1 #change (Feb23) from 0.2
	betaWc <<-  0.15  	
	betaW <<- 0.25
	betaZc <<-  1
	betaZ <<- 0
}


#--------------------------------------------------#--------------------------------------------------#-----------------------------#-----------------------------
# Get the true value of the average treatment effect
if(effect){
	Pop<- get.full.data(J=nPop, n=n, deptU=deptU, trial=trial, effect=effect)
	PopC<- aggregate(Pop, by=list(Pop$id), mean) 
	truth <- mean(PopC$Y1 - PopC$Y0)	
} else{
	truth<- 0
	PopC<- nPop<- NA
}
save(nPop, PopC, truth, file=file.name)



#--------------------------------------------------#--------------------------------------------------#-----------------------------#-----------------------------
#--------------------------------------------------#--------------------------------------------------#-----------------------------#-----------------------------


# Repeat the data generating process nReps times and implement the relevant estimators

unadj<- tmle1a <- tmle1b <-  tmle2 <- NULL
verbose=F

for(r in 1:nReps){
	
	# Generate the full hierarchical data
	data<- get.full.data(J, n, deptU=deptU, trial=trial, effect=effect)

	prob.txt <- mean( aggregate(data, by=list(data$id), mean)$A)
	
	# unadjusted estimatior
	clustU <- do.estimation.inference(psi=truth, data=data, Qinit.Indv=F, QAdj='U',  work.model=F, gAdj=NULL, prob.txt=prob.txt, verbose=verbose)
	unadj<- rbind(unadj, clustU)

	# TMLEs with varying working assumptions
	QAdj <-  c('W', 'Z')
	gAdj<- c('W', 'Z')
	
	# TMLE-Ia - clust level
	clust1<- do.estimation.inference(psi=truth, data=data, Qinit.Indv=F, QAdj=QAdj, work.model=F, gAdj=gAdj, verbose=verbose)
	tmle1a <- rbind(tmle1a, clust1)
	
	# note this above is equivalent to using an individual-level regression that only adjusts for cluster=level summaries when estimating QbarC
	#		See Appendix C
	#  	do.estimation.inference(psi=truth, data=data, Qinit.Indv=T, QAdj=c('Wc', 'Zc'), work.model=F, gAdj=c('Wc','Zc'), verbose=verbose)
	
	# TMLE-Ib  - still under larger causal model but with working assumptions
	ind1<- do.estimation.inference(psi=truth, data=data, Qinit.Indv=T, QAdj=QAdj, work.model=F, gAdj=gAdj, verbose=verbose)
	tmle1b <- rbind(tmle1b, ind1)
	
	# Note an alternative not explored here is including both indiv & clust when fitting QbarAW 
	# i.e. setting 	QAdj <- c('W', 'Z', 'Wc', 'Zc')

	# under the sub-model
	# TMLE-II
	QAdj <-  c('W', 'Z')
	gAdj<- c('W', 'Z')
	ind2<- do.estimation.inference(psi=truth, data=data, Qinit.Indv=T, QAdj=QAdj, work.model=T, gAdj=gAdj, verbose=verbose)
	tmle2 <- rbind(tmle2, ind2)
	
				
	save(nPop, PopC, truth, unadj, tmle1a, tmle1b, tmle2, file=file.name)
	
	print(r)

}

save(nPop, PopC, truth, unadj, tmle1a, tmle1b, tmle2, file=file.name)

colMeans(unadj[, 1:7])
colMeans(tmle1a[, 1:7])
colMeans(tmle1b[, 1:7])
colMeans(tmle2[, 1:7])




