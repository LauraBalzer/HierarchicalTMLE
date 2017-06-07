#-----------------------------#-----------------------------
# R code to implement the hierarchical TMLEs in the Network-Based Simulations
#		described in Balzer et al: "A New Approach to Hierarchical Data Analysis: 
#		Targeted Maximum Likelihood Estimation of Cluster-Based Effects Under Interference
#
# Simulation 2 - Cluster Randomized Trial for HIV Prevention & Treatment
# 		Python code by Dr. Patrick Staples to simulate the data is available 
#		@ https://github.com/ctphoenix/HIV-PrEP-Simulation
#
# Programer: Laura Balzer (lbbalzer@gmail.com)
#
# Last update: June 2, 2017
#-----------------------------#-----------------------------

#-------------------------------------------
rm(list=ls())

set.seed(1)

library('SuperLearner')
source('Simulation2_Functions_v2June2017.R')
source('ClusterTMLE_Functions_v2June2017.R')

# Indicator if there is an effect (vs. the null)
effect <- F

# Number of clusters in each study
J <- 32
# Number of individuals per cluster (constant)
n<- 75

# Number of repetitions of the data generating experiment
nReps <- 1000

#--------------------
# Read in true values of the risk under txt, control and population average treatment effect
if(effect){
	load('Simulation2_TRUTH_J32_n75_nReps1000_v04Apr2017.Rdata') 
	psi <- PATE
} else{
	psi <- 0 
} 
#----------------------
#file.start<- '~/Dropbox/Data For Laura/20160801/Design 1/'
file.start <- 'Data20160801/'

file2<- paste('Simulation2', '_effect', effect,  '_J',J, '_n', n,  '_nReps', nReps, '_v', format(Sys.time(), "%d%b%Y"), '.Rdata', sep='')

save( psi, file=file2)


unadj <- Gcomp <- IPTW<- TMLE.WorkF <- TMLE.WorkT <- NULL

for(i in 1:nReps ){		

	trial.number<- i-1  
	file<- paste(file.start,"CRT_", trial.number, ".txt", sep="")
	
	output<- do.network.analysis(file=file, psi=psi, effect=effect, J=J,  n=n, verbose=F)
	unadj <- rbind(unadj, output['unadj', ] )	
	Gcomp <- rbind(Gcomp, output['Gcomp', ] )	
	IPTW <- rbind(IPTW, output['IPTW', ] )	

	TMLE.WorkF <- rbind(TMLE.WorkF, output['TMLE.WorkF', ] )	
	TMLE.WorkT <- rbind(TMLE.WorkT, output['TMLE.WorkT', ] )	

	save(psi, unadj, Gcomp, IPTW,  TMLE.WorkF, TMLE.WorkT,  file=file2)

	print(file)

}

save(psi,  unadj, Gcomp, IPTW, TMLE.WorkF, TMLE.WorkT,  file=file2)

