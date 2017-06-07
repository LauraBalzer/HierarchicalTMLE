#-----------------------------#-----------------------------
# R code to implement the hierarchical TMLEs in the Network-Based Simulations
#		described in Balzer et al: "A New Approach to Hierarchical Data Analysis: 
#		Targeted Maximum Likelihood Estimation of Cluster-Based Effects Under Interference"
#
# In this trial setting, implements adaptive prespecification as described in 
#		Balzer et al. "Adaptive pre-specification in randomized trials with and wihtout pair-matching"
#
# Programer: Laura Balzer
# Contact: lbbalzer@gmail.com
#
# Last update: June 2, 2017
#-----------------------------#-----------------------------

#-----------------------------#-----------------------------
# do.network.analysis: master function read in simulated trial data
#		and do an analysis using adaptive pre-specification
#	input: file, true value of the causal parameter, 
#		number of clusters, number of individuals per cluster, 
#		indicator of whether to print the output (verbose)
# output: data.frame of point estimates and inference
#-----------------------------#-----------------------------
do.network.analysis<- function(file, psi, effect, J, n, verbose=T){
	
	
	# assign Test-and-Treat intervention (A=1) vs. control (A=0)
	A<- rbinom(J/2, 1, 0.5)
	A.2<- ifelse(A==1, 0, 1)
	A <- sample( c(A,A.2))
	
	data<- read.full.data.ind(file=file,  txt='1_3', con='1_1', effect=effect, A=A, n=n)
	
	# TMLE code requires the treatment probability
	prob.txt <- 0.5

	# Unadjusted - controls for dummy variable 'U' 
	unadj<- do.estimation.inference(psi=psi, data=data,  
		Qinit.Indv=F, QAdj='U', 	work.model=F,
		gAdj=NULL, prob.txt=prob.txt, verbose=verbose)
		
	# cluster-level Gcomp adjusting in outcome regression for average degree
	# in a trial, this is a TMLE; See Rosenblum & van der Laan 2010	
	Gcomp <- do.estimation.inference(psi=psi, data=data,  
		Qinit.Indv=F, QAdj='degree', 	work.model=F,
		gAdj=NULL, prob.txt=prob.txt, verbose=verbose)
	
	# cluster-level Augmented-IPTW adjusting in pscore regression for average degree
	IPTW <- do.estimation.inference(psi=psi, data=data,  
		Qinit.Indv=F, QAdj='U', 	work.model=F,
		gAdj='degree', prob.txt=prob.txt, verbose=verbose)

	# specify set of candidate adjustment variables
	ind.adj  <- c('U','degree', 'sex.work', 'neigh.prev')
	clust.adj <- c('prev.base', 'assort', 'num.component')

	# TMLE under the larger non-parametric model (i.e. treating the assumptions as working)
	TMLE.WorkF <- do.estimation.inference(psi=psi, data=data, 
		ind.adj=ind.adj, clust.adj=clust.adj, work.model=F,
		prob.txt=prob.txt, Do.AdaptivePrespec=T, verbose=verbose)
	
	# TMLE under the smaller sub-model (i.e. treating the assumptions as true)
	TMLE.WorkT<- do.estimation.inference(psi=psi, data=data, 
		ind.adj=ind.adj, clust.adj=clust.adj, work.model=T,
		prob.txt=prob.txt, Do.AdaptivePrespec=T, verbose=verbose)
			
	RETURN<- data.frame(rbind(
			unadj=unadj, Gcomp=Gcomp, IPTW=IPTW, 
			TMLE.WorkF=TMLE.WorkF,  TMLE.WorkT=TMLE.WorkT))
	
	RETURN
}


#-----------------------------#-----------------------------
# read.full.data.ind  - function to read and process in the network-based simulation data
#	input: file name, treatment and control conditions, indicator of effect (or not)
#		observed exposure, number of individual per cluster
#	output: stacked data frame with id=indpt unit
#
# For further details on the underlying network-based + epidemic simulation, see
#		https://github.com/ctphoenix/HIV-PrEP-Simulation
#-----------------------------#-----------------------------

read.full.data.ind<- function(file, txt, con,  effect, A,  n){
	
	data<- read.table(file, header=T)

	OUT<- NULL
	
	for(v in 0:31){
		
		data.village<- data[data$Cluster==v, ]	
		
		# incidence cohort - everyone not 'prevalent' at yr0
		HIVneg <- data.village$Prevalences==0
		
		# restrict data set to incidence cohotrt
		data.red<- data.village[HIVneg, ]
		
		# randomly sample n who are HIVneg
		these<- sample(1:nrow(data.red), size=n)
		data.red<- data.red[these, ]
		# create weight as 1/n
		alpha <- rep(1/n, n)
	
		# indv-level variables
		degree<- data.red[, 'Degree']
		component <- data.red[, 'Node_Component_Size']
		mins <- data.red[, 'Mins']
		sums <- data.red[, 'Sums']
		sex.work <- as.numeric(data.red[, 'Sex_Work']=='True')
	
		# variables that are a function of friends		
		# average number of partners of your friends
		neigh.degree <- data.red[, 'Mean_Neighbor_Degree']
	    # number of neighbors infected at BL
		neigh.prev <- data.red[, 'Total_Neighbor_Seeds']

		# community-level variable - these are constant within a communitty
		assort<- data.red[,  'Assortativity']
		LCC <- data.red[,  'LCC_Size']
		mean.component <- data.red[, 'Mean_Component_Size']
		num.component <- data.red[, 'Number_Of_Components']
		# baseline HIV prevalence
		prev.base <- mean(!HIVneg) 

		
		# when infected - 
		# baseline: prev=1
		# yr0-yr1: outcome==(7 - prev)
		# yr1-y2: outcome=6
		# yr2-yr3: outcome=5
		# yr3-yr4: outcome=4
		# yr4-yr5: outcome=3
		# yr5-yr6: outcome=2
		# yr6-yr7: outcome=1
		# never 	outcome=0

		# indicator of infected for 5,6,7 yrs among those neg at base-yr0
		# counterfactual individual-level and cluster-level outcome under the control 
		Y0 <- as.numeric( data.red[, paste('outcome', con, sep='_')] >4) 
		Yc0<- mean(Y0)
		if(effect){
			# counterfactual individual-level and cluster-level outcome under the intervention 
			Y1<- as.numeric( data.red[, paste('outcome', txt, sep='_')] >4) 
			Yc1<- mean(Y1)
		} else{
			Y1<- Y0
			Yc1 <- Yc0
		}
		
		obs.exp <- A[v+1]
	
		# now get observed outcome.
		if(obs.exp){
			Y <- Y1
			Yc <- Yc1
		} else{
			Y <- Y0
			Yc <- Yc0
		}
		
		# include an indicator of the cluster, and a dummy variable U =1
		new.matrix <- cbind(id=(v+1), U=1, 
			degree, component, mins, sums,  sex.work, 
			neigh.degree, neigh.prev, 
			assort,  LCC, mean.component, num.component, prev.base, 
			Y1, Yc1, Y0, Yc0, A=obs.exp, Y, Yc, alpha)
		OUT<- rbind(OUT, new.matrix)
	}
	OUT<- as.data.frame(OUT)
	OUT
}



