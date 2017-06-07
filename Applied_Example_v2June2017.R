#=========================================
# R Code to estimate the association of household economic status (binary)
#	on testing for HIV by the hybrid scheme in SEARCH (NCT01864603)
#
# Programmer: Laura Balzer
# lbbalzer@gmail.com
##
# Last update: June 2, 2017
#
#==========================================
rm(list=ls())

# Requires the SuperLearner v2.0-21
library('SuperLearner')

# load in necessary functions
source('ClusterTMLE_Functions_v2June2017.R')

# set the seed so that each time run all code; get same answer
set.seed(1)

load('Applied_Example_v13Apr2017.Rdata')

data<- data.new

# Primary analysis is running with SuperLearner with a large library
# * If false, then runs with small & simple library *** use for debugging ONLY
Do.SL.Big <-  T


#-------------------------------------#---------------------------------------------------

# Specify SuperLearner library
if(Do.SL.Big){ 
	 SL.library<-  c('SL.glm', 	'SL.glm.interaction', 'SL.gam',	'SL.glmnet','SL.mean')
}else{
	# If not, specified run with a small library for QC-ing
	SL.library<- c('SL.glm', 'SL.mean')
}
print(SL.library)

	
file<- paste( 'BaselineTest', paste('SuperLearnerBig', Do.SL.Big, sep=''), paste('v', format(Sys.time(), "%d%b%Y"),sep=''), 'Rdata', sep='.')
save(SL.library, file=file)


ind.adj<-  c("age.25to34", "age.35to44", "age.45plus", "male", 
	"primary", "secondary.plus", "formal.hi.occup", "informal.hi.occup", "informal.low.occup", "jobless", "single", "moAway1.plus")   
clust.adj<- c("Bugamba", "Kameke", "Kamuge", "Kazo", "Magunga", "Merikit", "Mitooma", "Muyembe", "Nankoma", 
	       "Nyatoto", "Ogongo", "Othoro", "Rubaare", "Ruhoko",  "Sena", "hh.size", "hh.male")         
QAdj<- c(ind.adj, clust.adj)
gAdj<- c(ind.adj, clust.adj)


# unadjusted; 
# need to add a dummy variable
data.temp <- cbind(data, U=1)
# also need to calculate the P(Ac=1)
prob.txt<- mean( aggregate(data.frame(A=data$A), by=list(data$id), mean)[,'A'])
unadj <- do.estimation.inference(psi=NA, data=data.temp,  
	Qinit.Indv=F, QAdj='U', work.model=F, 
	gAdj='U', prob.txt=prob.txt, Do.SuperLearner=F,  verbose=T)

# TMLE under the sub-model	
TMLE <- do.estimation.inference(psi=NA, data=data,  
	Qinit.Indv=T, QAdj=QAdj, work.model=T, 
	gAdj=gAdj, 	Do.SuperLearner=T,  SL.library.Q=SL.library, SL.library.g=SL.library,  
	verbose=T)
	
save(unadj, TMLE, file=file)