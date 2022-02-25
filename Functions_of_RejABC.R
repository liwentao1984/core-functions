source('Divide_Recombine_Parallelise.R')
source('Auxiliary_Functions.R')

##########################
# Functions	
##########################
rejABC_OneDataset<-function(observed_summary,reference_tables,tolerance_percentile=NULL,tolerance_value=NULL){
### What is it for? ###
# Accept/reject the simulated parameter values according to the weighted Euclidean distances with a single observed summary statistic, and output the accepted simulations.

### Specification of Implementation ###
# 1. Weights in the Euclidean distance are reciprocals of sample MADs of simulated summaries.

### In and Out ###
# Input: reference_tables - a list containing two elements, the N*p matrix 'parameters' containing the p-dimension parameter values and the N*d matrix 'summaries' containing the d-dimension summary statistics. observed_summary - a d*1 vector containing the summary statistic of the observed data set. tolerance_percentile - a scalar giving the proportion to be accepted in the simulated parameter values. tolerance_value - a scalar giving the value of tolerance level, and only one bewtween this and 'tolerance_percentile' should be given
# Output: a five-element list containing a N0*p matrix 'parameters', a N0*d vector 'summaries', a N0*1 vector 'acceptance_ind' with TRUE and FALSE elements, where N0 is the number of acceptance, a scalar 'tolerance' and a scalar 'acceptance_rate'

	weights_of_summaries_for_distance<-1/apply(reference_tables$summaries,2,mad) # a d*1 vector
	if(length(weights_of_summaries_for_distance)>1) weights_matrix<-diag(weights_of_summaries_for_distance)
	if(length(weights_of_summaries_for_distance)==1) weights_matrix<-as.matrix(weights_of_summaries_for_distance)
	weighted_Euclid_distances<-cal_distance_between_matrixcol_vector(t(reference_tables$summaries),observed_summary,weights_matrix^2) # a N*1 vector
	weighted_Euclid_distances[is.na(weighted_Euclid_distances)]<-max(weighted_Euclid_distances[!is.na(weighted_Euclid_distances)])
	if(!is.null(tolerance_percentile)) tolerance_value<-quantile(weighted_Euclid_distances,tolerance_percentile) # a scalar
	acceptance_rate<-mean(weighted_Euclid_distances<=tolerance_value)
	acceptance_ind<-(weighted_Euclid_distances<=tolerance_value)
	return(list(parameters=as.matrix(reference_tables$parameters[acceptance_ind,]),summaries=reference_tables$summaries[acceptance_ind,],acceptance_ind=acceptance_ind,tolerance=tolerance_value,acceptance_rate=acceptance_rate))
}

rejABC_MultiDatasets<-function(running_set_indicators,reference_tables,observed_summaries,tolerance_percentile=NULL,tolerance_value=NULL){
### What is it for? ###
# Modify the input of rejectionABC to adapt the parallelisation of running multiple sessions of rejection ABC for multiple data sets. 

### Specification of Implementation ###

### In and Out ###
# Input: running_set_indicators - a l*1 vector, giving the indicators for which observed data sets the rejection ABC is run for. reference_tables - a list containing two elements, the N*p matrix 'parameters' containing the p-dimension parameter values and the N*d matrix 'summaries' containing the d-dimension summary statistics. observed_summaries - a L*d matrix with rows containing the summary statistics of all observed data sets. tolerance_percentile - a scalar giving the proportion to be accepted in the simulated parameter values. tolerance_value - a scalar giving the value of tolerance level, and only one bewtween this and 'tolerance_percentile' should be given.
# Output: an l-element list containing the output lists of rejABC_OneDataset

	l<-length(running_set_indicators)
	if(l>1) observed_summaries_in_use<-observed_summaries[running_set_indicators,]
	if(l==1) observed_summaries_in_use<-t(observed_summaries[running_set_indicators,])
	accepted_results_for_all_set<-list(0)
	for(set_i in 1:l){
		accepted_results_for_all_set[[set_i]]<-rejABC_OneDataset(observed_summaries_in_use[set_i,],reference_tables,tolerance_percentile=tolerance_percentile,tolerance_value=tolerance_value)
		names(accepted_results_for_all_set)[set_i]<-paste0('set_',running_set_indicators[set_i])
	}
	return(accepted_results_for_all_set)
}

rejABC<-function(reference_tables,observed_summaries,tolerance_percentile=NULL,tolerance_value=NULL,paral=FALSE,ncores=if(paral) max(detectCores()-1,1) else 1, platform='Mac',packages=NULL){
### What is it for? ###
# For each of L observed datasets, accept/reject the simulated parameter values according to the corresponding weighted Euclidean distances, and output the accepted simulations.

### Specification of Implementation ###
# 1. The total dimension of the output is L*N0*(p+d), where N0 is the number of acceptance. If this is too large, L should be divided into smaller tasks. 

### In and Out ###
# Input: reference_tables - a list containing two elements, the N*p matrix 'parameters' containing the p-dimension parameter values and the N*d matrix 'summaries' containing the d-dimension summary statistics. observed_summaries - a L*d matrix with rows containing the summary statistics of all observed data sets. tolerance_percentile - a scalar giving the proportion to be accepted in the simulated parameter values. tolerance_value - a scalar giving the value of tolerance level, and only one bewtween this and 'tolerance_percentile' should be given. paral - if true, the algorithms for multiple datasets are run in parallel. ncores - the number of cores used for parallelisation. platform - 'Mac' or 'Win'. packages - If platform is 'Win', this vector of characters specifies packages to be loaded in each core
# Output: an L-element list, each element of which is the output lists of rejABC_OneDataset

	L<-nrow(observed_summaries)
	if(paral==TRUE){
		results_rejABC<-parallel_lapply(1:L,rejABC_MultiDatasets,ncores=ncores,platform=platform,packages=packages,reference_tables=reference_tables,observed_summaries=observed_summaries,tolerance_percentile=tolerance_percentile,tolerance_value=tolerance_value)
		names(results_rejABC)<-NULL
		results_rejABC<-do.call('c',results_rejABC) 
	}
	if(paral==FALSE) results_rejABC<-rejABC_MultiDatasets(1:L,reference_tables,observed_summaries,tolerance_percentile=tolerance_percentile,tolerance_value=tolerance_value)
	return(results_rejABC)
}
