source('Divide_Recombine_Parallelise.R')
source('Auxiliary_Functions.R')
library(glmnet)

##########################
# Functions	
##########################
local_linear_adjust<-function(accepted_reference_tables,observed_summary,weights_of_sample=NULL){
### What is it for? ###
# Post-adjust the ABC output by local linear regression.  

### Specification of Implementation ###
# 1. If the ABC output is a weighted sample, e.g. from the importance-sampling ABC, the weights are used in the calculation of the regression coeffcients.
# 2. For the Hierarchical normal mean example, the 4th-7th summaries are removed because they are linear combinations of the 1st-3rd summaries.

### In and Out ###
# Input: accepted_reference_tables - a list containing two elements, the N0*p matrix 'parameters' containing the p-dimension parameter values and the N0*d matrix 'summaries' containing the d-dimension summary statistics. observed_summary - a d*1 vector containing the summary statistic of the observed data set. weights_of_sample - a N0*1 vector containing the weights of ABC output if the output is a weighted sample
# Output: a two-element list containing a N0*p matrix 'parameters' and a N0*1 vector 'weights' 

	N<-nrow(accepted_reference_tables$parameters)
	if(is.null(weights_of_sample)) weights_of_sample<-rep(1/N,N)
	regression_coefficents<-cal_regress_coef(accepted_reference_tables$summaries,accepted_reference_tables$parameters,weights_of_sample)
	adjusted_params<-accepted_reference_tables$parameters-t(t(accepted_reference_tables$summaries)-observed_summary)%*%regression_coefficents
	return(list(parameters=adjusted_params,weights=weights_of_sample,coefs=regression_coefficents))
}

regABC_MultiDatasets<-function(running_set_indicators,reference_tables,observed_summaries,tolerance_percentile=NULL,tolerance_value=NULL,weights_of_sample=NULL){
### What is it for? ###
# Run rejection ABC for multiple observed datasets and post-adjust each output by local linear regression.

### Specification of Implementation ###
# 1. For the local linear adjustment, either the function 'abc' in <abc> or the function 'local_linear_adjust' can be used. Currently 'abc' is in use. 

### In and Out ###
# Input: running_set_indicators - a l*1 vector, giving the indicators for which observed data sets the rejection ABC is run for. reference_tables - a list containing two elements, the N*p matrix 'parameters' containing the p-dimension parameter values and the N*d matrix 'summaries' containing the d-dimension summary statistics. observed_summaries - a L*d matrix with rows containing the summary statistics of all observed data sets. tolerance_percentile - a scalar giving the proportion to be accepted in the simulated parameter values. weights_of_sample - weights for weighted sample (currently not in used)
# Output: an l-element list where each element is a two-element list containing the accepted values of parameters and their weights

	l<-length(running_set_indicators)
	if(l>1) observed_summaries_in_use<-observed_summaries[running_set_indicators,]
	if(l==1) observed_summaries_in_use<-t(observed_summaries[running_set_indicators,])
	N<-nrow(reference_tables$parameters)
	if(is.null(weights_of_sample)) weights_of_sample<-rep(1/N,N)
	adjusted_parameters_for_all_set<-list(0)
	for(set_i in 1:l){
		accepted_results<-rejABC_OneDataset(observed_summaries_in_use[set_i,],reference_tables,tolerance_percentile=tolerance_percentile,tolerance_value=tolerance_value)
		if(!is.vector(accepted_results$parameters)) N0<-nrow(accepted_results$parameters)
		if(is.vector(accepted_results$parameters)) N0<-length(accepted_results$parameters)
		accepted_weights<-weights_of_sample[accepted_results$acceptance_ind]
		adjusted_parameters_for_all_set[[set_i]]<-local_linear_adjust(list(parameters=accepted_results$parameters,summaries=accepted_results$summaries),observed_summaries_in_use[set_i,],accepted_weights)
		names(adjusted_parameters_for_all_set)[set_i]<-paste0('set_',running_set_indicators[set_i])
	}
	return(adjusted_parameters_for_all_set)
}

regABC<-function(reference_tables,observed_summaries,tolerance_percentile=NULL,tolerance_value=NULL,weights_of_sample=NULL,paral=FALSE,ncores=if(paral) max(detectCores()-1,1) else 1, platform='Mac',packages=NULL){
### What is it for? ###
# For each of L observed datasets, run regression ABC and output the post-adjusted accepted parameter values.

### Specification of Implementation ###
# 1. The total dimension of the output is L*N0*(p+d), where N0 is the number of acceptance. If this is too large, L should be divided into smaller tasks. 

### In and Out ###
# Input: reference_tables - a list containing two elements, the N*p matrix 'parameters' containing the p-dimension parameter values and the N*d matrix 'summaries' containing the d-dimension summary statistics. observed_summaries - a L*d matrix with rows containing the summary statistics of all observed data sets. tolerance_percentile - a scalar giving the proportion to be accepted in the simulated parameter values. paral - if true, the algorithms for multiple datasets are run in parallel. ncores - the number of cores used for parallelisation. platform - 'Mac' or 'Win'. packages - If platform is 'Win', this vector of characters specifies packages to be loaded in each core
# Output: an L-element list, each element of which is a two-element list containing a N0*p matrix 'parameters' and a N0*d vector 'summaries'

	L<-nrow(observed_summaries)
	if(paral==TRUE){
		results_regABC<-parallel_lapply(1:L,regABC_MultiDatasets,ncores=ncores,platform=platform,packages=packages,reference_tables=reference_tables,observed_summaries=observed_summaries,tolerance_percentile=tolerance_percentile,tolerance_value=tolerance_value,weights_of_sample=weights_of_sample)
		names(results_regABC)<-NULL
		results_regABC<-do.call('c',results_regABC) 
	}
	if(paral==FALSE) results_regABC<-regABC_MultiDatasets(1:L,reference_tables,observed_summaries,tolerance_percentile=tolerance_percentile,tolerance_value=tolerance_value,weights_of_sample=weights_of_sample)
	return(results_regABC)
}

net_coefficients<-function(y,x){
### What is it for? ###
# Calculate the regression coefficients using 'cv.glmnet', i.e. the elastic net with tuning paramters chosen by cross validation

### Specification of Implementation ###

### In and Out ###
# Input: y - n*p matrix; x - n*d matrix
# Output: d*p coefficient matrix

	d<-ncol(x); p<-ncol(y)
	regression_coefficents<-matrix(0,d,p)
	for(p_i in 1:p){
		model_glmnet<-cv.glmnet(x,y[,p_i])
		regression_coefficents[,p_i]<-as.numeric(coef(model_glmnet, s = "lambda.min"))[-1]			
	}
	return(regression_coefficents)
}

#### The option choosing between 'tolerance_percentile' and 'tolerance_value' has not been added in the rest of codes ####

regularizedABC_MultiDatasets<-function(running_set_indicators,reference_tables,observed_summaries,tolerance_percentile,weights_of_sample=NULL){

	p<-ncol(reference_tables[[1]]); d<-ncol(observed_summaries)
	N<-nrow(reference_tables[[1]])
	l<-length(running_set_indicators)
	if(l>1) observed_summaries_in_use<-observed_summaries[running_set_indicators,]
	if(l==1) observed_summaries_in_use<-t(observed_summaries[running_set_indicators,])

	adjusted_parameters_for_all_set<-list(0)
	for(set_i in 1:l){
		accepted_results<-rejABC_OneDataset(observed_summaries_in_use[set_i,],reference_tables,tolerance_percentile)
		N0<-nrow(accepted_results$parameters)
		accepted_weights<-weights_of_sample[accepted_results$acceptance_ind]
		if(!is.null(weights_of_sample)){
			ISResampled_ind<-sample(1:N0,N0,replace=TRUE,prob=accepted_weights)
		    resampled_sum_diff<-t(t(accepted_results$summaries[ISResampled_ind,1:100])-observed_summaries_in_use[set_i,1:100])
			regularised_coefficents<-net_coefficients(as.matrix(accepted_results$parameters[ISResampled_ind,]),resampled_sum_diff)
		} 
	    acc_sum_diff<-t(t(accepted_results$summaries)-observed_summaries_in_use[set_i,])
		# regularised_coefficents<-net_coefficients(accepted_results$parameters,acc_sum_diff)
		adjusted_params<-accepted_results$parameters-acc_sum_diff%*%regularised_coefficents
		adjusted_parameters_for_all_set[[set_i]]<-list(parameters=adjusted_params,weights=accepted_weights)
		names(adjusted_parameters_for_all_set)[set_i]<-paste0('set_',running_set_indicators[set_i])
	}
	return(adjusted_parameters_for_all_set)
}

regularizedABC<-function(reference_tables,observed_summaries,tolerance_percentile,weights_of_sample=NULL){

	L<-nrow(observed_summaries)
	results<-regularizedABC_MultiDatasets(1:L,reference_tables,observed_summaries,tolerance_percentile,weights_of_sample=weights_of_sample)
	return(results)
}

twofoldCV_regABC_MultiDatasets<-function(running_set_indicators,reference_tables,observed_summaries,tolerance_percentile,weights_of_sample=NULL,regularisation=FALSE){

	p<-ncol(reference_tables[[1]])
	N<-nrow(reference_tables[[1]])
	l<-length(running_set_indicators)
	if(l>1) observed_summaries_in_use<-observed_summaries[running_set_indicators,]
	if(l==1) observed_summaries_in_use<-t(observed_summaries[running_set_indicators,])
	adjusted_parameters_for_all_set<-list(0)
	for(set_i in 1:l){
		accepted_results<-rejABC_OneDataset(observed_summaries_in_use[set_i,],reference_tables,tolerance_percentile)
		N0<-nrow(accepted_results$parameters)
		if(is.null(weights_of_sample)) weights_of_sample<-rep(1/N0,N0)
		adj_param_all<-matrix(0,N0,p)

		acc_param1<-accepted_results$parameters[1:(N0/2),]
		acc_sum1<-accepted_results$summaries[1:(N0/2),]
		acc_weights1<-weights_of_sample[1:(N0/2)]/sum(weights_of_sample[1:(N0/2)])
		acc_param2<-accepted_results$parameters[-(1:(N0/2)),]
		acc_sum2<-accepted_results$summaries[-(1:(N0/2)),]
		acc_weights2<-weights_of_sample[-(1:(N0/2))]/sum(weights_of_sample[-(1:(N0/2))])

		if(!regularisation){
			regression_coefficent1<-cal_regress_coef(acc_sum1,acc_param1,acc_weights1)
			regression_coefficent2<-cal_regress_coef(acc_sum2,acc_param2,acc_weights2)			
		}
		if(regularisation){
			acc_sum1_diff<-t(t(acc_sum1)-observed_summaries_in_use[set_i,])
			acc_sum2_diff<-t(t(acc_sum2)-observed_summaries_in_use[set_i,])
			regression_coefficent1<-net_coefficients(acc_param1,acc_sum1_diff)
			regression_coefficent2<-net_coefficients(acc_param2,acc_sum2_diff)			
		}

		acc_sum1_diff<-t(t(acc_sum1)-observed_summaries_in_use[set_i,])
		acc_sum2_diff<-t(t(acc_sum2)-observed_summaries_in_use[set_i,])
		adj_param_all[1:(N0/2),]<-acc_param1-acc_sum1_diff%*%regression_coefficent2
		adj_param_all[-(1:(N0/2)),]<-acc_param2-acc_sum2_diff%*%regression_coefficent1
		adjusted_parameters_for_all_set[[set_i]]<-list(parameters=adj_param_all,weights=weights_of_sample)
		names(adjusted_parameters_for_all_set)[set_i]<-paste0('set_',running_set_indicators[set_i])
	}
	return(adjusted_parameters_for_all_set)
}

twofoldCV_regABC<-function(reference_tables,observed_summaries,tolerance_percentile,regularisation=FALSE){

	L<-nrow(observed_summaries)
	results<-twofoldCV_regABC_MultiDatasets(1:L,reference_tables,observed_summaries,tolerance_percentile,regularisation=regularisation)
	return(results)
}

nfoldCV_regABC_MultiDatasets<-function(running_set_indicators,reference_tables,observed_summaries,tolerance_percentile,weights_of_sample=NULL){

	p<-ncol(reference_tables[[1]])
	N<-nrow(reference_tables[[1]])
	l<-length(running_set_indicators)
	if(l>1) observed_summaries_in_use<-observed_summaries[running_set_indicators,]
	if(l==1) observed_summaries_in_use<-t(observed_summaries[running_set_indicators,])
	adjusted_parameters_for_all_set<-list(0)
	for(set_i in 1:l){
		accepted_results<-rejABC_OneDataset(observed_summaries_in_use[set_i,],reference_tables,tolerance_percentile)
		N0<-nrow(accepted_results$parameters)
		if(is.null(weights_of_sample)) weights_of_sample<-rep(1/N0,N0)

		adj_param_all<-matrix(0,N0,p)
		regression_mean_all<-matrix(0,N0,p)
		for(MC_i in 1:N0){
			acc_param<-accepted_results$parameters[MC_i,]
			acc_sum<-accepted_results$summaries[MC_i,]
			acc_weights_others<-weights_of_sample[-MC_i]/sum(weights_of_sample[-MC_i])
			acc_param_others<-accepted_results$parameters[-MC_i,]
			acc_sum_others<-accepted_results$summaries[-MC_i,]
			regression_coefficent<-cal_regress_coef(acc_sum_others,acc_param_others,acc_weights_others)

			adj_param<-acc_param-t(regression_coefficent)%*%(acc_sum-observed_summaries_in_use[set_i,])
			adj_param_all[MC_i,]<-adj_param
			regression_mean<-colMeans(acc_param_others-t(t(acc_sum_others)-observed_summaries_in_use[set_i,])%*%regression_coefficent)
			regression_mean_all[MC_i,]<-regression_mean
		}
		adjusted_parameters_for_all_set[[set_i]]<-list(parameters=adj_param_all,weights=weights_of_sample,regression_mean=regression_mean_all)
		names(adjusted_parameters_for_all_set)[set_i]<-paste0('set_',running_set_indicators[set_i])
	}
	return(adjusted_parameters_for_all_set)
}

nfoldCV_regABC<-function(reference_tables,observed_summaries,tolerance_percentile){

	L<-nrow(observed_summaries)
	results<-nfoldCV_regABC_MultiDatasets(1:L,reference_tables,observed_summaries,tolerance_percentile)
	return(results)
}
