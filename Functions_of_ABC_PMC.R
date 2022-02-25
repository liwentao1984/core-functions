library(matrixStats) # for 'weightedMad'

##########################
# ABC-PMC Functions	
##########################
ABC_PMC<-function(observations=NULL,observed_summary=NULL,simulator,initial_proposal_func,prior_density,proposal_func,N,bandwidth_percentile,iter_max=100,accept_method='unif',summary_normalising=TRUE,iter_stopnormalising=10,divide_counts=1,...){
### What is it for? ###
# Apply popultion Monte Carlo type of ABC algorithm on given observed data/summary and simulator. See the 'ABC-PMC algorithm' note on 10 Mar, 2018 for the algorithm. 

### Specification of Implementation ###

### In and Out ###	
# Input: 
# observations - a n*M matrix containing the M-dimension observed data, with row being each data point. This is not necessary if observed data is not needed for simulating pseudo data set and calculating summary statistic. 
# observed_summary - a d*1 vector containing the summary statistic of the observed data set. 
# simulator - the function simulating the pseudo summary statistics; it should use a N*p matrix, the parameter values, as the input, and give a N*d matrix, simulated summaries, as the output. 
# initial_proposal_func - a function generating parameter values at the first iteration and calculating the proposal densities; it should use a scalar, the simulation size N, as the input, and give a list containing a N*p matrix, the simulated parameter values, and a N*1 vector, the proposal density, as output.
# prior_density - the prior density function; it should use a N*p matrix, the parameter values, as the input, and give a N*1 vector, the prior densities, as the output. 
# proposal_func - a function generating parameter values for each PMC iteration and calculating the proposal densities; it should use a scalar, the simulation size N, a N0*p matrix, the latest set of particles, a N0*1 vector, the importance weights of these particles, and a p*p matrix, the , as inputs, and give a list containing a N*p matrix, the simulated parameter values, and a N*1 vector, the proposal density, as output.

# N - a scalar, the simulation size at each iteration.
# bandwidth_percentile - a scalar, the percentile used to calculate the bandwidth at each iteration.
# iter_max - a scalar, the maximum number of iteration in PMC.
# accept_method - the probability function used to accept/reject pseudo datasets; it should use a N*1 vector, the distances between summaries scaled by the bandwidth, as input, and give a N*1 TRUE/FALSE vector as output.
# summary_normalising - if TRUE, the summary statistics are standardised by their weighted median absolute deviation (MAD). 
# divide_counts - the number of division when dividing-recombining strategy is applied on 

# Output: 

	# Settings
	for(par_i in 1:(p+2)) dev.new()
	tmp<-initial_proposal_func(1)$sample
	names_param<-colnames(tmp)
	p<-length(c(tmp))
	d<-length(observed_summary)
	sd_inflat<-sqrt(2) 

	N_acc<-0; theta_acc<-NULL; weights_acc<-NULL; sumstat_acc<-NULL; linadj_theta_acc<-NULL
	divide_counts<-1
	sumstat_all_sd<-rep(1,d)
	stat_normalisation<-diag(1/sumstat_all_sd)
	summary_stat_target<-observed_summary

	# Population Monte Carlo
	iter_i<-1; end_iteration<-FALSE; 
	while(!end_iteration){
		if(iter_i==1){
	# acc_density<-list(list(0),list(0),list(0))
	# theta_mean_all<-list(parameter_true=parameter_true)		
			tmp<-initial_proposal_func(N)
			theta_all<-tmp$sample
			weights_all<-tmp$densities
			sumstat_all<-simulator(theta_all)
		}
		if(iter_i>1){
			N_newsample<-N-nrow(theta_acc)
			tmp<-proposal_func(N_newsample,means=proposal_means,SDs=proposal_sd,weights=proposal_weights)
			theta_new<-tmp$sample
			propdens_new<-tmp$densities
			priordens_new<-prior_density(theta_new)
			weights_new<-priordens_new/propdens_new
			sumstat_new<-simulator(theta_new)
			theta_all<-rbind(theta_acc,theta_new)
			weights_all<-c(weights_acc,weights_new)
			sumstat_all<-rbind(sumstat_acc,sumstat_new)
		}
		normweights_all<-weights_all/sum(weights_all)
		if(summary_normalising){
			sumstat_all_sd<-apply(sumstat_all,2,weightedMad,w=normweights_all)
			sumstat_sd_zero_ind<-(sumstat_all_sd<1e-8)
			if(sum(sumstat_sd_zero_ind)>0) sumstat_all_sd[sumstat_sd_zero_ind]<-1
			stat_normalisation<-diag(1/sumstat_all_sd)
		}
		obs_summary_normalised<-observed_summary/sumstat_all_sd
		sumstat_all_norm<-t(t(sumstat_all)/sumstat_all_sd)
		MC_dist<-dist_matr_vec(t(sumstat_all),c(observed_summary),stat_normalisation^2) # N*1 vector

		MC_dist[is.na(MC_dist)]<-Inf
		MC_dist_finite<-MC_dist
		if(sum(!is.finite(MC_dist))>0) MC_dist_finite<-MC_dist[is.finite(MC_dist)]

		bandwidth<-rootsearch_1dim(expected_acc_rate,val_wanted=bandwidth_percentile,interval=c(min(MC_dist_finite),quantile(MC_dist_finite,.9)),grid_size=1000,v=MC_dist,method=accept_method)
		acc_ind<-accept_kernel(MC_dist/bandwidth,method=accept_method)
		N_acc<-sum(acc_ind)
		N_acc_ind<-(1:N)[acc_ind]
		theta_acc<-as.matrix(theta_all[N_acc_ind,])
		colnames(theta_acc)<-names_param
		sumstat_acc<-sumstat_all[N_acc_ind,]
		sumstat_norm_acc<-sumstat_all_norm[N_acc_ind,]
		weights_acc<-weights_all[N_acc_ind]
		normweights_acc<-weights_acc/sum(weights_acc)

		proposal_means<-theta_acc
		proposal_weights<-normweights_acc 
		proposal_sd<-sd_inflat*sqrtmatrix(wvar(theta_acc,normweights_acc))

		# Plotting PMC progress
		accepted_reference_tables<-list(parameters=theta_acc,summaries=sumstat_norm_acc)
		resamp_ind<-sample(1:N_acc,size=N_acc,replace=T,prob=normweights_acc) # Multinomial resampling
		reg_adj_param<-local_linear_adjust(accepted_reference_tables,obs_summary_normalised,weights_of_sample=normweights_acc)$parameters # Regression-adjusted sample
		aa<-dev.list()
		if(iter_i==1){
			dev.set(which=aa[1])
			plot(1,1,xlim=c(1,iter_max),ylim=c(1,iter_max),pch=20,xlab='iter',ylab='iter',main='PMC_iteration')
			for(p_i in 1:p){
				dev.set(which=aa[p_i+1])
				plot(rep(1,N_acc),theta_acc[resamp_ind,p_i],xlim=c(1,iter_max),ylim=c(min(c(params_true[p_i],theta_acc[,p_i])),max(c(params_true[p_i],theta_acc[,p_i]))),pch=20,ylab=paste0('param',p_i),xlab='iter',main=paste0('Accepted parameter',p_i))
				abline(h=params_true[p_i])
			}
			dev.set(which=aa[p+2])
			plot(1,bandwidth,xlim=c(1,iter_max),ylim=c(0,2*bandwidth),pch=20,xlab='iter',ylab='tolerance',main='PMC_tolerance')
		}
		if(iter_i>1){
			dev.set(which=aa[1])
			points(iter_i,iter_i,pch=20)
			for(p_i in 1:p){
				dev.set(which=aa[p_i+1])
				points(rep(iter_i,N_acc),reg_adj_param[resamp_ind,p_i],pch=20)
			}
			dev.set(which=aa[p+2])
			points(iter_i,bandwidth,pch=20)
			if(iter_i==iter_stopnormalising){ # Plot the range of sumstats distance generated from the truth
				dist_obs_all<-dist_matr_vec(t(obs_sumstat_others),obs_sumstat,stat_normalisation^2)
				abline(h=quantile(dist_obs_all,acc_rate),lty=2,col=4)
				abline(h=quantile(dist_obs_all,2*acc_rate),lty=2,col=4)
			}
		}
		
		# par(mfrow=c(2,2))
		# for(param_i in 1:3){
		# 	acc_density[[param_i]][[iter_i]]<-density(theta_acc[,param_i],weights=normweights_acc)
		# 	names(acc_density[[param_i]])[iter_i]<-paste0('iter',iter_i)
		# 	theta_mean_all[[iter_i+1]]<-wmean(theta_acc,normweights_acc)
		# 	names(theta_mean_all)[iter_i+1]<-paste0('iter',iter_i)
		# 	draw_densities(acc_density[[param_i]],main=names(parameter_true)[param_i])
		# 	abline(v=parameter_true[param_i])
		# 	for(line_i in 1:iter_i) abline(v=theta_mean_all[[line_i]][param_i],col=line_i,lty=2)
		# }	
		iter_i<-iter_i+1
		if(iter_i>iter_max) end_iteration<-TRUE
		if(iter_i>=iter_stopnormalising) summary_normalising<-FALSE
		# flush.console() 
	}
	return(list(theta=theta_acc,weights=normweights_acc))
}

nonparamproposal<-function(N,means,SDs,sample_weights,method='normal',check_prior){
### What is it for? ###
# It contains the parameter value simulator and its density function.

# Input: N - the size of simulation. means,sample_weights - a N0*p matrix and a N0*1 vector, forming a weighted sample, and the sample_weights should be nonnegative. SDs - a p*p matrix, the standard deviation matrix of each component of the mixture proposal. 
# Output: a list containing a N*p matrix, the simulated parameter values, and a N*1 vector, the proposal density

	weights_normalised<-sample_weights/sum(sample_weights)
	sample_and_constant<-rnonparamproposal(N,means,SDs,weights_normalised,method=method,check_prior=check_prior,normalising_constant_output=TRUE)
	sample_simulated<-sample_and_constant$sample
	normalising_constant_est<-sample_and_constant$normalising_constant
	densities<-dnonparamproposal(sample_simulated,means,SDs,weights_normalised,method=method,constant=normalising_constant_est)
	return(list(sample=sample_simulated,densities=densities))
} 

rnonparamproposal<-function(N,means,SDs,weights_normalised,method='normal',check_prior,normalising_constant_output=FALSE){
### What is it for? ###
# A nonparametric-type proposal distribution constructed using a given sample, where the sample are the means and the mixture weights are from KDE. It can be chosen whether it's a mixture of normal or t_4 distribution. The distribution is truncated in the domain of prior distribution, so only simulated values falling in of the domain are accepted. The normalising constant due to this truncation is estimated by the acceptance proportion. 

# Input: N - the size of simulation. means,weights_normalised - a N0*p matrix and a N0*1 vector, forming a weighted sample, and the weights_normalised should sum to one. SDs - a p*p matrix, the standard deviation matrix of each component of the mixture proposal. method - currently two nonparametric kernels, 'normal' and 't', can be used. check_prior - a function giving TRUE/FALSE, checking whether a parameter value is in the domain of the prior. normalising_constant_output - if TRUE, the normalising constant will be output
# Output: a list containing N*p matrix of MC sample and a scalar if the normalising constant estiamted is output

	N0<-length(weights_normalised); p<-nrow(SDs)	
	N_sample<-0; N_need<-N-N_sample; N_total<-0
	results<-NULL; df_t<-4
	while(N_need>0){
		N_total<-N_total+N_need
		group_size<-rep(0,N0)
		mixture_ind<-sample(1:N0,N_need,weights_normalised,replace=T)
		mixture_ind<-table(mixture_ind); group_size[as.numeric(names(mixture_ind))]<-mixture_ind
		if(method=='normal') tmp<-matrix(rnorm(N_need*p),nrow=N_need)%*%SDs+means[rep(1:N0,group_size),] # N*p matrix
		if(method=='t') tmp<-matrix(rt(N_need*p,df=df_t),nrow=N_need)%*%SDs+means[rep(1:N0,group_size),] # N*p matrix
		prior_satisfied<-check_prior(tmp) # N*1 vector
		results<-rbind(results,tmp[prior_satisfied,])
		N_sample<-N_sample+sum(prior_satisfied==1)
		N_need<-N-N_sample
	}
	normalising_constant<-N_sample/N_total
	return(list(sample=results,normalising_constant=normalising_constant))
}

dnonparamproposal<-function(x,means,SDs,weights_normalised,method='normal',constant=1){
### What is it for? ###
# Calculating the density of a nonparametric-type proposal distribution constructed using a given sample, where the sample are the means and the mixture weights are from KDE. It can be chosen whether it's a mixture of normal or t_4 distribution. The distribution is truncated in the domain of prior distribution, and the normalising constant of this truncated distribution can be estimated by the acceptance proportion when simulating from it. 

# Input: x - a N*p matrix, the parameter values; means,weights_normalised - a N0*p matrix and a N0*1 vector, forming a weighted sample, and the weights_normalised should sum to one. SDs - a p*p matrix, the standard deviation matrix of each component of the mixture proposal. method - currently two nonparametric kernels, 'normal' and 't', can be used. constant - normalising constant of the mixture proposal truncated in the prior domain.
# Output: N*1 vector of proposal density

	sd_inv<-solve(SDs)
	p<-ncol(x)
	tmp<-t(sqdist_matr_matr(x%*%sd_inv,means%*%sd_inv)) # N0*N pairwise squared distance matrices between rows
	if(method=='normal'){
		tmp<-exp(-tmp/2) # N0*N matrix of unnormalised normal density
		norm_constant<-det(sd_inv)/(2*pi)^(p/2) 
		tmp<-colSums(tmp*weights_normalised)*norm_constant # N*1 vector of mixture normal density		
	}
	if(method=='t'){
		df_t<-4
		t_constant<-(df_t^(p/2)*pi^(p/2)*det(sd)*gamma(df_t/2)/gamma((df_t+p)/2))^(-1)
		tmp<-(1+tmp/df_t)^(-(df_t+p)/2)
		tmp<-colSums(tmp*weights_normalised)*t_constant
	}
	return(tmp/constant)
}


