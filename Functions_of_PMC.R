library(matrixStats) # for 'weightedMad'

##########################
# Functions	
##########################
PMC<-function(N,initial_proposal_func,target_density,proposal_func,iter_max,...){
### What is it for? ###

### Specification of Implementation ###

### In and Out ###	
# Input: 
# initial_proposal_func - a function generating parameter values at the first iteration and calculating the proposal densities; it should use a scalar, the simulation size N, as the input, and give a list containing a N*p matrix, the simulated parameter values, and a N*1 vector, the proposal density, as output.
# target_density - the target density function; it should use a N*p matrix, the parameter values, as the input, and give a N*1 vector, the prior densities, as the output. 
# proposal_func - a function generating parameter values for each PMC iteration and calculating the proposal densities; it should use a scalar, the simulation size N, a N0*p matrix, the latest set of particles, a N0*1 vector, the importance weights of these particles, and a p*p matrix, the , as inputs, and give a list containing a N*p matrix, the simulated parameter values, and a N*1 vector, the proposal density, as output.
# ... - parameters for 'target_density'

# N - a scalar, the simulation size at each iteration.
# iter_max - a scalar, the maximum number of iteration in PMC.


	tmp<-initial_proposal_func(1)$sample
	names_param<-colnames(tmp)
	sd_inflat<-sqrt(2) 

	iter_i<-1; end_iteration<-FALSE; 
	while(!end_iteration){
		if(iter_i==1){
			tmp<-initial_proposal_func(N)
			theta_new<-tmp$sample
			weights_new<-rep(1,N)
		}
		if(iter_i>1){
			tmp<-proposal_func(N,means=proposal_means,SDs=proposal_sd,weights_normalised=proposal_weights)
			theta_new<-tmp$sample
			colnames(theta_new)<-names_param
			propdens_new<-tmp$densities
			targdens_new<-target_density(theta_new,iter_i,...) # For gk example, using synthetic and empirical likelihood requires the additional parameter 'iter_i' in the target_density
			weights_new<-targdens_new/propdens_new
		}
		normweights_new<-weights_new/sum(weights_new)

browser()
targ_dens<-0
target_density(obs_pars[2],iter_i,...)
zquant<-normquantile_pool[[3]][((iter_i-1)*simul_m+1):(iter_i*simul_m),summary_sel]

ELobj_val(obs_pars,simul_m,zquant)

plot(density(proposal_means,weights=proposal_weights),col=2)
points(c(theta_new),targdens_new/mean(weights_new),pch=20)

dev.new()
plot(c(theta_new),targdens_new,pch=20)
# max(targdens_new)/min(targdens_new)


# gkdens<-sumstatloglik_gk(theta_new,iter_i)
# plot(c(theta_new),exp(gkdens-quantile(gkdens,0.9)),pch=20)

# new_density<-list(list(density(sample_posterior[,1])),list(density(sample_posterior[,2])),list(density(sample_posterior[,3])),list(density(sample_posterior[,4])))
# theta_mean_all<-list(parameter_true=obs_pars)
# par(mfrow=c(2,2))
# for(paramToEst_i in 1:length(paramToEst)){
# 	param_i<-paramToEst[paramToEst_i]
# 	new_density[[param_i]][[iter_i+1]]<-density(theta_new[,paramToEst_i],weights=normweights_new)
# 	names(new_density[[param_i]])[iter_i+1]<-paste0('iter',iter_i)
# 	theta_mean_all[[iter_i+1]]<-wmean(theta_new,normweights_new)
# 	names(theta_mean_all)[iter_i+1]<-paste0('iter',iter_i)
# 	draw_densities(new_density[[param_i]][c(1,iter_i:(iter_i+1))],main=names(obs_pars)[param_i])
# 	abline(v=obs_pars[param_i])
# 	for(line_i in iter_i:(iter_i+1)) abline(v=theta_mean_all[[line_i]][param_i],col=line_i-iter_i+2,lty=2)
# 	if(param_i==1) legend('topright',c('Posterior_19sumstat','iter_i','iter_i+1'),lty=c(1,1,1),col=1:3)
# }	
# c
# c

		proposal_means<-theta_new
		proposal_weights<-normweights_new
		proposal_sd<-sd_inflat*sqrtmatrix(wvar(theta_new,normweights_new))

		theta_old<-theta_new
		normweights_old<-normweights_new
		iter_i<-iter_i+1
		if(iter_i>iter_max) end_iteration<-TRUE
	}
	return(list(theta=theta_new,weights=normweights_new))
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
		results<-rbind(results,as.matrix(tmp[prior_satisfied,]))
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

