library(pbapply) # For 'pblapply', adding progress bar
# library(pbmcapply) # For 'parallel_lapply', adding progress bar to parallel version of 'lapply'
library(snowfall) # For parallelisation
library(parallel) # For parallelisation

##########################
# Functions	
##########################
sf_ExportEnv<-function(package_names){
### What is it for? ###
# When doing parallisation in Windows, load the required packages and export all variables into each core.

	l<-length(package_names)
	for(i in 1:l) sfLibrary(package_names[i],character.only = TRUE)
	sfExportAll()
}

divide<-function(obj,group_no,ini_index=0){
# Input: obj - the set to be divided into groups; group_no - number of groups
# Output: A list with divisions of obj as the elements 	
	obj_ind<-1:length(obj)
	obj_ind_split<-split(obj_ind,ceiling(obj_ind/(length(obj_ind)/group_no))+ini_index)
	return(lapply(obj_ind_split,function(x)obj[x]))
}

parallel_lapply<-function(arguments_for_cores,FUN,ncores=max(detectCores()-1,1),platform='Mac',packages=NULL,export_var='all',progress_bar=FALSE,...){
### What is it for? ###
# Parallelising lapply in Mac or Win using 'mclapply' in <parallel>, 'pbmclapply' in <pbmcapply> or 'parLapply' in <parallel>. 

### Specification of Implementation ###
# 1. In MAC, beware of the memory issue. 'mclapply' can not be used when the memory used in the master session is too large, since it export the entire environment of master session to all worker sessions, and the total memory needed is the memory in master session multiplied the number of cores. 
# 2. For MAC, it requires <pbmcapply> package to show progress bar. But the function 'pbmclapply' doesn't seem to work as it shows the warning that 'mc.cores>1 is not supported on Windows ...'
# 3. In Win, since the user can select variables exported to worker sessions, the memory issue could possibly be avoided.
# 4. For Win, 'pblapply' may be used for parallisation. See the manual for <pbapply> for guide. 
# 5. In Win, if the script is stopped by pressing the 'ESC' key in the middle of parallised running, the clusters need to be closed manually (via task manager). For unknown reason, 'sfStop()' command doesn't work.

### In and Out ###
# Input: arguments_for_cores - an L-element list with each element containing the first argument for the function to run. FUN - function to run, with the first argument being each element in arguments_for_cores. ncores - number of CPU cores to use. platform - 'Win' or 'Mac'. packages - If platform is 'Win', this vector of characters specifies packages to be loaded in each core. ... - arguments to be passed to FUN
# Output: an L-element list, each element of which is the output of one run of FUN.

	if(platform=='Mac'){
		if(progress_bar) results<-pbmclapply(arguments_for_cores,FUN,mc.cores=ncores,...) # This doesn't seem to work in Mac. It shows the warning that 'mc.cores>1 is not supported on Windows ...'
		if(!progress_bar) results<-mclapply(arguments_for_cores,FUN,mc.cores=ncores,...)
	} 
	if(platform=='Win'){
		if(progress_bar){
			cl <- makeCluster(ncores)
			if(export_var=='all'){
				clusterExport(cl, ls(),envir=environment()) # Export all local variables
				clusterExport(cl, ls(envir=.GlobalEnv)) # Export all global variables
			}
			for(i in 1:length(packages)){ # Load all required packages
				clusterExport(cl,"i",envir=environment())
				clusterEvalQ(cl, library(packages[[i]],character.only = TRUE))
			}
			results<-pblapply(arguments_for_cores,FUN,cl=cl,...)
			stopCluster(cl)
		}
		if(!progress_bar){
			sfInit(parallel=TRUE,cpus=ncores)
			if(!is.null(packages)) for(i in 1:length(packages)) sfLibrary(packages[i],character.only = TRUE)
			if(export_var=='all') sfExportAll()
			if(export_var!='all') sfExport(export_var)
			results<-parLapply(sfGetCluster(),arguments_for_cores,FUN,...)
			sfStop()
		}
	}
	return(results)
}

do_recombing<-function(x,out_type){
### What is it for? ###
# This is used for recombining the results after running algorithm on division or in parallel. 

### Specification of Implementation ###
# 1. How the results are recombed depends on the required structure of output. 

### In and Out ###	
# Input: x - an L-element list containing subsets of results to be combined. out_type - choose from 'combined_vectors', 'stacked_matrices', 'fixed_list_combined_vectors', 'fixed_list_stacked_matrices', 'fixed_list_combined_list' and 'stacked_list', and see 'Recombing_Diagram.JPG' for illustration.
# Output: The recombined results with structure depending on 'out_type'.

	if(out_type %in% c('combined_vectors','stacked_list')) x<-do.call(c,x)
	if(out_type=='stacked_matrices') x<-do.call(rbind,x)
	if(out_type %in% c('fixed_list_combined_vectors','fixed_list_stacked_matrices','fixed_list_combined_list')){
		if(out_type=='fixed_list_combined_vectors') combine<-c
		if(out_type=='fixed_list_stacked_matrices') combine<-rbind
		if(out_type=='fixed_list_combined_list') combine<-list
		k<-length(x[[1]]); x_list<-list(0)
		for(list_i in 1:k) x_list[[list_i]]<-do.call(combine,lapply(x,'[[',list_i))
		names(x_list)<-names(x[[1]])
		x<-x_list
	}
	return(x)
}

simpleDR_in_paral<-function(counts,x,FUN,in_type,out_type,...,paral=TRUE,ncores=max(detectCores()-1,1),platform='Mac',packages=NULL){
### What is it for? ###
# It is used to circumvent the memory limitation and also use all cores in parallel. This is done by a two-level division of the input x. First x is divided into smaller parts, each of which uses limited memory to run FUN. Then each part is divided into 'ncore' parts, on which FUN is run in parallel. 

### Specification of Implementation ###
# 1. User only needs to calculate the number of parts that x is divided into, by considering the memory usage of operating each part. Then the division for parallelisation is done automatically. 
# 2. 'pblapply' is used instead of 'lapply' to add progress bar when running.

### In and Out ###	
# Input: counts	- how many parts the rows of x are divided to. x - can be a scalar N, if N objects to be simulated, an N*1 vector, if the vector is to be divided, and an N*p matrix, if its rows are to be divided. FUN - the function to be run. in_type - choose from 'scalar', 'vector' and 'matrix'. out_type - see the arguments of 'do_recombing'. paral - if this is TRUE, FUN in each part is run in parallel. packages - If platform is 'Win', this vector of characters specifies packages to be loaded in each core. ... - arguments to be passed to FUN
# Output: Results after dividing and recombining, with structure depending on 'out_type'.

	if(paral){
		FUN_for_paral<-function(x_part,...){
			if(in_type=='scalar'){
				K<-x_part; tmp<-floor(K/ncores)
				groups<-as.list(c(rep(tmp,ncores-1),K-tmp*(ncores-1)))
				FUN_for_lapply<-function(ind,FUNC,...) FUNC(ind,...)
			}
			if(in_type=='vector'){
				K<-length(x_part)
				groups<-divide(1:K,ncores); names(groups)<-NULL
				FUN_for_lapply<-function(ind,FUNC,...) FUNC(x_part[ind],...)
			}
			if(in_type=='matrix'){ 
				K<-nrow(x_part)
				groups<-divide(1:K,ncores); names(groups)<-NULL
				FUN_for_lapply<-function(ind,FUNC,...) FUNC(x_part[ind,],...)
			}
			results_part<-parallel_lapply(groups,FUN_for_lapply,ncores=ncores,platform=platform,packages=packages,FUNC=FUN,...)
			tmp<-FUN_for_lapply(groups[[1]],FUN,...)
			results_part<-do_recombing(results_part,out_type)
			return(results_part)
		}		
	}
	if(!paral) FUN_for_paral<-FUN

	if(in_type=='scalar'){
		N<-x; tmp<-floor(N/counts)
		groups<-as.list(c(rep(tmp,counts-1),N-tmp*(counts-1))) # Create division indices
		output_test<-FUN(2,...)
		results<-pblapply(groups,function(ind) FUN_for_paral(ind,...)) # apply FUN_for_paral to each group
	}
	if(in_type=='vector'){
		N<-length(x) 
		groups<-divide(1:N,counts); names(groups)<-NULL
		output_test<-FUN(x[1:2],...)
		results<-pblapply(groups,function(ind) FUN_for_paral(x[ind],...))		
	}
	if(in_type=='matrix'){
		N<-nrow(x)
		groups<-divide(1:N,counts); names(groups)<-NULL; names(groups)<-NULL
		output_test<-FUN(x[1:2,],...)
		results<-pblapply(groups,function(ind) FUN_for_paral(x[ind,],...))
	}
	results<-do_recombing(results,out_type)
	return(results)
}
