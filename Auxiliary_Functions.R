library(matrixStats) 
library(ks) # For function 'kde'
library(spatstat) # For function 'ewcdf'

##########################
# Functions	
##########################
is.error <- function(x) inherits(x, "try-error")


as_matrix<-function(x,nrow=10){
  if(is.vector(x)&nrow==1) return(t(x))
  as.matrix(x)
}

rmNA<-function(x,type=c('matrix_row','matrix_col','vector')){
### What is it for? ###
# Remove NA elements for x, rows containing NA for 'matrix_row', columns cotaining NA for 'matrix_col', NA elements for 'vector',  
  if(type=='matrix_row'){
    NA_ind<-apply(x,1,function(z) sum(is.na(z))>0)
    x<-x[!NA_ind,]
    if(is.vector(x)) x<-matrix(x,nrow=1)
  }
  if(type=='matrix_col'){
    NA_ind<-apply(x,2,function(z) sum(is.na(z))>0) 
    x<-x[,!NA_ind]
    if(is.vector(x)) x<-as.matrix(x)
  }
  if(type=='vector'){x<-x[!is.na(x)]}
  return(x)
}

accept_kernel<-function(v,method='normal'){ 
### What is it for? ###
# The density function used to accept/reject the simulated summary statistic based on the distance measure. Currently 'normal' and 'unif' kernel can be used.

# Input: v - a N*1 vector, the distance between summaries scaled by the bandwidth; 
# Output: a N*1 TRUE/FALSE vector, indicating the acceptance/rejection of each simulation

  N0<-length(v)
  if(method=='normal'){
    acc_prob<-exp(-v^2*pi) # Use N(0,1/(2pi)) kernel density
    if(is.vector(v)) acc_unif<-runif(length(v)) # Generate uniform r.v. for accept/reject
    if(is.matrix(v)) acc_unif<-matrix(runif(length(v)),nrow(v),ncol(v))
    acc_ind<-acc_unif<acc_prob
  }
  if(method=='unif') acc_ind<-(v<1)
  return(acc_ind)
}

rootsearch_1dim<-function(func,val_wanted,interval,grid_size,...){
### What is it for? ###
# A quick and dirty root finding in one dimension for equation func==val_wanted by grid searching.

# Input: func - the target function, with the first arguemnt being the parameter to be searched. val_wanted - a scalar, the target value of func. interval - a 2*1 vector, where the grid is to be drawn. grid_size - a scalar
# Output: a scalar, taking value inside the interval

  grids<-seq(interval[1],interval[2],length=grid_size)
  val_grid<-sapply(grids,func,...)
  minimum<-grids[which.min(abs(val_grid-val_wanted))]
  return(minimum)
}

expected_acc_rate<-function(tol,v,method='normal'){
### What is it for? ###
# Given a vector of summary distance and a value of bandwidth, estimate the expectaion of the acceptance proportion for some acceptance kernel.

# Input: tol - a scalar, the bandwidth. v - a N*1 vector, the distance between summary statistics. method - currently two acceptance kernel can be used, 'normal' and 'unif'.
# Output: a scalar, an estimate of the expected accepted proportion
 
  scaled_v<-v/tol
  if(method=='normal') acc_prob<-exp(-scaled_v^2*pi) # Use N(0,1/(2pi)) kernel density
  if(method=='unif') acc_prob<-(scaled_v<1)
  exp_acc_rate<-sum(acc_prob)/length(v)
  return(exp_acc_rate)
}


HDR_2d_coverred<-function(tested_point,x,y,alpha,weights=NULL,xextend = 0.15, yextend = 0.15){
### What is it for? ###
# Decide whether the given point is coverred by the highest density region(HDR) of the data in two dimension. The HDR is calculated by revising the 'hdr.2d' function in the 'hdrcde' package. The revision considers weighted data. The package 'ks' is required to run the function 'kde' for calculating multivariate-dimension kernel density estimate.

### In and Out ###
# Input: tested_point - a length 2 vector. x,y - vectors of the two-dimension data. alpha - HDR with alpha probability is calculated. weights - a vector with the same length as x/y if the data are weighted. xextend, yextend - two scalars which are proportions of range of x and y; the density is estimated on a grid extended by xextend and yextend beyond the range of x and y.
# Output: 
    in_range<-(tested_point[1]>min(x))*(tested_point[1]<max(x))*(tested_point[2]>min(y))*(tested_point[2]<max(y))
    if(!in_range) {coverred<-FALSE; return(coverred)}
    den <- den.estimate.2d(x,y,weights=weights)
    fxy <- interp.2d(den$x, den$y, den$z, c(tested_point[1],x), c(tested_point[2],y))
    falpha <- quantile(fxy[-1], 1-alpha)
    coverred<-(fxy[1]>=falpha)

    return(coverred)
}

ESS<-function(weights){
### What is it for? ###
# Calculate the effective sample size of importance sampling weights, N/(1+cv(w)^2) where cv is coefficient of variation, i.e. sd(w)/mean(w). An alternative expression of ESS is n*mean(w)^2/mean(w^2)

### In and Out ###
# Input: weights is a N*1 vector 
# Output: a scalar smaller than N

  N<-length(weights)
  N/(1+sd(weights)^2/mean(weights)^2)  
}

signif2<-function(x){
### What is it for? ###
# Round the values of x to have two significant digits 

### In and Out ###
# Input: x is a vector 
# Output: output the rounded vector x
	signif(x,2)
}

savepdf <- function(file, width=16, height=10)
# https://robjhyndman.com/hyndsight/crop-r-figures/
{
  fname <- paste(file,".pdf",sep="")
  pdf(fname, width=width/2.54, height=height/2.54,
      pointsize=10)
  par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
}

savepdf("Gaussian")
# Plotting commands here
dev.off()


cal_MAE<-function(estimations,true_values,normalised=TRUE){
### What is it for? ###
# Calculate the mean absolute error (MAE) with/without normalisation of k p-dimension estimators 

### In and Out ###
# Input: estimations - a k*p matrix with rows containing the estimations. true_values - a k*p matrix with rows containing the true values. normalised - if true, the standard deviations of k true values are used to normalise the MAEs.
# Output: a p*1 vector containing the MAEs/NMAEs for all coordinates

	# sd_of_truth<-sqrt(colVars(true_values)) # a p*1 vector
	MAE<-colMeans(abs(estimations-true_values)) # a p*1 vector
	tmp1<-abs(estimations-true_values)/abs(true_values)
	# tmp_order1<-order(tmp1[,1],decreasing=T)
	# tmp_order2<-order(tmp1[,2],decreasing=T)
	# ind_rm<-c(tmp_order1[1:10],tmp_order2[1:10])
	# NMAE<-colMeans(tmp1[-ind_rm,]) # a p*1 vector
	NMAE<-colMeans(tmp1) # a p*1 vector	
	if(normalised) return(NMAE)
	if(!normalised) return(MAE)
}

rowVars<-function(x) {
### What is it for? ###
# Calculate the variance for each row of x

### In and Out ###
# Input: x - a p*q matrix
# Output: a p*1 vector

  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

wmean<-function(X,weights,normalised=FALSE){
# Input: X - N*d matrix; weights - N*1 vector; normalised - indicate whether the weights sum up to 1
# Output: d*1 weighted mean for X
	if(normalised==FALSE) weights<-weights/sum(weights)
	c(t(X)%*%weights)
}

wvar<-function(X,weights,diagonal=FALSE,normalised=FALSE){
# Input: X - N*d matrix; weights - N*1 vector; normalised - indicate whether the weights sum up to 1
# Output: d*d weighted variance matrix for X
	if(normalised==FALSE) weights<-weights/sum(weights)
	if(diagonal) return(wmean(X^2,weights)-(wmean(X,weights))^2)
	if(!diagonal) return(crossprod(weights*X,X)-tcrossprod(wmean(X,weights)))
}

wcov<-function(X,Y,weights,normalised=FALSE){
# Input: X - N*d matrix; Y - N*p matrix; weights - N*1 vector; normalised - indicate whether the weights sum up to 1
# Output: d*p weighted covriance matrix for X and Y
	if(normalised==FALSE) weights<-weights/sum(weights)
	crossprod(weights*X,Y)-tcrossprod(wmean(X,weights),wmean(Y,weights))
}

wquant<-function(x,q,weights,normalised=FALSE){
# Input: x - N*1 vector; weights - N*1 vector; normalised - indicate whether the weights sum up to 1
# Output: a scalar quantile for the weighted sample
	if(normalised==FALSE) weights<-weights/sum(weights)
	tmp_cdf<-ewcdf(x,weights)
	quantile.ewcdf(tmp_cdf,q)
}

sqdist_matr_matr <- function(A,B){
# Input: A,B - m*p and n*p matrices
# Output: m*n matrix, pairwise squared distance between rows of A and B
 	an = rowSums(A^2)
	bn = rowSums(B^2) 
    
    m = nrow(A)
    n = nrow(B)
 
    tmp = matrix(rep(an, n), nrow=m) 
    tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
    tmp - 2 * tcrossprod(A,B)
}

dist_matr_vec<-function(A,v,Sigma=0){ 
# Input: A - d*N matrix; v - d*1 vector; Sigma - d*d diagonal scale matrix
# Output: the Euclidean distane between each column of A and v
  if(sum(Sigma)==0) Sigma<-diag(rep(1,length(v)))
  distance<-c(sqrt(colSums(Sigma%*%(A-v)^2)))
  return(distance)
}

sqrtmatrix<-function(A){
  if(nrow(A)==1&ncol(A)==1) return(sqrt(A)) 
	a.eig <- eigen(A)
	a.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
	return(a.sqrt)
}

cal_regress_coef<-function(X,Y,weights,lm_fit=TRUE){
# Input: X - N*d matrix; Y - N*p matrix; weights - for weighted least square
# Output: d*p regression coefficient matrix

  p<-ncol(Y)
  if(lm_fit){
    lm_results<-lm(Y~X,weights=weights)
    if(p>1) lm_coef<-lm_results$coef[-1,] else lm_coef<-lm_results$coef[-1]
    if(sum(is.na(lm_coef))>0) lm_coef[is.na(lm_coef)]<-0
    return(lm_coef)
  }
	if(!lm_fit){
    if(is.vector(X)) X<-t(t(X))
    if(is.vector(Y)) Y<-t(t(Y))
    return(solve(wvar(X,weights))%*%wcov(X,Y,weights))
  }
}

cal_distance_between_matrixcol_vector<-function(A,v,Sigma=0){
### What is it for? ###
# Calculate the distance bewteen each column of matrix A and vector v. 

### Specification of Implementation ###
# 1. If Sigma is non-zero, the weighted distance is calculated.

### In and Out ###
# Input: A - a d*N matrix. v - d*1 vector. Sigma - d*d diagonal weight matrix.
# Output: a N*1 vector containing the Euclidean distane between each column of A and v

	if(sum(Sigma)==0) Sigma<-diag(rep(1,length(v)))
	distance<-c(sqrt(colSums(Sigma%*%(A-v)^2)))
	return(distance)
}

draw_densities<-function(density_list,color_vec=1:length(density_list),lty_vec=rep(1,length(density_list)),main=NULL,xlim=NULL,with_legend=FALSE){
### What is it for? ###
# Draw the given kernel density objects in a single plot

### Notes ###

### In and Out ###	
# Input: density_list - a list with elements being p kerndel density objects, and the names will be used in the plot lengend
# Output: A plot object

	p<-length(density_list)
	x_all<-unlist(lapply(density_list,function(obj)obj$x))
	y_all<-unlist(lapply(density_list,function(obj)obj$y))
	if(is.null(xlim)) {xmin<-min(x_all); xmax<-max(x_all)}
  if(!is.null(xlim)) {xmin<-xlim[1]; xmax<-xlim[2]}
	ymin<-min(y_all); ymax<-max(y_all)
	for(plot_i in 1:p){
		if(plot_i==1) plot(density_list[[1]],xlim=c(xmin,xmax),ylim=c(ymin,ymax),col=color_vec[1],lty=lty_vec[1],main=main)
		else lines(density_list[[plot_i]],col=color_vec[plot_i],lty=lty_vec[plot_i])
	}
  if(with_legend) if(!is.null(names(density_list))) legend('topright',names(density_list),lty=lty_vec,col=color_vec)
}

draw_lines<-function(x_vec,lines_matrix,line_col=1:ncol(lines_matrix),line_type=rep(1,ncol(lines_matrix)),line_width=rep(1,ncol(lines_matrix)),main=NULL,sub=NULL,ylim=NULL,xname=NULL,yname=NULL,with_legend=FALSE,legend_position='bottomright',line_names=NULL,with_CI=FALSE,lines_sd_matrix=NULL,min_max_quantile=1){
### What is it for? ###
# Draw lines connecting points in given matrix in a single plot, and also draw the confidence band for each line if the standard deviation for each value in 'lines_matrix' is specified.

### Notes ### 
# 1. For each value X in lines_matrix, if its standard deviation d is specified, then the confidence interval is calculated by (X-1.96*d,X+1.96*d).

### In and Out ###	
# Input: lines_matrix - an n*d matrix where each column is drawn as a line; lines_sd_matrix - an n*d matrix
# Output: A plot object

  CI_line_type<-rep(max(line_type)+1,ncol(lines_matrix))
	d<-ncol(lines_matrix); n<-nrow(lines_matrix)
	y_all<-c(lines_matrix)
	if(is.null(ylim)){ymin<-quantile(y_all,1-min_max_quantile); ymax<-quantile(y_all,min_max_quantile)}
  if(with_CI){
    CI_lower_matrix<-lines_matrix-1.96*lines_sd_matrix 
    CI_upper_matrix<-lines_matrix+1.96*lines_sd_matrix
    value_all<-c(c(y_all),c(CI_lower_matrix),c(CI_upper_matrix))
    if(is.null(ylim)){
      ymin<-quantile(value_all,1-min_max_quantile)
      ymax<-quantile(value_all,min_max_quantile)
    }
  }
  if(!is.null(ylim)){ymin<-ylim[1];ymax<-ylim[2]}
  for(plot_i in 1:d){
    if(plot_i==1) plot(x_vec,lines_matrix[,1],ylim=c(ymin,ymax),col=line_col[1],lty=line_type[1],lwd=line_width[1],main=main,sub=sub,xlab=xname,ylab=yname,type='l')
    else lines(x_vec,lines_matrix[,plot_i],col=line_col[plot_i],lty=line_type[plot_i],lwd=line_width[plot_i])
    if(with_legend){if(plot_i==1) legend(legend_position,line_names,lty=line_type,col=line_col,lwd=line_width)}
    if(with_CI){
      lines(x_vec,CI_lower_matrix[,plot_i],col=line_col[plot_i],lty=CI_line_type[plot_i],lwd=line_width[plot_i])
      lines(x_vec,CI_upper_matrix[,plot_i],col=line_col[plot_i],lty=CI_line_type[plot_i],lwd=line_width[plot_i])
    }
  }
}


Process_outputlist<-function(x,FUNC1,FUNC2,...){
### What is it for? ###
# Given a list x, this apply FUNC1 to each list element, then use FUNC2 to combine the list of FUNC1 output

### Specification of Implementation ###
# 1. It is required that FUNC2 uses only one argument.

### In and Out ###
# Input: x - a list; FUNC1 - function using data structure of element of x as arguments; FUNC2 - function using list as the argument; ... - other arguments for FUNC1
# Output: has the data structure same as output of FUNC2 


	output1<-lapply(x,FUNC1,...) # a length(x)-element list, each element of which is the output of FUNC1
	output2<-do.call(FUNC2,output1)
	return(output2)
}

Monte_Carlo_stability_check<-function(output_processing_func,N,...){
### What is it for? ###
# This is to check whether the Monte Carlo size of the experiment is large enough so that the criterions for comparison have converged.

### Specification of Implementation ###

### In and Out ###
# Input: output_processing_func - a function summarising the output of experiments to a list, and its first argument should be the vector of indices of Monte Carlo sample used to calculate the criterions. N - the Monte Carlo size used in the experiment. ... - arguments for output_processing_func 
# Output: an k*l table where k is the number of checkpoints and l is the number of criterions to be compared.

	number_checkpoints<-10
	checkpoints<-cumsum(rep(N/number_checkpoints,number_checkpoints))
	self_check_table<-NULL
	for(check_i in 1:number_checkpoints){
		self_check_table<-rbind(self_check_table,unlist(output_processing_func(1:checkpoints[check_i],...)))
	}
	rownames(self_check_table)<-paste0('N=',checkpoints)
	return(self_check_table)
}


####################################################
# Functions adapted from other packages
####################################################

ewcdf <- function(x, weights=rep(1/length(x), length(x)))
#
#     ewcdf.R from R package 'spatstat'
#
#     $Revision: 1.12 $  $Date: 2017/08/30 02:04:46 $
#
#  With contributions from Kevin Ummel
#
{
  nw <- length(weights)
  nx <- length(x)
  if(nw == 0) {
    weights <- rep(1/nx, nx)
  } else if(nw == 1) {
    weights <- rep(weights, nx)
  } else if(nw != nx) stopifnot(length(x) == length(weights))
  # remove NA's together
  nbg <- is.na(x) 
  x <- x[!nbg]
  weights <- weights[!nbg]
  n <- length(x)
  if (n < 1)
    stop("'x' must have 1 or more non-missing values")
  stopifnot(all(weights >= 0))
  # sort in increasing order of x value
  ox <- order(x)
  x <- x[ox]
  w <- weights[ox]
  # find jump locations and match
  vals <- sort(unique(x))
  xmatch <- factor(match(x, vals), levels=seq_along(vals))
  # sum weight in each interval
  wmatch <- tapply(w, xmatch, sum)
  wmatch[is.na(wmatch)] <- 0
  cumwt <- cumsum(wmatch)
  # make function
  rval <- approxfun(vals, cumwt,
                    method = "constant", yleft = 0, yright = sum(wmatch),
                    f = 0, ties = "ordered")
  class(rval) <- c("ewcdf", "ecdf", "stepfun", class(rval))
  assign("w", w, envir=environment(rval))
  attr(rval, "call") <- sys.call()
  return(rval)
}

quantile.ewcdf <- function(x, probs=seq(0,1,0.25), names=TRUE, ...,
                           normalise=TRUE, type=1) {
  #trap.extra.arguments(..., .Context="quantile.ewcdf")
  if(!(type %in% c(1,2)))
    stop("Only quantiles of type 1 and 2 are implemented", call.=FALSE)
  env <- environment(x)
  xx <- get("x", envir=env)
  n <- length(xx)
  Fxx <- get("y", envir=env)
  maxFxx <- max(Fxx)
  eps <- 100 * .Machine$double.eps
  if(normalise) {
    Fxx <- Fxx/maxFxx
    maxp <- 1
  } else {
    maxp <- maxFxx
  }
  if(any((p.ok <- !is.na(probs)) &
         (probs/maxp < -eps | probs/maxp > 1 + eps))) {
    allowed <- if(normalise) "[0,1]" else
               paste("permitted range", prange(c(0, maxp)))
    stop(paste("'probs' outside", allowed), call.=FALSE)
  }
  if (na.p <- any(!p.ok)) {
    o.pr <- probs
    probs <- probs[p.ok]
    probs <- pmax(0, pmin(maxp, probs))
  }
  np <- length(probs)
  if (n > 0 && np > 0) {
    qs <- numeric(np)
    if(type == 1) {
      ## right-continuous inverse
      for(k in 1:np) qs[k] <- xx[min(which(Fxx >= probs[k]))]
    } else {
      ## average of left and right continuous
      for(k in 1:np) {
        pk <- probs[k]
        ik <- min(which(Fxx >= probs[k]))
        qs[k] <- if(Fxx[ik] > pk) (xx[ik] + xx[ik-1L])/2 else xx[ik]
      }
    }
  } else {
    qs <- rep(NA_real_, np)
  }
  if (names && np > 0L) {
    dig <- max(2L, getOption("digits"))
    if(normalise) {
      probnames <-
        if(np < 100) formatC(100 * probs, format="fg", width=1, digits=dig) else
        format(100 * probs, trim = TRUE, digits = dig)
      names(qs) <- paste0(probnames, "%")
    } else {
      names(qs) <-
        if(np < 100) formatC(probs, format="fg", width=1, digits=dig) else
        format(probs, trim=TRUE, digits=dig)
    }
  }
  if (na.p) {
    o.pr[p.ok] <- qs
    names(o.pr) <- rep("", length(o.pr))
    names(o.pr)[p.ok] <- names(qs)
    o.pr
  } else qs
}

den.estimate.2d<-function (x, y, weights=NULL, xextend = 0.15, 
    yextend = 0.15) 
{
    N<-length(x)
    if(is.null(weights)) weights<-rep(1,N)
    if(!is.null(weights)) weights<-weights/sum(weights)*N 
    xr <- diff(range(x, na.rm = TRUE))
    yr <- diff(range(y, na.rm = TRUE))
    xr <- c(min(x) - xr * xextend, max(x) + xr * xextend)
    yr <- c(min(y) - yr * yextend, max(y) + yr * yextend)
    X <- cbind(x, y)
    den <- kde(x = X, w = weights, xmin = c(xr[1], yr[1]), 
        xmax = c(xr[2], yr[2]))
    den <- list(x = den$eval.points[[1]], y = den$eval.points[[2]], 
        z = den$estimate)
    return(den)
}

interp.2d<-function (x, y, z, x0, y0) 
#     interp.2d from R package 'hdrcde'
#     
{
    nx <- length(x)
    ny <- length(y)
    n0 <- length(x0)
    z0 <- numeric(length = n0)
    xr <- diff(range(x))
    yr <- diff(range(y))
    xmin <- min(x)
    ymin <- min(y)
    j <- ceiling(((nx - 1) * (x0 - xmin))/xr)
    k <- ceiling(((ny - 1) * (y0 - ymin))/yr)
    j[j == 0] <- 1
    k[k == 0] <- 1
    j[j == nx] <- nx - 1
    k[k == ny] <- ny - 1
    v <- (x0 - x[j])/(x[j + 1] - x[j])
    u <- (y0 - y[k])/(y[k + 1] - y[k])
    AA <- z[cbind(j, k)]
    BB <- z[cbind(j + 1, k)]
    CC <- z[cbind(j + 1, k + 1)]
    DD <- z[cbind(j, k + 1)]
    z0 <- (1 - v) * (1 - u) * AA + v * (1 - u) * BB + v * u * 
        CC + (1 - v) * u * DD
    return(z0)
}

### From https://github.com/cran/ENmisc ###
wtd.boxplot <-
function(x, ...) UseMethod("wtd.boxplot")

wtd.boxplot.default <-
function(x, weights=NULL, ..., range = 1.5, width = NULL, varwidth = FALSE,
         notch = FALSE, outline = TRUE, names, plot = TRUE,
         border = par("fg"), col = NULL, log = "",
         pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5),
         horizontal = FALSE, add = FALSE, at = NULL)
{
args <- list(x, ...)
    namedargs <-
if(!is.null(attributes(args)$names))
    attributes(args)$names != ""
else
    rep(FALSE, length.out = length(args))
    pars <- c(args[namedargs], pars)
    groups <- if(is.list(x)) x else args[!namedargs]
    if (!is.null(weights)){
    if(!is.list(weights)) weights<-list(weights)
    datasize<-sapply(groups,length)
    wtsize<-sapply(weights,length)
    if (length(datasize)!=length(wtsize))
      stop("number of groups for data and weights are different")
    if (any(datasize != wtsize))
        stop("group sizes for data and weights are different")
    groupwts<-weights
    }
    else groupwts<-NULL
    if(0 == (n <- length(groups)))
stop("invalid first argument")
    if(length(class(groups)))
groups <- unclass(groups)
    if(!missing(names))
attr(groups, "names") <- names
    else {
if(is.null(attr(groups, "names")))
    attr(groups, "names") <- 1:n
        names <- attr(groups, "names")
    }
    for(i in 1:n) {
  if(is.null(groupwts[[i]]))
groups[i] <- list(wtd.boxplot.stats(groups[[i]],
   weights=NULL,
   coef=range)) # do.conf=notch)
   else
groups[i] <- list(wtd.boxplot.stats(groups[[i]],
  weights=groupwts[[i]],
   coef=range)) # do.conf=notch)
   }
    stats <- matrix(0,nrow=5,ncol=n)
    conf  <- matrix(0,nrow=2,ncol=n)
    ng <- out <- group <- numeric(0)
    ct <- 1
    for(i in groups) {
stats[,ct] <- i$stats
        conf [,ct] <- i$conf
        ng <- c(ng, i$n)
        if((lo <- length(i$out))) {
            out   <- c(out,i$out)
            group <- c(group, rep.int(ct, lo))
        }
        ct <- ct+1
    }
    z <- list(stats = stats, n = ng, conf = conf, out = out, group = group,
              names = names)
    if(plot) {
bxp(z, width, varwidth = varwidth, notch = notch, log = log,
#            border = border, col = col, pars = pars,
            border = border, boxfill = col, pars = pars,
            outline = outline, horizontal = horizontal, add = add, at = at)
invisible(z)
    }
    else z
}


wtd.boxplot.stats <-
function(x, weights=NULL, coef = 1.5, do.conf=TRUE, do.out=TRUE)
{

    nna <- !is.na(x)
    n <- sum(nna)                       # including +/- Inf
#   stats <- stats::fivenum(x, weights=weights, na.rm = TRUE) # is the new call
# the previous lines needs to be uncommented fot inclusion in the R distribution
# and the next line has to be deleted
    stats <- wtd.fivenum(x, weights=weights, na.rm = TRUE) # is the call for the test version
    iqr <- diff(stats[c(2, 4)])
    if(coef < 0) stop("'coef' must not be negative")
    if(coef == 0)
do.out <- FALSE
    else {                              # coef > 0
out <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
if(any(out[nna])) stats[c(1, 5)] <- range(x[!out], na.rm = TRUE)
    }
    conf <- NULL
    if(do.conf && is.null(weights)) conf <- stats[3] + c(-1.58, 1.58) * iqr / sqrt(n)
    if(do.conf&& !is.null(weights))
    conf <- stats[3] + c(-1.58, 1.58) * iqr * sqrt(sum((weights/sum(weights))^2))
    nn<-ifelse(is.null(weights),n,sum(weights))
    list(stats = stats, n = nn, conf = conf,
 out = if(do.out) x[out & nna] else numeric(0))
}


wtd.fivenum <-
function(x, weights=NULL, na.rm=TRUE)
{
    interpolatedindex<-function(myval,weights){
    indices<-1:length(weights)
    n<-sum(weights)
    weightsbelow<-rep(0,length(weights))
  for (i in 2:length(weights))
        weightsbelow[i] <- weightsbelow[i-1]+weights[i-1]
    weightsabove<-n-weightsbelow-weights
    lowcands<-weightsbelow<myval
    highcands<-weightsabove<n-myval
    (ifelse(any(lowcands),max(indices[lowcands]),1)+
     ifelse(any(highcands),min(indices[highcands]),length(x)))/2
    }
    if (is.null(weights)) weights<-rep(1,length(x))
    if (length(x)>1)
    equalweights<- all((weights[2:length(weights)]-
          weights[1:length(weights)-1])==0)
    else
    equalweights<-TRUE
    xna <- (is.na(x) | weights==0)
    if(na.rm) x <- x[!xna]
    else if(any(xna)) return(rep.int(NA,5))
    sortorder<-order(x)
    x <- x[sortorder]
    weights<-weights[sortorder]
    n <- sum(weights)
    if(n == 0) rep.int(NA,5)
    else {
    if (equalweights){
  d <- c(1, 0.5*floor(0.5*(n+3)), 0.5*(n+1),
       n+1-0.5*floor(0.5*(n+3)), n)
      }
  else {
    if(length(x)>1)
  d<-c(1,sapply(c(0.25*n,0.5*n,0.75*n),
       function(xxx)interpolatedindex(xxx,weights)),
       length(x))
  else
  d<-rep(1,5)
  }
    0.5*(x[floor(d)]+x[ceiling(d)])
  }
}
