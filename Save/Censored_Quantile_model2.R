

rm(list=ls())

library(mnormt)
#library(metaSEM)
library(gtools)
library(ggplot2)
library(scales)
library(parallel)
dyn.load('lib/sqmc.so')	
source('src/sqmc.R')
source('algo.R')
source('algo_time.R')
source('define_models.R')
numCores <- detectCores() 
#####################################################################################
### Choose the example
#####################################################################################
d<-20
datalength<-10^6
#####################################################################################
### Generate observations
#####################################################################################
#/setting a seed for the RNG
set.seed(9585)
link<-function(X, theta_star){
      return(X%*%matrix(theta_star,nrow=length(theta_star)))
}

#/sample covariates
cov_x<-matrix(rWishart(1,d-1,diag(d-1)),d-1)
cov_x<-solve(cov_x)
#cov_x<-diag(1,d-1)
#/sample covariates
X <- cbind(rep(1,datalength),rmnorm(datalength, rep(0, d-1) , cov_x))
#/true parameter value
theta_star <- c(0, rnorm(d-1,mean=0))
#theta_star <- c(-5, rnorm(d-1,mean=0))
#/sample response variable
sigma_y<-2
Y <- rnorm(datalength, link(X,theta_star),sd=sigma_y)
Y[Y<0]<-0
sum(Y==0)/datalength
# concatenate Y and X into one matrix
observations <- cbind(Y, X)
rm(X,Y)
#####################################################################################
##Choose (h_t) where h_0=0 and h_t=c_h t^{-alpha} when h_t>0
#####################################################################################
alpha<-0.5
h_seq<-(1:datalength)^(-alpha)
#####################################################################################
#Generate (t_p) 
#####################################################################################
scale1<-10
diff_min<-10
f_m<-function(x){return (max(0.5,log(x)))} #need f_m(x)<=x^beta, beta<alpha
f_m<-function(x){return (max(0.5,x^{alpha*0.8}))}
stud_YN<-rep(0,datalength)
t_0<-10
t_1<-1
while(t_1 <= datalength){
  t_1<-t_0+max(diff_min,scale1*floor(f_m(t_0)*log(1/h_seq[t_0])))
  if(t_1 <= datalength) stud_YN[t_1]<-1
  t_0<-t_1
}
#####################################################################################
#Choose when perturbations are applied (4 options)
#####################################################################################
#1.Every steps after START observations
start<-1
pert_times<-start:datalength
stud2_YN<-stud_YN
#4.Every GAP setps after START observations
gap<-5000
start<-5000
pert_times<-seq(start, datalength, gap)
work_times<-c(1:datalength)[stud_YN==1]
stud2_YN<-rep(0,datalength)
for(i in 1:(length(pert_times)-1) ){
     if(max(work_times)>=pert_times[i]){
       work<-min(work_times[work_times>=pert_times[i]])
       if(work<pert_times[i+1] && work>=pert_times[i]){
          stud2_YN[pert_times[i]]<-1
       }else{
         stud2_YN[pert_times[i]]<-0
       }
   }else{
      stud2_YN[pert_times[i]]<-0
   }
}
#####################################################################################
## Choose Sigma_t in the definition of mu_t 
#####################################################################################
c<-sqrt(10)
Sigma_t<-diag(c^2,d)

#####################################################################################
## Build the target
#####################################################################################

target_parameters <- list(quantile=0.5, alpha_val=c(c,alpha),  prior_mean = theta_star+10,
                          prior_covariance =  diag(2, d, d), prior_sd=sqrt(2),
                          h2_t=h_seq^2,T_p=pert_times, student=stud2_YN, Sigma=Sigma_t,
                          nu=50)

target <- list(dimension = d,
               rprior = rprior_quantile,
               dprior = dprior_quantile,
               loglikelihood =loglikelihood_Cquantile, 
               rmutation=rmutation_pi,
               dmutation=dmutation_pi,
               parameters = target_parameters)
#####################################################################################
## Estimation  
#####################################################################################

#Set tuning and target parameters
tuning_parameters <- list(N=10000, K=10000, J=20,  ESS_bound=0.7, resampling=SSP_Resampling) 

target2<-target
target2$parameters$T_p<-datalength+1
##Estimate tilde{pi} and bar{pi}####################
#T_end<-20000#10^5  #datalength
# set.seed(30485)
start_time <- Sys.time()
theta<- online(tuning_parameters, target2  , observations[1:T_end,], numCores-1)
end_time <- Sys.time()
##running time is
print(end_time-start_time)

T_end<-20000
start_time <- Sys.time() 
theta<-Online_censored(tuning_parameters, target , observations[1:T_end,], cl=1, seed=sample(1:10^5,1))
end_time <- Sys.time()
##running time is
print(end_time-start_time)



##Check output####################

#plot ess
plot(theta$ESS/tuning_parameters$N,type='l')

#plot acceptance rate
plot(theta$ACC,type='l')

#mean of tilde{pi}
plot(theta$MEAN[,k],type='l')
abline(h=theta_star[k])



#mean of bar{pi}
bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))


plot(bar_mean[,k],type='l')
abline(h=theta_star[k])


error <-apply((t(bar_mean)-theta_star)^2,2,sum)^{1/2}
error <-apply(abs(t(bar_mean)-theta_star),2,max) 


error2 <-apply((t(theta$MEAN)-theta_star)^2,2,sum)^{1/2}

tr<-1:T_end
nn<-length(tr)
MC<-tr^(-1/2)
C<-error[nn]/MC[nn]
work<-C*MC
points<-seq(1,T_end,length.out=10^5)


plot(log(tr[points]), log(error[points]), type='l')
lines(log(tr[points]), log(work[points]),type='l', col='red')








#####################################################################################
## Estimate tilde{pi}_t and pi_t with same running time 
#####################################################################################

T_end<- datalength
set.seed(30485)
start_time <- Sys.time()
theta<- online(tuning_parameters, target, observations[1:T_end,], cl=numCores-2)
end_time <- Sys.time()
print(end_time-start_time)

time_max<-as.numeric(end_time)-as.numeric(start_time)

target2<-target 
target2$parameters$T_p<-datalength+1

set.seed(30485) 
start_time2 <- Sys.time()
theta_post<- online_time(tuning_parameters, time_max, target2, observations[1:T_end,], cl=numCores-2)
end_time2 <- Sys.time()
##running time is
print(end_time2-start_time2)

target2$parameters$prior_mean<- theta_star
set.seed(30485) 
start_time2 <- Sys.time()
theta_post2<- online_time(tuning_parameters, time_max, target2, observations[1:T_end,], cl=numCores-2)
end_time2 <- Sys.time()
##running time is
print(end_time2-start_time2)


#write.table(theta$MEAN, 'Data/Non_lin_quantile_mean_tildeStart.txt')
#write.table(theta$MEAN2, 'Data/Non_lin_quantile_mean2_tildeStart.txt')
#write.table(error, 'Data/Non_lin_quantile_error_barStart.txt')
#write.table(theta$BAR, 'Data/Non_lin_quantile_theta_barStart.txt')
#write.table(theta$ESS, 'Data/Non_lin_quantile_ess_tildeStart.txt')
#write.table(var_bar, 'Data/Non_lin_quantile_var_barStart.txt')
#write.table(post_var, 'Data/Non_lin_quantile_var_post.txt')
######################################
##Check output####################

plot(theta$ESS/tuning_parameters$N,type='l')

#compare the estimates
plot( theta $MEAN[,k],type='l')
lines(theta_post$MEAN[,k],type='l', col='red')
abline(h=theta_star[k])


bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))

plot(bar_mean[,k],type='l')
#lines(theta_post$MEAN[,k], type='l', col='red')
abline(h=theta_star[k])












