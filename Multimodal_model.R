rm(list=ls())
library(mnormt)
library(parallel)
dyn.load('lib/gpfso.so')	
source('src/function.R')
source('gpfso.R')
source('opsmc.R')
numCores <- detectCores()
#####################################################################################
### Choose the example
#####################################################################################
d<-20  
datalength<- 10^7
#####################################################################################
### Generate observations
#####################################################################################
#/setting a seed for the RNG
set.seed(9585)
link<-function(X, beta_star){
 return( apply(exp(-t(t(X)*(beta_star^2))),1,sum)+X%*%as.matrix(beta_star[seq(length(beta_star),1,-1)]))
}
theta_star <- rep(-1,d)
sigma_y<-2
#/sample covariates
X <- matrix(runif((d)*datalength, min=-1,max=1), datalength,d, byrow=T)
#/sample response variable
Y <- rnorm(datalength, link(X, theta_star),sd=sigma_y)
# concatenate Y and X into one matrix
observations <- cbind(Y, X)
rm(X,Y)
observations<-observations[1:10^6,]
datalength<-10^6
#####################################################################################
## Choose target parameters
#####################################################################################
c_Sigma<-10
alpha<-0.5
#####################################################################################
## Build target
#####################################################################################
rprior<- function(nparticles, target_parameters){
    d<-target_parameters$dim
    return(matrix(runif(d*nparticles, min=-21, max=19),nrow=nparticles, byrow=TRUE))

}
comp_seq<-h_Tp(alpha,datalength)
target_parameters <- list(learning_rate=c(sqrt(c_Sigma),alpha), dim=d, 
                          prior_var = rep(2,d), student=comp_seq$Tp, nu=50)
rm(comp_seq)
target <- list(dimension = d,
               rprior = rprior,
               loglikelihood =loglikelihood_multimodal, 
               parameters = target_parameters)
#####################################################################################
## Tunning parameters
#####################################################################################
tuning_parameters <- list(N=2000,   c_ess=0.7) 
#####################################################################################
## G-PFSO: Escape from local mode
#####################################################################################             
set.seed(30485)
seed_vec=sample(1:10^5,10)
T_end<-30000 
start_time <- Sys.time()
theta<-GPFSO_multi(tuning_parameters, target , observations[1:T_end,], cl=1, seed=seed_vec[8])
end_time <- Sys.time()
print(end_time-start_time)


k<-14
T_end2<-30000
plot(theta$MEAN[100:T_end2,k],type='l')
abline(h=theta_star[k])



##save data for Figure 3(a)
write.table(theta$MEAN[,k],"Data/Multimodal/traj_theta14_d20.txt",col.names=FALSE, row.names=FALSE)
#####################################################################################
## G-PFSO: 10 runs with T=10^6
#####################################################################################
tuning_parameters <- list(N=2000, K=10000,  c_ess=0.7, resampling=SSP_Resampling) 
set.seed(30485)
T_end<-10^6
M<-10
n_points<-10^5
points<-seq(1,T_end,length.out=n_points)
res1<-matrix(0,M,n_points)
res2<-matrix(0,M,n_points)
for(m in 1:M){
   theta<-GPFSO_multi(tuning_parameters, target , observations[1:T_end,], cl=1, seed=sample(1:10^5,1))
   bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
   bar_mean<-bar_mean[points,]
   error <-apply((t(bar_mean)-theta_star)^2,2,sum)^{1/2}
   res1[m,]<-error
   error <-apply(abs(t(bar_mean)-theta_star),2,max) 
   res2[m,]<-error
   rm(theta,bar_mean)
   print(m)
}
tr<-(1:T_end)[points]

errorE<-apply(res1,2,mean)
errorS<-apply(res2,2,mean)

##save data for Figure 3(c)
write.table(errorE,"Data/Multimodal/errorE_d20_M10.txt",col.names=FALSE, row.names=FALSE)
write.table(tr,"Data/Multimodal/tr_d20.txt",col.names=FALSE, row.names=FALSE)
#####################################################################################
## G-PFSO: 100 runs with T=10^5 and with c_sigma=10
#####################################################################################
set.seed(30485)
T_end<-10^5
M<-100
res1<-rep(0,M)
res2<-rep(0,M)
for(m in 1:M){
   theta<-GPFSO_multi(tuning_parameters, target , observations[1:T_end,], cl=1, seed=sample(1:10^5,1))
   bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
   res1[m]<-sum((bar_mean[T_end,]-theta_star)^2)^{1/2} 
   res2[m]<-max(abs(bar_mean[T_end,]-theta_star))
   print(m)
}

##save data for Figure 3(a)
write.table(res2,"Data/Multimodal/boxplotS_d20_M100_T10p5.txt",col.names=FALSE, row.names=FALSE)
#####################################################################################
## G-PFSO: 100 runs with T=10^5 and with c_sigma=3
#####################################################################################
c_Sigma<-3
target$parameters$learning_rate[1]<-sqrt(c_Sigma)
set.seed(30485)
T_end<-10^5
M<-100
res1<-rep(0,M)
res2<-rep(0,M)
for(m in 1:M){
   theta<-GPFSO_multi(tuning_parameters, target , observations[1:T_end,], cl=1, seed=sample(1:10^5,1))
   bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
   res1[m]<-sum((bar_mean[T_end,]-theta_star)^2)^{1/2} 
   res2[m]<-max(abs(bar_mean[T_end,]-theta_star))
   print(m)
   print(res2[m])
}
##save data for Figure 3(a)
write.table(res2,"Data/Multimodal/boxplotS_d20_M100_T10p5_var3.txt",col.names=FALSE, row.names=FALSE)
#####################################################################################
## G-PFSO: 100 runs with T=10^5 and with c_sigma=1
#####################################################################################
c_Sigma<-1
target$parameters$learning_rate[1]<-sqrt(c_Sigma)
set.seed(30485)
T_end<-10^5
M<-100
res1<-rep(0,M)
res2<-rep(0,M)
for(m in 1:M){
   theta<-GPFSO_multi(tuning_parameters, target , observations[1:T_end,], cl=1, seed=sample(1:10^5,1))
   bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
   res1[m]<-sum((bar_mean[T_end,]-theta_star)^2)^{1/2} 
   res2[m]<-max(abs(bar_mean[T_end,]-theta_star))
   print(m)
   print(res2[m])
}
##save data for Figure 3(a)
write.table(res2,"Data/Multimodal/boxplotS_d20_M100_T10p5_var1.txt",col.names=FALSE, row.names=FALSE)
#####################################################################################
## G-PFSO: 100 runs with T=10^5, c_sigma=1 and heavier Student's tails
####################################################################################

c_Sigma<-1
target$parameters$nu<-1.5
target$parameters$learning_rate[1]<-sqrt(c_Sigma)
 

set.seed(30485)
T_end<-10^5
M<-100
res1<-rep(0,M)
res2<-rep(0,M)
for(m in 1:M){
   theta<-GPFSO_multi(tuning_parameters, target , observations[1:T_end,], cl=1, seed=sample(1:10^5,1))
   bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
   res1[m]<-sum((bar_mean[T_end,]-theta_star)^2)^{1/2} 
   res2[m]<-max(abs(bar_mean[T_end,]-theta_star))
   print(m)
   print(res2[m])
}
##save data for Figure 3(a)
write.table(res2,"Data/Multimodal/boxplotS_d20_M100_T10p5_var1_stud.txt",col.names=FALSE, row.names=FALSE)
#####################################################################################
## OPSMC: 100 runs with T=10^5 and iota as in the paper.
#####################################################################################                         
set.seed(30485)
T_end<-10^5
M<-100
res1<-rep(0,M)
res2<-rep(0,M)
for(m in 1:M){
   theta<-OPSMC(tuning_parameters,  target, observations[1:T_end,])
   res1[m]<- sqrt(sum((theta$MEAN[T_end,]-theta_star)^2)) 
   res2[m]<-max(abs(theta$MEAN[T_end,]-theta_star))
   print(m)
}
##save data for Figure 3(a)
write.table(res2,"Data/Multimodal/OPSMC_boxplotS_d20_M100_T10p5.txt",col.names=FALSE, row.names=FALSE)













