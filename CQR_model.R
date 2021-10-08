rm(list=ls())
library(mnormt)
library(parallel)
dyn.load('lib/gpfso.so')	
source('src/function.R')
source('gpfso.R')
numCores <- detectCores() 
#####################################################################################
### Choose the example
#####################################################################################
d<-5
datalength<-10^7
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
#/sample covariates
X <- cbind(rep(1,datalength),rmnorm(datalength, rep(0, d-1) , cov_x))
#/true parameter value
theta_star <- c(3, rnorm(d-1,mean=0))
#/sample response variable
sigma_y<-2
Y <- rnorm(datalength, link(X,theta_star),sd=sigma_y)
Y[Y<0]<-0
sum(Y==0)/datalength
# concatenate Y and X into one matrix
observations <- cbind(Y, X)
rm(X,Y)
#####################################################################################
## Choose target parameters
#####################################################################################
c_Sigma<-1
alpha<-0.5
#####################################################################################
## Build target
#####################################################################################
comp_seq<-h_Tp(alpha,datalength)

target_parameters <- list(quantile=0.5, learning_rate=c(sqrt(c_Sigma),alpha), 
                           prior_mean = theta_star+10, prior_var =rep(2, d), 
                             student=comp_seq$Tp, nu=50)
rm(comp_seq)                      
target <- list(dimension = d, parameters = target_parameters)

theta_star[1]<-  3+ sigma_y*qnorm(target$parameters$quantile)

#####################################################################################
## Tunning parameters
#####################################################################################
tuning_parameters <- list(N=1000,  c_ess=0.7) 
##########################################################################################################################
#G-PFSO: 20 runs with T=10^7 and quantile=0.5
##########################################################################################################################
target$parameters$quantile<-0.5
theta_star[1]<-  3
set.seed(30485)
T_end<-10^7
M<-20
points<-seq(1,T_end,length.out=10^5)
res1<-matrix(0,M,10^5)
res2<-matrix(0,M,10^5)
for(m in 1:M){
   theta<-GPFSO_cqr(tuning_parameters, target , observations[1:T_end,], cl=1, seed=sample(1:10^5,1))
   bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
   error <-apply((t(bar_mean)-theta_star)^2,2,sum)^{1/2}
   res1[m,]<-error[points]
   error <-apply(abs(t(bar_mean)-theta_star),2,max) 
   res2[m,]<-error[points]
   print(m)
}
errorE<-apply(res1,2,mean)
errorS<-apply(res2,2,mean)
tr<-(1:T_end)[points]

#save data for Figure 2(b)
write.table(errorE,"Data/CQR/errorE_d5_q5_M20.txt",col.names=FALSE, row.names=FALSE)
write.table(tr,"Data/CQR/tr_d5.txt",col.names=FALSE, row.names=FALSE)

##########################################################################################################################
#Jittering: 20 runs with T=10^7 and quantile=0.5
##########################################################################################################################
target$parameters$quantile<-0.5
theta_star[1]<-  3
set.seed(30485)
T_end<-10^7
M<-20
points<-seq(1,T_end,length.out=10^5)
res1<-matrix(0,M,10^5)
res2<-matrix(0,M,10^5)
for(m in 1:M){
   theta<-jittering_cqr(tuning_parameters, target, iota=1, observations[1:T_end,], cl=1, seed=sample(1:10^5,1))
   bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
   error <-apply((t(bar_mean)-theta_star)^2,2,sum)^{1/2}
   res1[m,]<-error[points]
   error <-apply(abs(t(bar_mean)-theta_star),2,max) 
   res2[m,]<-error[points]
   print(m)
}
errorE<-apply(res1,2,mean)
errorS<-apply(res2,2,mean)

#save data for Figure 2(b)
write.table(errorE,"Data/CQR/errorE_d5_q5_jit_M20.txt",col.names=FALSE, row.names=FALSE)

##########################################################################################################################
#G-PFSO: 20 runs with T=10^7 and quantile=0.99
##########################################################################################################################
target$parameters$quantile<-0.99
theta_star[1]<-  3+ sigma_y*qnorm(target$parameters$quantile)
set.seed(30485)
T_end<-10^7
M<-20
points<-seq(1,T_end,length.out=10^5)
res1<-matrix(0,M,10^5)
res2<-matrix(0,M,10^5)
for(m in 1:M){
   theta<-GPFSO_cqr(tuning_parameters, target , observations[1:T_end,], cl=1, seed=sample(1:10^5,1))
   bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
   error <-apply((t(bar_mean)-theta_star)^2,2,sum)^{1/2}
   res1[m,]<-error[points]
   error <-apply(abs(t(bar_mean)-theta_star),2,max) 
   res2[m,]<-error[points]
   print(m)
}
errorE<-apply(res1,2,mean)
errorS<-apply(res2,2,mean)
tr<-(1:T_end)[points]

#save data for Figure 2(c)
write.table(errorE,"Data/CQR/errorE_d5_q99_M20.txt",col.names=FALSE, row.names=FALSE)

##########################################################################################################################
#G-PFSO: 20 runs with T=10^7, quantile=0.99 and alpha=0.3
##########################################################################################################################
target$parameters$quantile<-0.99
theta_star[1]<-  3+ sigma_y*qnorm(target$parameters$quantile)
alpha<-0.3
comp_seq<-h_Tp(alpha,datalength)

target_parameters <- list(quantile=0.99, learning_rate=c(c_Sigma,alpha),  prior_mean = theta_star+10,
                          prior_var = rep(2, d), student=comp_seq$Tp, nu=50)
rm(comp_seq)
target$parameters<-target_parameters
set.seed(30485)
T_end<-10^7
M<-20
points<-seq(1,T_end,length.out=10^5)
res1<-matrix(0,M,10^5)
res2<-matrix(0,M,10^5)
for(m in 1:M){
   theta<-GPFSO_cqr(tuning_parameters, target , observations[1:T_end,], cl=1, seed=sample(1:10^5,1))
   bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
   error <-apply((t(bar_mean)-theta_star)^2,2,sum)^{1/2}
   res1[m,]<-error[points]
   error <-apply(abs(t(bar_mean)-theta_star),2,max) 
   res2[m,]<-error[points]
   print(m)
}
errorE<-apply(res1,2,mean)
errorS<-apply(res2,2,mean)
tr<-(1:T_end)[points]

#save data for Figure 2(c)
write.table(errorE,"Data/CQR/errorE_d5_q99_M20_a3.txt",col.names=FALSE, row.names=FALSE)

##########################################################################################################################
#ADA GRAD: 1000 runs with T=10^7 and quantile=0.5
##########################################################################################################################
quantile<-0.5
target$parameters$quantile<-quantile
theta_star[1]<- 3+ sigma_y*qnorm(target$parameters$quantile)

T_end<-10^7
M<-1000
val<-matrix(0,M,d)
c_val<-0.5
set.seed(30485)
for(m in 1:M){
       val[m,]<- rnorm(d,theta_star+10,sd=sqrt(2))
}
est<- ADA_cqr(observations[1:T_end,], quantile, c=c_val,  theta0=val)

errorE<- sqrt(apply( (t(est$TILDE)-theta_star)^2,2,sum))
errorS<- (apply( abs(t(est$TILDE)-theta_star),2,max))

#save data for Figure 2(a)
write.table(errorS,"Data/CQR/ADA_BoxplotS_d5_q50.txt",col.names=FALSE, row.names=FALSE)

##########################################################################################################################
#ADA GRAD: 1000 runs with T=10^7 and quantile=0.99
##########################################################################################################################

quantile<-0.99
target$parameters$quantile<-quantile
theta_star[1]<-  3+ sigma_y*qnorm(target$parameters$quantile)

T_end<-10^7
M<-1000
val<-matrix(0,M,d)
c_val<-0.5
set.seed(30485)
for(m in 1:M){
       val[m,]<- rnorm(d,theta_star+10,sd=sqrt(2))
}
est<- ADA_cqr(observations[1:T_end,], quantile, c=c_val,  theta0=val)

errorE<- sqrt(apply( (t(est$TILDE)-theta_star)^2,2,sum))
errorS<- (apply( abs(t(est$TILDE)-theta_star),2,max))

#save data for Figure 2(a)
write.table(errorS,"Data/CQR/ADA_BoxplotS_d5_q99.txt",col.names=FALSE, row.names=FALSE)



