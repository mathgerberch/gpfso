rm(list=ls())
library(mnormt)
library(parallel)
dyn.load('lib/gpfso.so')	
source('src/function.R')
source('gpfso.R')
numCores <- detectCores()
#####################################################################################
### Simulation set-up
#####################################################################################
datalength<- 10^6
dx<-4
a1<-1:dx
m1<-1:dx
m2<-1:dx
s1<-1:dx
s2<-1:dx
d<-length(c(a1,m1,m2,s1,s2))

#sample covariates
cov_x<-diag(1,dx-1)
set.seed(30485)
X <- cbind(rep(1,datalength),rmnorm(datalength, rep(0, dx-1) , cov_x)) 
 
#true parameter
theta_star<- c(c(1,0.1,0.1,-0.1), c(1, 1, 1,-1), c(-1,1,1,1), c(0,1,1,1), c(0.5,-1,-1,1))

#sample the response variables
proba<- matrix(X[,a1], ncol=length(a1))%*%matrix(theta_star[1:length(a1)],ncol=1)
proba<- exp(-proba)/(1+exp(-proba))
mu1<- matrix(X[,m1], ncol=length(m1))%*%matrix(theta_star[(length(a1)+1):length(c(a1,m1))],ncol=1)
mu2<- matrix(X[,m2], ncol=length(m2))%*%matrix(theta_star[(length(c(a1,m1))+1):length(c(a1,m1,m2))],ncol=1)
sigma1<- exp(-matrix(X[,s1], ncol=length(s1))%*%matrix(theta_star[(length(c(a1,m1,m2))+1):length(c(a1,m1,m2,s1))],ncol=1))
sigma2<- exp(-matrix(X[,s2], ncol=length(s2))%*%matrix(theta_star[(length(c(a1,m1,m2,s1))+1):length(c(a1,m1,m2,s1,s2))],ncol=1))
ind<-rbinom(datalength, 1,proba)
Y<-rep(0,datalength)
Y[ind==1]<- rnorm(sum(ind==1), mu1[ind==1], sigma1[ind==1])
Y[ind==0]<- rnorm(sum(ind==0), mu2[ind==0], sigma2[ind==0])
observations <- cbind(Y, X)
rm(proba,mu1,mu2,sigma1,sigma2,ind,X,Y)

#####################################################################################
## Choose target parameters
#####################################################################################
c_Sigma<-1
alpha<-0.5
#####################################################################################
## Build target
#####################################################################################
comp_seq<-h_Tp(alpha,datalength)
target_parameters <- list(mu=c(1,rep(0,d-1)), sd=rep(1,d), prob=a1, mu1=m1, mu2=m2,
                          sig1=s1, sig2=s2, student=comp_seq$Tp, nu=50, 
                          learning_rate=c(sqrt(c_Sigma),alpha))

rm(comp_seq)
target <- list(dimension = d, parameters = target_parameters)            
#####################################################################################
##Generate seeds
#####################################################################################
set.seed(30485) 
M<-100
seed_vec<-sample(1:10^5,M)
#####################################################################################
## G-PFSO: 100 runs with T=10^5 iterations
##################################################################################### 
T_end<- 10^5
T_end2<-2*10^4
res1<-rep(0,M)
res2<-rep(0,M)
res3<-rep(0,M)
res4<-rep(0,M)
res11<-rep(0,M)
res22<-rep(0,M)
res33<-rep(0,M)
res44<-rep(0,M)
for(m in 1:M){
   theta<-GPFSO_mixture(N=5000, target, observations[1:T_end,], cl=1, seed=seed_vec[m])
   bar_mean<-apply(theta$MEAN,2,mean)
   bar_mean2<-apply(theta$MEAN[1:T_end2,],2,mean)
   est<-theta$MEAN[T_end,]
   est2<-theta$MEAN[T_end2,]
   res1[m]<- sqrt( sum(bar_mean-theta_star)^2)
   res2[m]<- max( abs(bar_mean-theta_star))
   res3[m]<- sqrt( sum(est-theta_star)^2)
   res4[m]<- max(abs(est-theta_star))
   res11[m]<- sqrt( sum(bar_mean2-theta_star)^2)
   res22[m]<- max( abs(bar_mean2-theta_star))
   res33[m]<- sqrt( sum(est2-theta_star)^2)
   res44[m]<- max(abs(est2-theta_star))
   print(m)
   print(cbind(res44[m],res4[m]))
}

##Save data for Figure 5(a)
#write.table(theta_star,"Data/Mixture/theta_star.txt",col.names=FALSE, row.names=FALSE)
#write.table(res2,"Data/Mixture/BoxplotS_M100_T10p5.txt",col.names=FALSE, row.names=FALSE)
#write.table(res4,"Data/Mixture/tBoxplotS_M100_T10p5.txt",col.names=FALSE, row.names=FALSE)
#write.table(res22,"Data/Mixture/BoxplotS_M100.txt",col.names=FALSE, row.names=FALSE)
#write.table(res44,"Data/Mixture/tBoxplotS_M100.txt",col.names=FALSE, row.names=FALSE)


#Example of run where G-PFSO Escapes from local mode
###############################################
T_end<-10^5
start_time <- Sys.time()
theta<-GPFSO_mixture(N=5000, target, observations[1:T_end,], cl=1, seed=seed_vec[14])
end_time <- Sys.time()
##running time is
print(end_time-start_time)
k<-3
plot(theta$MEAN[1000:T_end,k],type='l')
abline(h=theta_star[k])

k<-5
plot(theta$MEAN[1:50000,k],type='l')
abline(h=theta_star[k])

#Save data for Figure 4(b)
#write.table(theta$MEAN[,3],"Data/Mixture/traj_3.txt",col.names=FALSE, row.names=FALSE)
#write.table(theta$MEAN[,5],"Data/Mixture/traj_5.txt",col.names=FALSE, row.names=FALSE)


T_end<-10^5
start_time <- Sys.time()
theta<-GPFSO_mixture(N=5000, target, observations[1:T_end,], cl=1, seed=seed_vec[25])
end_time <- Sys.time()

k<-6
plot(theta$MEAN[,k],type='l')
abline(h=theta_star[k])

#Save data for Figure 4(b)
#write.table(theta$MEAN[,6],"Data/Mixture/traj_6.txt",col.names=FALSE, row.names=FALSE)


 
 
#####################################################################################
#Jittering estimation: 100 runs with T=10^5 iterations
##################################################################################### 
T_end<-10^5
T_end2<-2*10^4
M<-100
res1<-rep(0,M)
res2<-rep(0,M)
res3<-rep(0,M)
res4<-rep(0,M)
res11<-rep(0,M)
res22<-rep(0,M)
res33<-rep(0,M)
res44<-rep(0,M)
est_mat<-matrix(0,M,d)
for(m in 1:M){
   theta<-jittering_mixture(N=30000, target, observations[1:T_end,], iota= 0.05, cl=1, seed=seed_vec[m])
   bar_mean<-apply(theta$MEAN,2,mean)
   bar_mean2<-apply(theta$MEAN[1:T_end2,],2,mean)
   est<-theta$MEAN[T_end,]
   est2<-theta$MEAN[T_end2,]
   res1[m]<- sqrt( sum(bar_mean-theta_star)^2)
   res2[m]<- max( abs(bar_mean-theta_star))
   res3[m]<- sqrt( sum(est-theta_star)^2)
   res4[m]<- max(abs(est-theta_star))
   res11[m]<- sqrt( sum(bar_mean2-theta_star)^2)
   res22[m]<- max( abs(bar_mean2-theta_star))
   res33[m]<- sqrt( sum(est2-theta_star)^2)
   res44[m]<- max(abs(est2-theta_star))
   est_mat[m,]<-theta$MEAN[T_end,]
   print(m)
   print(res4[m])
}
   
##Save data for Figure 5(a)
#write.table(res4,"Data/Mixture/PSMCO_tBoxplotS_M100_iota005_N30000_T10p5.txt",col.names=FALSE, row.names=FALSE)
#write.table(res44,"Data/Mixture/PSMCO_tBoxplotS_M100_iota005_N30000.txt",col.names=FALSE, row.names=FALSE)


#####################################################################################
## G-PFSO: 5 runs with T=2*10^6 iterations
#####################################################################################
#double the sample size
set.seed(9051985)
cov_x<-diag(1,dx-1)
X <- cbind(rep(1,datalength),rmnorm(datalength, rep(0, dx-1) , cov_x)) 
#sample the response variables
proba<- matrix(X[,a1], ncol=length(a1))%*%matrix(theta_star[1:length(a1)],ncol=1)
proba<- exp(-proba)/(1+exp(-proba))
mu1<- matrix(X[,m1], ncol=length(m1))%*%matrix(theta_star[(length(a1)+1):length(c(a1,m1))],ncol=1)
mu2<- matrix(X[,m2], ncol=length(m2))%*%matrix(theta_star[(length(c(a1,m1))+1):length(c(a1,m1,m2))],ncol=1)
sigma1<- exp(-matrix(X[,s1], ncol=length(s1))%*%matrix(theta_star[(length(c(a1,m1,m2))+1):length(c(a1,m1,m2,s1))],ncol=1))
sigma2<- exp(-matrix(X[,s2], ncol=length(s2))%*%matrix(theta_star[(length(c(a1,m1,m2,s1))+1):length(c(a1,m1,m2,s1,s2))],ncol=1))
ind<-rbinom(datalength, 1,proba)
Y<-rep(0,datalength)
Y[ind==1]<- rnorm(sum(ind==1), mu1[ind==1], sigma1[ind==1])
Y[ind==0]<- rnorm(sum(ind==0), mu2[ind==0], sigma2[ind==0])
observations2 <- cbind(Y, X)
rm(proba,mu1,mu2,sigma1,sigma2,ind,X,Y)
observations<-rbind(observations,observations2)
rm(observations2)

## Reduce the value of C_Sigma to 0.5, to be able
## to visualize well the convergence rate
#target$parameters$learning_rate[1]<-0.5

T_end<-nrow(observations)
M1<-5
n_points<-10^5
points<-seq(1,T_end, length.out=n_points)
res1<-matrix(0,M,n_points)
res2<-matrix(0,M,n_points)
use_seed<-seed_vec[1:M1]
for(m in 1:M1){
   theta<-GPFSO_mixture(N=5000, target, observations[1:T_end,], cl=1, seed=use_seed[m])
   bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
   bar_mean<-bar_mean[points,]
   est<-theta$MEAN[T_end,]
   error <-apply((t(bar_mean)-theta_star)^2,2,sum)^{1/2}
   res1[m,]<-error
   error <-apply(abs(t(bar_mean)-theta_star),2,max) 
   res2[m,]<-error
   rm(theta,bar_mean)
   print(m)
   print(res2[m,n_points])
}
tr<-(1:T_end)[points]
errorE<-apply(res1[1:M1,],2,mean)
errorS<-apply(res2[1:M1,],2,mean)

##Save data for Figure 5(c)
#write.table(errorE,"Data/Mixture/errorE_M5.txt",col.names=FALSE, row.names=FALSE)
#write.table(errorS,"Data/Mixture/errorS_M5.txt",col.names=FALSE, row.names=FALSE)
#write.table(tr,"Data/Mixture/tr.txt",col.names=FALSE, row.names=FALSE)

error<-errorE
t_in<-1#10^4
beta<-0.5
work<-(tr)^{-beta}
c<-error[length(error)]/work[length(error)]
work<-c*work
plot(log(tr[t_in:length(tr)]),log(error[t_in:length(tr)]), type='l')
lines(log(tr[t_in:length(tr)]), log(work[t_in:length(tr)]), type='l', col='red')


 
















