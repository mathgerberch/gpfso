rm(list=ls())
library(mnormt)
library(parallel)
dyn.load('lib/gpfso.so')	
source('src/function.R')
source('gpfso.R')
numCores <- detectCores() 
#####################################################################################
### Exchange rate data
#####################################################################################
d<-9
data2<-read.table("Data/g_and_k/exchange_rate.txt") 
x1<-data2$V4
x2<-data2$V5
observations<-cbind(diff(log(x1)),diff(log(x2)))
datalength<-nrow(observations)
observations<-100*observations
rm(data2,x1,x2)
#####################################################################################
## Choose target parameters
#####################################################################################
c_Sigma<-10
alpha<-0.8
#####################################################################################
## Build target
#####################################################################################
rprior_gk<- function(nparticles, target_parameters){
       bound<-5.5
       work<-matrix(0,nparticles,d)
       work[,c(1,5)]<-matrix(rnorm(2*nparticles, target_parameters$prior_mean[1], sd=target_parameters$prior_sd[1]),nrow=nparticles,ncol=2,byrow=T)
       work[,c(2,6)]<-matrix(rexp(2*nparticles, rate=1/target_parameters$prior_mean[2]),nrow=nparticles,ncol=2,byrow=T)
       work[,c(3,7)]<-matrix(runif(2*nparticles, min=-bound, max=bound), nrow=nparticles, ncol=2,byrow=T)
       f<-function(x){return(-0.045-0.01*x^2)}
       work[,c(4,8)]<- f(work[,c(3,7)])+matrix(rexp(2*nparticles, rate=1/target_parameters$prior_mean[3]),nrow=nparticles,ncol=2,byrow=T)
       work[,9]<- runif(nparticles,min=-1,max=1)  
       return(work)
}
n_it2<- 10^6
comp_seq<-h_Tp(alpha,n_it2)

target_parameters <- list(prior_mean=c(0,1,1), prior_sd=c(1),
                          h_t=comp_seq$H, student=comp_seq$Tp, Sigma=diag(c_Sigma,d),
                          nu=50)                        
rm(comp_seq)
target <- list(dimension = d,
               rprior = rprior_gk,
               loglikelihood=fast_gandk,
               parameters = target_parameters)
#####################################################################################
## Others
#####################################################################################    
##Objective function for quasi-Newton method            
lb<- -Inf
l<-c(lb,0.001,-5.5,-0.5,lb,0.001,-5.5,-0.5,-0.999)
u<-c(-lb,-lb, 5.5,-lb,-lb,-lb, 5.5,-lb,0.999)
lik_g_and_k<-function(x){
     if( min( x[c(4,8)]+0.045+0.01*x[c(3,7)]^2)<0){
       return(99999)
     }else{
        return(-fast_gandk_multi_obs(x,observations, cl=8))
     }
}
##Seed and set the number M of replications
set.seed(30485)
M<-100
seed_vec<-sample(1:10^5,M)
#########################################################################################
# G-PFSO: 100 runs with T=50000 iterations and with (alpha,c_Sigma,N)=(0.8,10,500)
#########################################################################################
n_it<-  50000
n_it1<- 10000
M<-100
res1<-rep(0,M)
res2<-rep(0,M)
est_mat1<-matrix(0,M,d)
est_mat11<-matrix(0,M,d)
est_mat2<-matrix(0,M,d)
est_mat22<-matrix(0,M,d)
for(m in 1:M){
   set.seed(seed_vec[m])
   use_obs<-observations[sample(1:datalength, n_it2, replace=T),]
   theta<- G_PFSO(N=500, target, use_obs[1:n_it,], cl=10)
   bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
   est_mat1[m,]<-bar_mean[nrow(bar_mean),]
   est_mat11[m,]<-bar_mean[n_it1,]
   res1[m]<-fast_gandk_multi_obs(est_mat1[m,], matrix(observations,ncol=2),cl=8)
   est_mat2[m,]<-theta$MEAN[nrow(theta$MEAN),]
   est_mat22[m,]<-theta$MEAN[n_it1,]
   res2[m]<- fast_gandk_multi_obs(est_mat2[m,], matrix(observations,ncol=2),cl=8)
   print(m) 
}

#Computation of the MLE
####################################
if(max(res1)>max(res1)){
   max_ind<-(1:M)[res1==max(res1)]
   mle_est<-est_mat1[max_ind,]
   lik_est<-max(res1)
}else{
   max_ind<-(1:M)[res2==max(res2)]
   mle_est<-est_mat2[max_ind,]
   lik_est<-max(res2)
}
print(lik_est)
est_qN<-optim(mle_est,lik_g_and_k, method="L-BFGS-B", lower=l, upper=u)
lik_mle<-est_qN$val
print(lik_mle)
theta_star<-est_qN$par
##Save MLE, used for Figure 5(a)
write.table(theta_star,"Data/g_and_k/mle.txt",col.names=FALSE, row.names=FALSE)

#Compute estimation errors
####################################
res3<-apply(abs( t(est_mat1)-theta_star),2,max)
res33<-apply(abs( t(est_mat11)-theta_star),2,max)
res4<-apply(abs( t(est_mat2)-theta_star),2,max)
res44<-apply(abs( t(est_mat22)-theta_star),2,max)

##Save data for Figure 5(b)
write.table(res33, "Data/g_and_k/errorS_Boxplot_nit10000.txt",col.names=FALSE, row.names=FALSE)
write.table(res44, "Data/g_and_k/errorS_tBoxplot_nit10000.txt",col.names=FALSE, row.names=FALSE)
write.table(res3, "Data/g_and_k/errorS_Boxplot.txt",col.names=FALSE, row.names=FALSE)
write.table(res4, "Data/g_and_k/errorS_tBoxplot.txt",col.names=FALSE, row.names=FALSE)


#convergence to a local maximum
####################################
ind<-(1:M)[res4>4]

#The different local optima to which
#the algorithm has converged
#(5 elements of Theta_{loc,star} are:
l_max1<-est_mat1[ind[1],]
l_max2<-est_mat1[ind[2],]
l_max3<-est_mat1[ind[4],]
l_max4<-est_mat1[ind[5],]
##Save these local optima, the first one being used for Figure 5(a)
write.table(cbind(l_max1,l_max2,l_max3,l_max4), "Data/g_and_k/local_maxima.txt",col.names=FALSE, row.names=FALSE)
   


#####################################################################################
# BFGS randomly initialized. We observe a large range of values for the resulting
# estimation error, suggesting the existence of a large number of stationnary points,
# and potentially the existence of more local maxima than those indentified in the paper.
# (Warning: numerical error for some values of m)
######################################################################################
#theta_star<-read.table("Data/g_and_k/mle.txt")$V1
#M<-100
#res1<-rep(0,M)
#res2<-rep(0,M)
#est_mat<-matrix(0,M,d)
#for(m in 1:M){
#   val<-target$rprior(1, target$parameters)
#   est<-optim(val,lik_g_and_k, method="L-BFGS-B", lower=l, upper=u))
#   est_mat[m,]<-est$par
#   res1[m]<- est$val
#   res2[m]<-max(abs(est$par-theta_star))
#   print(m)
#   print(res1[m])
#   print(res2[m])
#}

#########################################################################################
# G-PFSO: 100 runs with T=50000 iterations and with (alpha,c_Sigma,N)=(0.5,10,500)
#########################################################################################
theta_star<-read.table("Data/g_and_k/mle.txt")$V1
n_it2<-5*10^6
alpha<-0.5
comp_seq<-h_Tp(alpha,n_it2)
target$parameters$h_t<-comp_seq$H                     
target$parameters$student<-comp_seq$Tp
n_it<-  50000
n_it1<- 10000
M<-100
res1<-rep(0,M)
res2<-rep(0,M)
est_mat1<-matrix(0,M,d)
est_mat11<-matrix(0,M,d)
est_mat2<-matrix(0,M,d)
est_mat22<-matrix(0,M,d)
for(m in 1:M){
   set.seed(seed_vec[m])
   use_obs<-observations[sample(1:datalength, n_it2, replace=T),]
   theta<- G_PFSO(N=500, target, use_obs[1:n_it,], cl=10)
   bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
   est_mat1[m,]<-bar_mean[nrow(bar_mean),]
   est_mat11[m,]<-bar_mean[n_it1,]
   res1[m]<- fast_gandk_multi_obs(est_mat1[m,], observations, cl=10)
   est_mat2[m,]<-theta$MEAN[nrow(theta$MEAN),]
   est_mat22[m,]<-theta$MEAN[n_it1,]
   res2[m]<- fast_gandk_multi_obs(est_mat2[m,], observations, cl=10)
   print(m) 
   print(max(abs(est_mat2[m,]-theta_star)))
}
#Compute estimation errors
####################################
res3<-apply(abs( t(est_mat1)-theta_star),2,max)
res33<-apply(abs( t(est_mat11)-theta_star),2,max)
res4<-apply(abs( t(est_mat2)-theta_star),2,max)
res44<-apply(abs( t(est_mat22)-theta_star),2,max)

##save data for figure 5(b)
write.table(res3, "Data/g_and_k/errorS_Boxplot_alpha05.txt",col.names=FALSE, row.names=FALSE)


#####################################################################################
# G-PFSO: 100 runs with T=50000 iterations, (alpha,c_Sigma,N)=(0.5,1,500)
######################################################################################
theta_star<-read.table("Data/g_and_k/mle.txt")$V1
n_it2<-5*10^6
alpha<-0.5
c_Sigma<-1
comp_seq<-h_Tp(alpha,n_it2)
target$parameters$h_t<-comp_seq$H                   
target$parameters$student<-comp_seq$Tp
target$parameters$Sigma<-diag(c_Sigma,d)
rm(comp_seq) 
n_it<-  50000
n_it1<- 10000
M<-100
res1<-rep(0,M)
res2<-rep(0,M)
est_mat1<-matrix(0,M,d)
est_mat11<-matrix(0,M,d)
est_mat2<-matrix(0,M,d)
est_mat22<-matrix(0,M,d)
for(m in 1:M){
   set.seed(seed_vec[m])
   use_obs<-observations[sample(1:datalength, n_it2, replace=T),]
   theta<- G_PFSO(N=500, target, use_obs[1:n_it,], cl=8)
   bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
   est_mat1[m,]<-bar_mean[nrow(bar_mean),]
   est_mat11[m,]<-bar_mean[n_it1,]
   res1[m]<- fast_gandk_multi_obs(est_mat1[m,], observations, cl=8)
   est_mat2[m,]<-theta$MEAN[nrow(theta$MEAN),]
   est_mat22[m,]<-theta$MEAN[n_it1,]
   res2[m]<- fast_gandk_multi_obs(est_mat2[m,], observations, cl=8)
   print(m) 
   print(max(abs(est_mat2[m,]-theta_star)))
}


#Compute estimation errors
####################################
res3<-apply(abs( t(est_mat1)-theta_star),2,max)
res33<-apply(abs( t(est_mat11)-theta_star),2,max)
res4<-apply(abs( t(est_mat2)-theta_star),2,max)
res44<-apply(abs( t(est_mat22)-theta_star),2,max)

##save data for figure 5(b)
write.table(res3, "Data/g_and_k/errorS_Boxplot_alpha05_c1_N500.txt",col.names=FALSE, row.names=FALSE)


#####################################################################################
# G-PFSO: 5 runs with 10^6 iterations (alpha,c_Sigma,N)=(0.5,1,500)
######################################################################################
theta_star<-read.table("Data/g_and_k/mle.txt")$V1
n_it2<-5*10^6
n_it<- 2*10^6
seed_use<-seed_vec[1:5]
M<-length(seed_use)
res<-matrix(0,M,n_it)
for(m in 1:M){
   set.seed(seed_use[m])
   use_obs<-observations[sample(1:datalength, n_it2, replace=T),]
   theta<- G_PFSO(N=500, target, use_obs[1:n_it,], cl=10)
   bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
   res[m,]<- sqrt(apply((t(bar_mean)-theta_star)^2,2,sum))
   print(m) 
}

#plot error rate
n_points<-10^5
points<-seq(10^4,n_it, length.out=n_points)
error<-apply(res,2,mean)

##save data for figure 5(c)
write.table(error[points], "Data/g_and_k/errorE_c1_500.txt",col.names=FALSE, row.names=FALSE)
write.table((1:n_it)[points], "Data/g_and_k/tr.txt",col.names=FALSE, row.names=FALSE)



