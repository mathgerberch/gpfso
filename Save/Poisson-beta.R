

rm(list=ls())
library(mnormt)
library(parallel)
dyn.load('lib/sqmc.so')	
source('src/sqmc.R')
source('algo_PB.R')
source('algo.R')
source('define_models.R')
numCores <- detectCores() 
#####################################################################################
### Choose the example
#####################################################################################
datalength<-10^5
#####################################################################################
### Generate observations
#####################################################################################
#/setting a seed for the RNG
set.seed(8424)
#theta_star<-c(5,2,8)
theta_star<-c(500,2,8)
d<-length(theta_star)
#sigma_y<-5/100
sigma_y<-5
observations<-rep(0,datalength)
for(i in 1:datalength){
   s<-rbeta(1, theta_star[2],theta_star[3])
   x<- rpois(1,s*theta_star[1])
   observations[i]<-rnorm(1, x, sd=sigma_y)
}
observations<-as.matrix(observations)

#####################################################################################
##Choose (h_t) where h_0=0 and h_t=c_h t^{-alpha} when h_t>0
#####################################################################################
alpha<-0.5
h_seq<- (1:datalength)^(-alpha)
#####################################################################################
#Generate (t_p) 
#####################################################################################
scale1<- 10
diff_min<-  10
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
#Every steps 
pert_times1<-1:datalength
stud_YN1<-stud_YN
gap1<-1
#4.Every GAP setps after START observations
gap<-5
start<-5
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
c<-sqrt(1)
Sigma_t<-diag(c,d)
#####################################################################################
## Build the target
#####################################################################################

target_parameters <- list(sd=sigma_y, alpha_val=c(c,alpha), M=100, 
                         prior=c(1/theta_star),
                          h2_t=h_seq^2,T_p=pert_times, student=stud2_YN, Sigma=Sigma_t,
                          nu=50)
                          
target_parameters$scale_M<-target_parameters$M


target <- list(dimension = d,
               rprior = rprior_PB,
               dprior = dprior_PB,
               loglikelihood =loglikelihood_PB, 
               rmutation=rmutation_pi,
               dmutation=dmutation_pi,
               parameters = target_parameters)
               
target2<-target
target2$parameters$scale_M<-100
target2$parameters$T_p<-pert_times1
target2$parameters$student<-stud_YN1
#####################################################################################
## Estimation  
# Very good: gap=10, J=100, M=10, N=10000 until T1=1000, then M=40 (running time is 1.5h
#for 2000 observations
#####################################################################################

#Set tuning and target parameters

tuning_parameters <- list(N=10000, K=10000, ESS_bound=0.7, resampling=SSP_Resampling, J=100) 

#set.seed(30485)
T_end<- 100
start_time <- Sys.time()
theta<-fast_PB(tuning_parameters, target2 , as.matrix(observations[1:T_end]),  seed=sample(1:10^5,1))
#theta<-online_PB(tuning_parameters, target,target2, T1=10^3,  as.matrix(observations[1:T_end]))#, PI_BAR=TRUE)
#theta<-online(tuning_parameters, target, as.matrix(observations[1:T_end]))
end_time <- Sys.time()
##running time is
print(end_time-start_time)



bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
error <-apply(abs(t(bar_mean)-theta_star),2,max)

error <-apply(abs(t(theta$MEAN)-theta_star),2,max)


points<-1:T_end
beta<-0.5
work<-points^{-beta}
cc<-error[length(error)]/work[length(error)]
work<-cc*work

##plot
plot(log(points), log(error), type='l')
lines(log(points), log(work), type='l', col='red')


plot(theta$ACC, type='l')
ESS<-1/sum(theta$TILDE[,d+1]^2)
print(ESS)


set.seed(30485)
T_end<- 10^4
start_time <- Sys.time()
theta<-online(tuning_parameters, target, as.matrix(observations[1:T_end]))
end_time <- Sys.time()
##running time is
print(end_time-start_time)


set.seed(30485)
T_end<-datalength
M<-10
n_points<- 10^5
points<-seq(1,T_end,length.out=n_points)
res1<-matrix(0,M,n_points)
res2<-matrix(0,M,n_points)
tres1<-matrix(0,M,n_points)
tres2<-matrix(0,M,n_points)
var1<-matrix(0,M,n_points)
var2<-matrix(0,M,n_points)
seed_mat<-sample(1:10^5,M)
for(m in 4:4){
   set.seed(seed_mat[m])
   theta<-online_PB(tuning_parameters, target,target2, T1=10^3, as.matrix(observations[1:T_end]))
   bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
   error <-apply((t(bar_mean)-theta_star)^2,2,sum)^{1/2}
   res1[m,]<-error[points]
   error <-apply(abs(t(bar_mean)-theta_star),2,max) 
   res2[m,]<-error[points]
   error <-apply((t(theta$MEAN)-theta_star)^2,2,sum)^{1/2}
   tres1[m,]<-error[points]
   error <-apply(abs(t(theta$MEAN)-theta_star),2,max) 
   tres2[m,]<-error[points]
   var<-theta$MEAN2-(theta$MEAN)^2
   var1[m,]<-apply(var,1,mean)
   var2[m,]<-apply(var,1,max)
   print(m)
}







tr<-(1:T_end)[points]
write.table(theta_star,"Data/PB/theta_star_d3.txt",col.names=FALSE, row.names=FALSE)
write.table(res1,"Data/PB/errorE_d3_M10.txt",col.names=FALSE, row.names=FALSE)
write.table(res2,"Data/PB/errorS_d3_M10.txt",col.names=FALSE, row.names=FALSE)
write.table(var1,"Data/PB/varianceA_d3_M10.txt",col.names=FALSE, row.names=FALSE)
write.table(var2,"Data/PB/varianceS_d3_M10.txt",col.names=FALSE, row.names=FALSE)
write.table(tr,"Data/PB/tr_d3.txt",col.names=FALSE, row.names=FALSE)


















