

rm(list=ls())

library(mnormt)
library(gtools)
library(ggplot2)
library(scales)
library(parallel)
dyn.load('lib/sqmc.so')	
source('src/sqmc.R')
#source('algo_logistic.R')
source('algo.R')
source('algo_time.R')
source('define_models.R')
numCores <- detectCores() 
#####################################################################################
### Choose the example
#####################################################################################
d<-2
datalength<-10^6
##############################
set.seed(9585)
a_star<--2


f0<-function(x) sin(10*x)
 
 
plot(f0(seq(0,1,0.01)), type='l', ylim=c(-1,1))
X<-cbind( rep(1,datalength),rmnorm(datalength, 1, 1))
tseq<-(1:datalength)/datalength
theta_star<-1+0.2*f0(tseq)
 

proba<-1/(1+exp(-a_star-theta_star*X[,2]))

Y<-rbinom(datalength,1,proba)
observations <- cbind(Y, X)
rm(X,Y,proba)
#####################################################################################
### Covertype data
#####################################################################################
set.seed(9585)
data<-read.csv("Data/covtype.data")
datalength<-nrow(data)
Y<-rep(0,datalength)
Y[data[,55]==1]<-1
X<-cbind(rep(1,datalength),data[,1:53])
d<-ncol(X)
observations<-cbind(Y,X)
rm(X,Y,data)
observations<-observations[sample(1:datalength, datalength,replace=FALSE),]
observations<-apply(observations,2,as.numeric)



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
#2.Every steps before END observations
end<-1000
pert_times<-1:end
stud2_YN<-stud_YN
#3. only when student noise after START obserations
start<-1000
pert_times<-c(1:datalength)[stud_YN==1]
pert_times<-pert_times[pert_times>start]
stud2_YN<-stud_YN
print(max(diff(pert_times))) #compute largest gap

#4.Every GAP setps after START observations
gap<-1000   
start<-1 
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
#5.Every setps before START and every GAP steps after START observations
gap<-10^6
start<-10000
pert_times<-c(1:start, seq(start+1, datalength, gap))
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
Sigma_t<-diag(1,d)

#####################################################################################
## Build the target
#####################################################################################

target_parameters <- list(prior_mean = rep(0,d),
                          prior_covariance =  diag(2, d, d),
                          h2_t=h_seq^2,T_p=pert_times, student=stud2_YN, Sigma=Sigma_t,
                          nu=50)
target_parameters$prior_sd<- sqrt(diag(target_parameters$prior_covariance))
target_parameters$alpha<-c(1,alpha)
 

target <- list(dimension = d,
               rprior = rprior_quantile,
               dprior = dprior_quantile,
               loglikelihood =loglikelihood_logistic, 
               rmutation=rmutation_pi,
               dmutation=dmutation_pi,
               parameters = target_parameters)
#####################################################################################
## Estimation  
#####################################################################################

#Set tuning and target parameters
tuning_parameters <- list(N=2000, K=1000,  ESS_bound=0.7, resampling=SSP_Resampling, J=100) 



##Estimate tilde{pi} and bar{pi}####################
T_end<-10000#  datalength
seed_val<-30485
set.seed(30485)
start_time <- Sys.time()
theta<- online(tuning_parameters, target, observations[1:T_end,], PI_BAR=FALSE, cl=6)
end_time <- Sys.time()
##running time is
print(end_time-start_time) #10.51655 mins

tuning_parameters <- list(N=5000, K=1000,  ESS_bound=0.7, resampling=SSP_Resampling, J=100) 

T_end<-20000
target2<-target 
target2$parameters$T_p<-datalength+1
#seed_val<-30485
#set.seed(30485)
start_time <- Sys.time()
theta2<- online(tuning_parameters, target, observations[1:T_end,], PI_BAR=FALSE, cl=6)
end_time <- Sys.time()
##running time is
print(end_time-start_time) #26.51655 mins






time_max<-as.numeric(end_time)-as.numeric(start_time)

target2<-target 
target2$parameters$T_p<-datalength+1

set.seed(30485) 
start_time2 <- Sys.time()
theta_post<- online_time(tuning_parameters, time_max, target, observations[1:T_end,], cl=numCores-2)
end_time2 <- Sys.time()
##running time is
print(end_time2-start_time2)

est<-compute_loss(theta2, observations[1:T_end,], target)

#write.table(theta$MEAN[,2], 'betaTilde_e2.txt')
#write.table(bar_mean[,2], 'betaBar_e2.txt')
#write.table(theta$MEAN[,1], 'aTilde_e2.txt')
#write.table(bar_mean[,1], 'a_e2.txt')

bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))

start<-100
plot(theta$MEAN[start:T_end,2], type='l')
lines(bar_mean[start:T_end,2], type='l', col='blue')
lines(theta_star, type='l', col='red')




start<-100
plot(theta$MEAN[start:T_end,1], type='l')
lines(bar_mean[start:T_end,1], type='l', col='blue')
abline(h=a_star,   col='red')



target2<-target 
target2$parameters$T_p<-datalength+1
set.seed(30485)
start_time <- Sys.time()
theta2<- online(tuning_parameters, target2, observations[1:T_end,], PI_BAR=FALSE, cl=numCores-2)
end_time <- Sys.time()
##running time is
print(end_time-start_time)


start<-100

plot(theta2$MEAN[start:T_end,2], type='l')
lines(theta_star, type='l', col='red')


write.table(theta2$MEAN[,2], 'betaPost_e2.txt')
write.table(theta2$MEAN[,1], 'aPost_e2.txt')
###################################################




T_end<-  datalength
start_time <- Sys.time()
est<-Online_logistic(tuning_parameters$N, observations[1:T_end,], target_parameters, tuning_parameters$ESS_bound, cl=8, seed=30485)
end_time <- Sys.time()
##running time is 26 minutes
print(end_time-start_time)





loss<-pred$LOSS
loss.b<-pred$LOSS.B

#loss<-pred2$LOSS
#loss.b<-pred2$LOSS.B



p1<-loss[observations[2:T_end,1]==k]
work<-rep(0,length(p1))
work[p1<log(2)]<-1
p1<-work
p1<-cumsum(p1)/(1:length(p1))

p2<-loss.b[observations[2:T_end,1]==k]
work<-rep(0,length(p2))
work[p2<log(2)]<-1
p2<-work
p2<-cumsum(p2)/(1:length(p2))

plot(p2, type='l', col='red')
lines(p1, type='l')



p1<-loss
work<-rep(0,length(p1))
work[p1<log(2)]<-1
p1<-work
p1<-cumsum(p1)/(1:length(p1))

p2<-loss.b
work<-rep(0,length(p2))
work[p2<log(2)]<-1
p2<-work
p2<-cumsum(p2)/(1:length(p2))

plot(p2, type='l', col='red')
lines(p1, type='l')



##Check output####################

#plot ess
plot(theta$ESS/tuning_parameters$N,type='l')

#plot acceptance rate
plot(theta$ACC,type='l')

k<-2
#mean of tilde{pi}
plot(theta$MEAN[,k],type='l')
lines(theta_star,type='l', col='red')
#abline(h=theta_star[k])


#mean of bar{pi}

loss.b<-rep(0,datalength-1)
thetab<-0
for(t in 1:(datalength-1)){
      thetab<-(t-1)*thetab/t+theta$MEAN[t,]/t
      loss.b[t]<- -loglikelihood_logistic(matrix(thetab,1,d), matrix(observations[t+1,], ncol=d+1),target$parameters, cl=1)
  }

#write.table(est$LOSS, "Data/Covertype/loss.txt", col.names=FALSE, row.names=FALSE)
#write.table(est$LOSS.B, "Data/Covertype/lossb.txt", col.names=FALSE, row.names=FALSE)
#write.table(est2$LOSS, "lossSA_low.txt")
#write.table(est2$LOSS.B, "lossSab_low.txt")



loss<-read.table("loss.txt")$x
loss.b<-read.table("lossb.txt")$x
loss.sa<-read.table('lossSA.txt')$x
loss.sab<-read.table('lossSAb.txt')$x

k=1

p1<-loss[observations[2:datalength,1]==k]
work<-rep(0,length(p1))
work[p1<log(2)]<-1
p1<-work
p1<-cumsum(p1)/(1:length(p1))

p2<-loss.b[observations[2:datalength,1]==k]
work<-rep(0,length(p2))
work[p2<log(2)]<-1
p2<-work
p2<-cumsum(p2)/(1:length(p2))



p3<-loss.sa[observations[2:datalength,1]==k]
work<-rep(0,length(p3))
work[p3<log(2)]<-1
p3<-work
p3<-cumsum(p3)/(1:length(p3))

p4<-loss.sab[observations[2:datalength,1]==k]
work<-rep(0,length(p4))
work[p4<log(2)]<-1
p4<-work
p4<-cumsum(p4)/(1:length(p4))
rm(work)

start<-500

plot(p2[start:length(p2)], type='l', col='red', ylim=c(0,1))
lines(p4[start:length(p4)], type='l', col='blue')
lines(p1[start:length(p1)], type='l')
lines(p3[start:length(p3)], type='l', col='green')


p1<-loss
work<-rep(0,length(p1))
work[p1<log(2)]<-1
p1<-work
p1<-cumsum(p1)/(1:length(p1))

p2<-loss.b
work<-rep(0,length(p2))
work[p2<log(2)]<-1
p2<-work
p2<-cumsum(p2)/(1:length(p2))



p3<-loss.sa
work<-rep(0,length(p3))
work[p3<log(2)]<-1
p3<-work
p3<-cumsum(p3)/(1:length(p3))

p4<-loss.sab
work<-rep(0,length(p4))
work[p4<log(2)]<-1
p4<-work
p4<-cumsum(p4)/(1:length(p4))
rm(work)






start<-10

plot(p2[start:length(p2)], type='l', col='red', ylim=c(0,1))
lines(p4[start:length(p4)], type='l', col='blue')
lines(p1[start:length(p1)], type='l')
lines(p3[start:length(p3)], type='l', col='green')



##SVB##############################



start_time <- Sys.time()
est<-SVB_logistic(N=2000, observations, sigma0=sqrt(2),alpha=0.5, c=0.001, cl=numCores-2, seed=30485)
end_time <- Sys.time()
print(end_time-start_time)


start_time <- Sys.time()
est<-NGVI_logistic(N=1000, observations, sigma0=sqrt(0.2),  alpha=0.5, c=0.001, eta=1, cl=numCores-2,  seed=30485)
end_time <- Sys.time()
print(end_time-start_time)



start_time <- Sys.time()
est2<-SA_logistic(observations,theta_in=rep(0,d), 0.5,0.001)
end_time <- Sys.time()
##running time is  
print(end_time-start_time)




p1<-est$LOSS.B
work<-rep(0,length(p1))
work[p1<log(2)]<-1
p1<-work
p1<-1-cumsum(p1)/(1:length(p1))




p2<-est2$LOSS.B
work<-rep(0,length(p2))
work[p2<log(2)]<-1
p2<-work
p2<-1-cumsum(p2)/(1:length(p2))
rm(work)


p11<-est$LOSS
work<-rep(0,length(p11))
work[p11<log(2)]<-1
p11<-work
p11<-1-cumsum(p11)/(1:length(p11))




p22<-est2$LOSS
work<-rep(0,length(p22))
work[p22<log(2)]<-1
p22<-work
p22<-1-cumsum(p22)/(1:length(p22))
rm(work)





start<-1
plot(p1[start:length(p1)], type='l', col='red', ylim=c(0.2,1))
#lines(p11[start:length(p11)], type='l', col='orange')
lines(p2[start:length(p2)], type='l', col='blue')
#lines(p22[start:length(p22)], type='l', col='green')

###



p3<-est$LOSS.B[observations[2:datalength,1]==k]
work<-rep(0,length(p3))
work[p3<log(2)]<-1
p3<-work
p3<-cumsum(p3)/(1:length(p3))




p4<-est2$LOSS.B[observations[2:datalength,1]==k] 
work<-rep(0,length(p4))
work[p4<log(2)]<-1
p4<-work
p4<-cumsum(p4)/(1:length(p4))
rm(work)


p33<-est$LOSS[observations[2:datalength,1]==k]
work<-rep(0,length(p33))
work[p33<log(2)]<-1
p33<-work
p33<-cumsum(p33)/(1:length(p33))




p44<-est2$LOSS[observations[2:datalength,1]==k] 
work<-rep(0,length(p44))
work[p44<log(2)]<-1
p44<-work
p44<-cumsum(p44)/(1:length(p44))
rm(work)



start<-1

plot(p3[start:length(p3)], type='l', col='red', ylim=c(0,1))
lines(p33[start:length(p33)], type='l', col='orange')
lines(p4[start:length(p4)], type='l', col='blue')
lines(p44[start:length(p44)], type='l', col='green')


 
##SA##############################



start_time <- Sys.time()
est<-SA_logistic(observations,theta_in=rep(0,d), 0.5,0.001)
end_time <- Sys.time()
##running time is  
print(end_time-start_time)





p3<-est$LOSS
work<-rep(0,length(p3))
work[p3<log(2)]<-1
p3<-work
p3<-cumsum(p3)/(1:length(p3))




p4<-est$LOSS.B
work<-rep(0,length(p4))
work[p4<log(2)]<-1
p4<-work
p4<-cumsum(p4)/(1:length(p4))
rm(work)


start<-500

plot(p2[start:length(p2)], type='l', col='red', ylim=c(0,1))
lines(p4[start:length(p4)], type='l', col='blue')
lines(p1[start:length(p1)], type='l')
lines(p3[start:length(p3)], type='l', col='green')






plot(p4, type='l', col='red')
lines(p3, type='l')



#################################
#mean of tilde{pi}
plot(est$THETA[,2],type='l')
lines(theta_star,type='l', col='red')

plot(est$THETA[,1],type='l')
abline(h=a_star)

plot(theta$MEAN[,2],type='l')
lines(theta_star,type='l', col='red')

plot(theta$MEAN[,1],type='l')
abline(h=a_star)

e1<-(est$THETA[,2]-theta_star)^2+(est$THETA[,1]-a_star)^2
e2<-(theta$MEAN[,2]-theta_star)^2+(theta$MEAN[,1]-a_star)^2
plot(e1[(datalength-10000):datalength],type='l', ylim=c(0,0.1))
lines(e2[(datalength-10000):datalength],type='l', col='red')


#mean of tilde{pi}
plot(theta$MEAN[,2],type='l')
lines(theta_star,type='l', col='red')
lines(est$THETA[,2],type='l', col='blue')
#abline(h=theta_star[k])



plot(cumsum(est$LOSS/theta$LOSS)/(1:length(est$LOSS)), type='l')
abline(h=1)


lines(cumsum(theta$LOSS)/(1:length(theta$LOSS)), type='l', col='red')




