

rm(list=ls())
library(mnormt)
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
#theta_star <- c(3, rnorm(d-1,mean=0)) ##Version simple (15%) d=20
#theta_star <- c(3, rnorm(d-1,mean=0)) ##Version simple d=5
theta_star <- c(-2, rnorm(d-1,mean=0)) ##Difficult: 75%, d=20, quantile=0.99 (Gap=1, start=15000)
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
gap<-1#15000
start<-15000
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

target_parameters <- list(quantile=0.99, alpha_val=c(c,alpha),  prior_mean = theta_star+10,
                          prior_covariance =  diag(2, d, d),
                          h2_t=h_seq^2,T_p=pert_times, student=stud2_YN, Sigma=Sigma_t,
                          nu=50)

#target_parameters <- list(quantile=0.5, alpha_val=c(c,alpha),  prior_mean = theta_star,
#                          prior_covariance =  diag(2, d, d),
#                          h2_t=h_seq^2,T_p=pert_times, student=stud2_YN, Sigma=Sigma_t,
#                          nu=50)


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
theta_star[1]<- -2+ sigma_y*qnorm(target$parameters$quantile)
#Set tuning and target parameters
tuning_parameters <- list(N=10000, K=10000, J=10,  ESS_bound=0.7, resampling=SSP_Resampling) 

#tuning_parameters <- list(N=20000, K=10000, J=20,  ESS_bound=0.7, resampling=SSP_Resampling) 


target2<-target
target2$parameters$T_p<-datalength+1
##Estimate tilde{pi} and bar{pi}####################
T_end<-50000#10^6 #datalength
 set.seed(30485)
start_time <- Sys.time()
theta<- online(tuning_parameters, target   , observations[1:T_end,], cl=numCores-2)
end_time <- Sys.time()
##running time is
print(end_time-start_time)



#######################################
tuning_parameters <- list(N=10000, K=10000, J=20,  ESS_bound=0.7, resampling=SSP_Resampling) 
T_end<-  datalength

T_end1<-start
 set.seed(30485)
start_time <- Sys.time()
theta1<- online(tuning_parameters, target   , observations[1:T_end1,], cl=numCores-2)
end_time <- Sys.time()
start_time <- Sys.time()
theta<- Online_censored2(T_end1+1,log(theta1$TILDE[,d+1]),theta1$TILDE[,(1:d)], tuning_parameters, target, observations[1:T_end,], cl=numCores-3,seed=sample(1:10^5,1))
end_time <- Sys.time()
##running time is
print(end_time-start_time)

theta$MEAN<-rbind(theta1$MEAN,theta$MEAN)

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


#plot(bar_mean[,k],type='l')
#abline(h=theta_star[k])


error <-apply((t(bar_mean)-theta_star)^2,2,sum)^{1/2}
#error <-apply(abs(t(bar_mean)-theta_star),2,max) 


#error2 <-apply((t(theta2$MEAN)-theta_star)^2,2,sum)^{1/2}

tr<-2:T_end
nn<-length(tr)
MC<-tr^(-1/2)
C<-error[nn]/MC[nn]
work<-C*MC
points<-seq(200000,T_end,length.out=10^5)


plot(log(tr[points]), log(error[points]), type='l')
lines(log(tr[points]), log(work[points]),type='l', col='red')




lines(log(tr[points]), log(error2[points]),type='l', col='red')








####################################

I<-SSP_Resampling(theta2$TILDE[,(d+1)])
theta_post<-theta2$TILDE[I,1:d]

k<-1

df<- data.frame(Method= factor(c(rep("post", tuning_parameters$N),rep("app.post", tuning_parameters$K) )),
		weights	 	= c(theta_post[,k],theta3$BAR[,k])
		)

p1<-ggplot(data=df,  aes(x=weights, linetype=Method)) + geom_density(size=1.5, stat = "density") +theme_bw()+geom_vline(xintercept = theta_star[k])+
            theme(axis.text=element_text(size=25, colour="black"), axis.title=element_text(size=25))+theme(legend.text=element_text(size=25))+
            theme(legend.position = 'none')+xlab("")+ylab("density")+scale_linetype_manual(values=c(1,2,3)) +xlim(c(2.5,3.5))#+ylim(c(0,0.0061))


df<- data.frame(Method= factor(c(rep("app. post", tuning_parameters$K) )),
		weights	 	= c(theta3$BAR[,k])
		)

p1<-ggplot(data=df,  aes(x=weights, linetype=Method)) + geom_density(size=1.5, stat = "density") +theme_bw()+geom_vline(xintercept = theta_star[k])+
            theme(axis.text=element_text(size=25, colour="black"), axis.title=element_text(size=25))+theme(legend.text=element_text(size=25))+
            theme(legend.position = 'none')+xlab("")+ylab("density")+scale_linetype_manual(values=c(1,2,3)) +xlim(c(-7,7))#+ylim(c(0,0.0061))



df<- data.frame(Method= factor(c(rep("app. post", tuning_parameters$K),rep("app. post2", tuning_parameters$K) )),
		weights	 	= c(theta$BAR[,k],theta3$BAR[,k])
		)

p1<-ggplot(data=df,  aes(x=weights, linetype=Method)) + geom_density(size=1.5, stat = "density") +theme_bw()+geom_vline(xintercept = theta_star[k])+
            theme(axis.text=element_text(size=25, colour="black"), axis.title=element_text(size=25))+theme(legend.text=element_text(size=25))+
            theme(legend.position = 'none')+xlab("")+ylab("density")+scale_linetype_manual(values=c(1,2,3)) +xlim(c(0,4))#+ylim(c(0,0.0061))













