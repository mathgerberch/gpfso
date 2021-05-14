

rm(list=ls())

library(mnormt)
library(metaSEM)
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
datalength<-10^5
n_items<-5
mu_z<-1
#####################################################################################
### Generate observations
#####################################################################################
#/setting a seed for the RNG
set.seed(9585)
d<-n_items*3
#theta_star<-matrix(rnorm(d,0,1),n_items,3)

theta_star<-matrix(rnorm(d,0,2),n_items,3)

observations<-matrix(0,datalength,n_items)
work<-(1/(1+exp(theta_star[,3])))
for(i in 1:datalength){ 
   prob.val<-work+(1-work)/(1+exp(theta_star[,1]*(theta_star[,2]-rnorm(1, mean=1))))
   observations[i,]<-rbinom(n_items, size = 1, prob.val)
   
}
theta_star<-c(t(theta_star))
 

apply(observations,2,mean)

#####################################################################################
##Choose (h_t) where h_0=0 and h_t=c_h t^{-alpha} when h_t>0
#####################################################################################
alpha<-0.3
h_seq<-(1:datalength)^(-alpha)
##number of times each h_t is repeated
#num_rep<-10
#h_seq2<-{}
#i<-1
#while(length(h_seq2)<=(datalength-num_rep)){
# h_seq2<-c(h_seq2,rep(h_seq[i],num_rep))
# i<-i+1
#}
#h_seq<-h_seq2
#####################################################################################
#Generate (t_p) 
#####################################################################################
scale1<-  10
diff_min<-10
f_m<-function(x){return (max(0.5,log(x)))} #need f_m(x)<=x^beta, beta<alpha
#f_m<-function(x){return (max(0.5,x^{alpha*0.8}))}
stud_YN<-rep(0,datalength)
t_0<- 10
t_1<-1
while(t_1 <= datalength){
  t_1<-t_0+max(diff_min,scale1*floor(f_m(t_0)*log(1/h_seq[t_0])))
  if(t_1 <= datalength) stud_YN[t_1]<-1
  t_0<-t_1
}

#alpha_in<-0.3
#c_h<-1
#t_in<-5000
#h_seq<-c(c_h*(1:t_in)^(-alpha_in), h_seq[(t_in+1):datalength])
#####################################################################################
#Choose when perturbations are applied (4 options)
#####################################################################################
#1.Every steps after START observations
start<-1
pert_times<-start:datalength
stud2_YN<-stud_YN
#2.Every steps before END observations
end<-1
pert_times<-1:end
stud2_YN<-stud_YN
#3. only when student noise after START obserations
start<-1000
pert_times<-c(1:datalength)[stud_YN==1]
pert_times<-pert_times[pert_times>start]
stud2_YN<-stud_YN
print(max(diff(pert_times))) #compute largest gap

#4.Every GAP setps after START observations
gap<- 5 
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
#5.Every setps before START and every GAP steps after START observations
gap<-10
start<-100
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
Sigma_t<-0.1*diag(rep(c(1,1,1), n_items),d)


#####################################################################################
## Build the target
#####################################################################################

target_parameters <- list(scale_M=30, prior_mean = rep(0,d), prior_covariance =  diag(1, d, d),
                          h2_t=h_seq^2,T_p=pert_times, student=stud2_YN, Sigma=Sigma_t,
                          nu=50)

target <- list(dimension = d,
               rprior = rprior_quantile,
               dprior = dprior_quantile,
               loglikelihood =loglikelihood_RI, 
               rmutation=rmutation_pi,
               dmutation=dmutation_pi,
               parameters = target_parameters)
#####################################################################################
## Estimation  
#####################################################################################

#Set tuning and target parameters
tuning_parameters <- list(N=10000, K=10000, ESS_bound=0.9, resampling=SSP_Resampling) 


target$parameters$T_p<-pert_times
target$parameters$h2_t<-h_seq^2
target$parameters$student<-stud2_YN
#target$parameters$T_p<-datalength+1
tuning_parameters$N<-20000 

##Estimate tilde{pi} and bar{pi}####################
T_end<-20000#datalength
 set.seed(30485)
start_time <- Sys.time()
theta<- online(tuning_parameters, target, as.matrix(observations[1:T_end,]), PI_BAR=TRUE)
end_time <- Sys.time()
##running time is
print(end_time-start_time)

#target$parameters$Sigma<-0.5*diag(diag(cov(theta$BAR)),d)
 

##Check output####################
#plot ess
plot(theta$ESS/tuning_parameters$N,type='l')

#mean of tilde{pi}
plot(theta$MEAN[,k],type='l')
abline(h=theta_star[k])#  , v=(1:length(theta$MEAN[,k]))[stud2_YN==1])


#mean of bar{pi}
bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))

plot(bar_mean[,k],type='l')
abline(h=theta_star[k])



##plot estimated marginal distributions
df<- data.frame(Method= factor(c(rep("BAR", tuning_parameters$K) )),
		weights	 	= c(theta$BAR[,k])
		)

p1<-ggplot(data=df,  aes(x=weights, linetype=Method)) + geom_density(size=1.5, stat = "density") +theme_bw()+geom_vline(xintercept =theta_star[k])+
            theme(legend.position = 'none')+xlab("")+ylab("")+scale_linetype_manual(values=c(1,2,3))#+xlim(c(mintheta$BAR[,k]))


#####################################################################################
## Estimate tilde{pi}_t and pi_t with same running time 
#####################################################################################
target$parameters$T_p<-pert_times
target$parameters$h2_t<-h_seq^2
target$parameters$student<-stud2_YN
tuning_parameters$N<-10000
target$parameters$prior_covariance<-  diag(c_h^2, d, d)

T_end<-200000#datalength
set.seed(30485)
start_time <- Sys.time()
theta<- online(tuning_parameters, target, observations[1:T_end,], cl=numCores-2)
end_time <- Sys.time()


time_max<-as.numeric(end_time)-as.numeric(start_time)
#tuning_parameters$N<-30000
target2<-target 
target2$parameters$T_p<-datalength+1
set.seed(30485) 
start_time2 <- Sys.time()
theta_post<- online_time(tuning_parameters, time_max, target2, as.matrix(observations),  cl=numCores-2)
end_time2 <- Sys.time()
##running time is
print(end_time2-start_time2)

plot(bar_mean[,k],type='l')
lines(theta_post$MEAN[,k], type='l', col='red')
abline(h=theta_star[k])


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

#####################################################################################
## MMD distance with true posterior
#####################################################################################


yy_mmd<-XY(X=theta_post$TILDE[,1:d], Wx=theta_post$TILDE[,(d+1)], Y=theta_post$TILDE[,1:d], Wy=theta_post$TILDE[,(d+1)],c=2)

xy_mmd<-XY(X=theta_post$TILDE[,1:d], Wx=theta_post$TILDE[,(d+1)], Y=theta$BAR,c=2)
xx_mmd<-XY(X=theta$BAR,   Y=theta$BAR,c=2)

#MMD distance is
sqrt(xx_mmd+yy_mmd-2*xy_mmd)

#####################################################################################
## Plot ESS
#####################################################################################

points2<-seq(100,length(theta_post$ESS),100) 

df<- data.frame(Method= factor(c(rep("2",  length(points2) ) )),
		axis	 	= c(points2/1000),
                res 		= c( theta_post$ESS[points2]/tuning_parameters$N)
		)



p1<-ggplot(data=df,  aes(x=axis, y=res, group=Method, linetype=Method)) + geom_line() +
        scale_x_continuous( breaks = c(50,100,150,200,250))+
	xlab(expression("")) +ylab("") +theme_bw()+
        theme(legend.title=element_blank())+theme(legend.position="bottom")+
	theme(axis.text=element_text(size=25, colour="black"),
        axis.title=element_text(size=25))+theme(legend.text=element_text(size=25))+
	scale_linetype_manual(values=c(1, 3,2))+
        scale_size_manual(values=c(0.1,2))+
	theme(legend.key.size = unit(1.2, "cm"))+ theme(legend.position = 'none')
pdf("ess.pdf")
p1
dev.off()


#####################################################################################
## Convergence rate of the means
#####################################################################################
bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
error <-apply(abs(t(bar_mean)-theta_star),2,max)
#error2 <-apply(abs(t(theta_post$MEAN)-theta_star),2,max)



##rate beta
beta<-0.5
work<-(1:length(error))^{-beta}
c<-error[length(error)]/work[length(error)]
work<-c*work

##plot
#T_end<-datalength
plot(log(2:length(error)), log(error[2:length(error)]), type='l', col="blue")
lines(log(2:length(error)), log(work[2:T_end]), type='l', col='red')


##nicer plot
points<-seq(1000,T_end,1000)
points2<-seq(1000,length(error2),1000) 

df<- data.frame(Method= factor(c(rep("1", length(points)),rep("2", length(points2)),rep("3", length(points)) )),
		axis	 	= c(points,points2,points),
                res 		= c(error[points],error2[points2],work[points])
		)



p1<-ggplot(data=df,  aes(x=axis, y=res, group=Method, linetype=Method)) + geom_line() +geom_hline(yintercept = error2[length(theta_post$ESS)])+
        geom_vline(xintercept =  length(theta_post$ESS))+# expand_limits(y=0)+
	scale_x_log10(limits = c(1000,T_end),breaks = c(10^{3}  ,10^{4},10^5,10^6  ), labels = trans_format("log10", math_format(10^.x)))+
	scale_y_log10(limits = c(  10^{-6}, 2), breaks = c(10^{-6},10^{-4},10^{-2},10   ),labels = trans_format("log10", math_format(10^.x)))+
	xlab(expression("")) +ylab("") +theme_bw()+
        theme(legend.title=element_blank())+theme(legend.position="bottom")+
	theme(axis.text=element_text(size=25, colour="black"),
        axis.title=element_text(size=25))+theme(legend.text=element_text(size=25))+
	scale_linetype_manual(values=c(1, 3,2))+
        scale_size_manual(values=c(0.1,2))+
	theme(legend.key.size = unit(1.2, "cm"))+ theme(legend.position = 'none')+
        annotate("text", x = 10^{3.1}, y = 10^{-1.4},  label = expression(t^{-0.5}), size=10)


pdf("quantile_d30_alpha03_rate_bar_q05_Delta100_N10000_N30000_comp1.pdf")
p1
dev.off()





##Compare ratio mean and variances
bar_mean2<-apply(theta$MEAN2,2,cumsum)/(1:nrow(theta$MEAN))
bar_var<-bar_mean2-bar_mean^2
post_var<-theta_post$MEAN2-theta_post$MEAN^2


post_sd<-sqrt(post_var[length(theta_post$ESS),])
post_bais<-abs(theta_post$MEAN[length(theta_post$ESS),]-theta_star)


#write.table(bar_mean, "Data_figs/bar_mean_d30_alpha03_q05_Delta100_N10000.txt")
#write.table(bar_var, "Data_figs/bar_var_d30_alpha03_q09_Delta100_N10000.txt")
#write.table(theta_post$MEAN, "Data_figs/post_mean_d30_q05_N30000.txt")
#write.table(post_var, "Data_figs/post_var_d30_q05_N30000.txt")
#write.table(theta_post$ESS, "Data_figs/post_ess_d30_q05_N30000.txt")






bar_sd<-sqrt(bar_var[length(theta$ESS),])
bar_bais<-abs(bar_mean[length(theta$ESS),]-theta_star)

##Boxplot
df<- data.frame(Method= factor(c(rep("1", d),rep("2", d) )),
		axis	 	= c(bar_bais/bar_sd,post_bais/post_sd) 
		)

p1 <- ggplot(df, aes(x = Method, y = axis)) +ylab("") +theme_bw()+ scale_x_discrete(name = " ")+
        geom_boxplot()+scale_y_log10(  breaks = c(10^{-3},10^{-2},10^{-1},10^{0}  , 10,10^2 ),labels = trans_format("log10", math_format(10^.x)))+
        theme(legend.title=element_blank())+theme(legend.position="bottom")+
	theme(axis.text=element_text(size=25, colour="black"),
        axis.title=element_text(size=25))+theme(legend.text=element_text(size=25),axis.text.x=element_blank()) 
pdf('boxplot.pdf')
p1
dev.off()



#####################################################################################
## Estimated marginal distributions
#####################################################################################

##posterior distribution
df <- data.frame(Method= factor(c(rep("1", nrow(theta_post$TILDE)))),
		values	 	= c(theta_post$TILDE[,k]),
                weight_val	= c(theta_post$TILDE[,(d+1)])
		)
p1<-ggplot(data=df,  aes(x=values,  weight=weight_val, linetype=Method)) + geom_density(size=1.5, stat = "density") +theme_bw()+geom_vline(xintercept =theta_star[k])+
            theme(legend.position = 'none')+xlab("")+ylab("")+scale_linetype_manual(values=c(1,2,3))+xlim(c(theta_star[k]-0.1,theta_star[k]+0.1))


##bar{pi}
df<- data.frame(Method= factor(c(rep("BAR", tuning_parameters$K) )),
		weights	 	= c(theta$BAR[,k])
		)

p1<-ggplot(data=df,  aes(x=weights, linetype=Method)) + geom_density(size=1.5, stat = "density") +theme_bw()+geom_vline(xintercept =theta_star[k])+
            theme(legend.position = 'none')+xlab("")+ylab("")+scale_linetype_manual(values=c(1,2,3))#+xlim(c(mintheta$BAR[,k]))


pdf('quantile_d5_5.pdf')
p1
dev.off()
 

f(pert_times))) #compute largest gap

#4.Every GAP setps after START observations
gap<-100   
start<-100 
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
gap<-50
start<-500
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

target_parameters <- list(quantile=0.9,   prior_mean = rep(0,d),
                          prior_covariance =  diag(4, d, d),
                          h2_t=h_seq^2,T_p=pert_times, student=stud2_YN, Sigma=Sigma_t,
                          nu=50)

target <- list(dimension = d,
               rprior = rprior_quantile,
               dprior = dprior_quantile,
               loglikelihood =loglikelihood_quantile, 
               rmutation=rmutation_quantile,
               dmutation=dmutation_quantile,
               parameters = target_parameters)
#####################################################################################
## Estimation  
#####################################################################################

#Set tuning and target parameters
tuning_parameters <- list(N=10000, K=10000, J=1, ESS_bound=0.7, resampling=SSP_Resampling) 


target$parameters$T_p<-pert_times
target$parameters$h2_t<-h_seq^2
target$parameters$student<-stud2_YN
target$parameters$quantile<-0.9

#tuning_parameters$N<-2000
#tuning_parameters$K<-10000

##Estimate tilde{pi} and bar{pi}####################
T_end<-100000# datalength
set.seed(30485)
start_time <- Sys.time()
theta<- online(tuning_parameters, target, observations[1:T_end,], theta_star, cl=numCores-2)
end_time <- Sys.time()
##running time is
print(end_time-start_time)


#mean of tilde{pi}
plot(theta$MEAN[,k],type='l')



work<-matrix(0,nrow(theta$MEAN), d)
for(i in 1:d) work[,i]<-theta$MEAN[,i]/scale_values[i]

work[,1]<-theta$MEAN[,1]
for(i in 2:d) work[,1]<-work[,1]+work[,i]*mean_values[i]




#mean of bar{pi}
bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))

plot(bar_mean[,k],type='l')
 











