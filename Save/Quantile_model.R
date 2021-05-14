

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
d<-5
#####################################################################################
### Generate observations
#####################################################################################
#/setting a seed for the RNG
set.seed(9585)
#/number of observations
datalength <- 10^7
#/covariance matrix of covariates
cov_x<-matrix(rWishart(1,d-1,diag(d-1)),d-1)
cov_x<-solve(cov_x)
#/sample covariates
X <- cbind(rep(1,datalength),rmnorm(datalength, rep(0, d-1) , cov_x))
#/true parameter value
theta_star <- rnorm(d)# 1+5*runif(d) #(d:1) / d
theta_starSave<-theta_star
#/sample response variables
sigma_y<-2
Y <- rnorm(datalength, X %*% theta_star,sd= sigma_y)
#// concatenate Y and X into one matrix
observations <- cbind(Y, X)
rm(X,Y)
#####################################################################################
##Choose (h_t) where h_0=0 and h_t=c_h t^{-alpha} when h_t>0
#####################################################################################
alpha<-0.6
c_h<-1

h_seq<-c(c_h*(1:datalength)^(-alpha))
#####################################################################################
#Generate (t_p) 
#####################################################################################
scale1<-10
diff_min<-10
#f_m<-function(x){return (max(0.5,log(x)))} #need f_m(x)<=x^beta, beta<alpha
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
end<-1
pert_times<-1:end
stud2_YN<-stud_YN
#3. only when student noise after START obserations
start<-100
pert_times<-c(1:datalength)[stud_YN==1]
pert_times<-pert_times[pert_times>start]
stud2_YN<-stud_YN
print(max(diff(pert_times))) #compute largest gap

#4.Every GAP setps after START observations
gap<-1
start<-500
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
gap<-2
start<-2
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

target_parameters <- list(quantile=0.9,   prior_mean = theta_star-10,
                          prior_covariance =  diag(2, d, d),
                          h2_t=h_seq^2,T_p=pert_times, student=stud2_YN, Sigma=Sigma_t,
                          nu=50)

target <- list(dimension = d,
               rprior = rprior_quantile,
               dprior = dprior_quantile,
               loglikelihood =loglikelihood_quantile, 
               rmutation=rmutation_pi,
               dmutation=dmutation_pi,
               parameters = target_parameters)
#####################################################################################
## Estimation  
#####################################################################################

#Set tuning and target parameters
tuning_parameters <- list(N=1000, K=5000,   ESS_bound=0.7, resampling=SSP_Resampling) 


#target$parameters$T_p<-pert_times
#target$parameters$h2_t<-h_seq^2
#target$parameters$student<-stud2_YN
target$parameters$quantile<-0.9

set.seed(30485)
T_end<- 1000# datalength
start_time <- Sys.time()
theta<- online(tuning_parameters, target, observations[1:T_end,], PI_BAR=TRUE, cl=numCores-2) 
end_time <- Sys.time()
##running time is
print(end_time-start_time)

#IBIS
target2<-target 
target2$parameters$T_p<-datalength+1

start_time2 <- Sys.time()
theta_post<- online(tuning_parameters,  target2, observations[1:T_end,],  cl=numCores-2)
end_time2 <- Sys.time()
##running time is
print(end_time2-start_time2)


##############################
# MMD computations
##############################
target$parameters$quantile<-0.9
tuning_parameters$N<-50000
tuning_parameters$K<-50000
set.seed(30485)
observations<-observations[1:100000,]
datalength<-nrow(observations)
target2<-target 
target2$parameters$T_p<-datalength+1

M<-1
t_seq<-c(1000,10000,100000)
alpha_seq<-c(0.3,0.5)
Delta_seq<-c(1,5000,10000)

res_bar<-array(0,dim=c(length(t_seq),length(alpha_seq),length(Delta_seq)))
res_tilde<-array(0,dim=c(length(t_seq),length(alpha_seq),length(Delta_seq)))

for(t in 1:length(t_seq)){
    for(i in 1:M){  
        theta_post<- online(tuning_parameters,  target2, observations[1:t_seq[t],],  cl=numCores-2)
        part_post<-XY_comupute(theta_post$TILDE[,1:d],theta_post$TILDE[,1:d],theta_post$TILDE[,d+1],theta_post$TILDE[,d+1], gamma=1)
        theta_post$ESS<-0
        theta_post$MEAN<-0
        theta_post$MEAN2<-0
        theta_post$ACC<-0
        for(j in 1:length(alpha_seq)){
           for(k in 1:length(Delta_seq)){
               param<-compute_param(alpha_seq[j], Delta_seq[k], datalength, c_h=1)
               target$parameters$h2_t<-param$H^2
               target$parameters$T_p<-param$T_seq
               target$parameters$student<-param$TP
               rm(param)
               theta<- online(tuning_parameters, target, observations[1:t_seq[t],], PI_BAR=TRUE, cl=numCores-2)
               part_bar<-XY_comupute(theta$BAR,theta$BAR,rep(1/tuning_parameters$K,tuning_parameters$K),rep(1/tuning_parameters$K,tuning_parameters$K), gamma=1)
               part_tilde<-XY_comupute(theta$TILDE[,1:d],theta$TILDE[,1:d],theta$TILDE[,d+1],theta$TILDE[,d+1], gamma=1)
               part_barPost<-XY_comupute(theta$BAR,theta_post$TILDE[,1:d],rep(1/tuning_parameters$K,tuning_parameters$K),theta_post$TILDE[,d+1], gamma=1) 
               part_tildePost<-XY_comupute(theta$TILDE[,1:d],theta_post$TILDE[,1:d],theta$TILDE[,d+1],theta_post$TILDE[,d+1], gamma=1)
               theta$ESS<-0
               theta$MEAN<-0
               theta$MEAN2<-0
               theta$ACC<-0
               res_bar[t,j,k]<-res_bar[t,j,k]+sqrt(part_bar+part_post-2*part_barPost)
               res_tilde[t,j,k]<-res_tilde[t,j,k]+sqrt(part_tilde+part_post-2*part_tildePost)
           }

        }

    }
  cat('.')
}



##########################################






#write.table(theta$MEAN, 'Data/Lin_quantile_mean_tilde_q09_alpha01.txt')
#write.table(theta$MEAN2, 'Data/Lin_quantile_mean2_tilde_q09_alpha01.txt')
theta_star[1]<-theta_starSave[1]+sigma_y*qnorm(target_parameters$quantile)
bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
error <-apply(abs(t(bar_mean)-theta_star),2,max)
#write.table(error, 'Data/Lin_quantile_error_bar_q09_alpha01_Delta10000.txt')




##Check output####################
#compute true parameter value
theta_star[1]<-theta_starSave[1]+ sigma_y*qnorm(target$parameters$quantile)
#plot ess
plot(theta$ESS/tuning_parameters$N,type='l')

#mean of tilde{pi}
plot(theta$MEAN[,k],type='l')
abline(h=theta_star[k])


#mean of bar{pi}
bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))

plot(bar_mean[,k],type='l')
abline(h=theta_star[k])


##QMC estimation######################
tuning_parameters <- list(N=2^10, K=10000,  resampling=SSP_Resampling) 
T_end<- datalength
set.seed(30485)
start_time <- Sys.time()
theta<- online_qmc(tuning_parameters, target, observations[1:T_end,], PI_BAR=TRUE, cl=1)
end_time <- Sys.time()
##running time is
print(end_time-start_time)

#compute true parameter value
theta_star<-theta_starSave+ scale_sd*qnorm(target$parameters$quantile)
#plot ess
plot(theta$ESS/tuning_parameters$N,type='l')

#mean of tilde{pi}
plot(theta$MEAN[,k],type='l')
abline(h=theta_star[k])

#########################


#####################################################################################
## Estimate tilde{pi}_t and pi_t with same running time 
#####################################################################################
target$parameters$T_p<-pert_times
target$parameters$h2_t<-h_seq^2
target$parameters$student<-stud2_YN
target$parameters$quantile<-0.9
tuning_parameters$N<-10000
tuning_parameters$K<-10000
#target$parameters$prior_covariance<-  diag(c_h^2, d, d)

T_end<-datalength
set.seed(30485)
start_time <- Sys.time()
theta<- online(tuning_parameters, target, observations[1:T_end,], PI_BAR=TRUE, cl=numCores-2)
end_time <- Sys.time()


time_max<-as.numeric(end_time)-as.numeric(start_time)
#tuning_parameters$N<-30000
target2<-target 
target2$parameters$T_p<-datalength+1
set.seed(30485) 
start_time2 <- Sys.time()
theta_post<- online_time(tuning_parameters, time_max, target2, observations[1:T_end,],  cl=numCores-2)
end_time2 <- Sys.time()
##running time is
print(end_time2-start_time2)


##IBIS
target2<-target 
target2$parameters$T_p<-datalength+1
set.seed(30485) 
start_time2 <- Sys.time()
theta_post<- online(tuning_parameters,  target2, observations[1:T_end,],  cl=numCores-2)
end_time2 <- Sys.time()
##running time is
print(end_time2-start_time2)

######################################
##Check output####################
theta_star[1]<-theta_starSave[1]+sigma_y*qnorm(target$parameters$quantile)


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
theta_star[1]<-theta_starSave[1]+sigma_y*qnorm(target_parameters$quantile)
bar_mean<-apply(theta$MEAN,2,cumsum)/(1:nrow(theta$MEAN))
error <-apply(abs(t(bar_mean)-theta_star),2,max)
bar_mean2<-apply(theta2$MEAN,2,cumsum)/(1:nrow(theta2$MEAN))
error2<-apply(abs(t(bar_mean2)-theta_star),2,max)

error2 <-apply(abs(t(theta_post$MEAN)-theta_star),2,max)



##rate beta
beta<-0.5
work<-(1:length(error))^{-beta}
c<-error[length(error)]/work[length(error)]
work<-c*work

##plot
T_end<-length(error)#datalength
#plot(log(2:length(error2)), log(error2[2:length(error2)]), type='l')
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
            theme(legend.position = 'none')+xlab("")+ylab("")+scale_linetype_manual(values=c(1,2,3))+xlim(c(theta_star[k]-1.5,theta_star[k]+1.5))


pdf('quantile_d5_5.pdf')
p1
dev.off()
 

#####################################################################################
## REAL DATA EXAMPLE
#####################################################################################

#####PREPARE DATA####################################################################
##download data
Data<-read.csv("Data/natl2017.csv", header=TRUE, sep=",")
#select data
weight<-Data$dbwt  #remove 9999 
race<-Data$mbrace  ##1=white
white<-rep(0,length(race))
white[race==1]<-1
black<-rep(0,length(race))
black[race==2]<-1
married<-Data$dmar 
married[married==2]<-0
age<-Data$mager #in year
education<-Data$meduc  
ed1<-rep(0,length(education))
ed1[education<3]<-1 #<=12
ed2<-rep(0,length(education))
ed2[education==3]<-1  
ed2[education==4]<-1 
ed3<-rep(0,length(education))
ed3[education>4]<-1  
ed3[education==9]<--1  #to remove  
weight_gain<-Data$wtgain #remove 99
weight2_gain<-rep(0,length(Data$wtgain))
weight2_gain[weight_gain==98]<-1


natal<-Data$precare  #remove 5
natal2<-rep(0,length(natal))
natal2[natal==2]<-1
natal3<-rep(0,length(natal))
natal3[natal==3]<-1
novisit<-rep(0,length(natal))
novisit[natal==4]<-1
cig1<-Data$cig_1 #remove 99
cig2<-Data$cig_2 #remove 99
cig3<-Data$cig_3 #remove 99
sex<-as.numeric(Data$sex)
sex[sex==1]<-3
sex[sex==2]<-1
sex[sex==3]<-0 #1=boy
rm(Data)
#Matrix of observations
observations<-cbind(weight,rep(1,length(weight)), married,age,ed1,ed2,ed3,weight_gain, natal2,natal3,novisit,cig1,cig2,cig3,sex, white,black, age^2,weight2_gain)
rm(weight,  married,age,ed1,ed2,ed3,weight_gain, natal2,natal3,novisit,cig1,cig2,cig3,sex, white,black,  weight2_gain)
#remove missing observations
observations<-observations[observations[,1]!=9999,]
observations<-observations[observations[,7]>-1,]
observations<-observations[observations[,8]!=99,]
observations<-observations[observations[,12]!=99,]
observations<-observations[observations[,13]!=99,]
observations<-observations[observations[,14]!=99,]
observations<-observations[is.na(observations[,3])==FALSE,]
nrow(observations)


#write.table(observations, "Data/application.txt")
#####Download data DATA####################################################################
observations<-read.table("Data/application.txt")
#rescale observations
scale_values<-rep(1,ncol(observations))
mean_values<-rep(0,ncol(observations))
for(i in 3:(ncol(observations))){
  mean_values[i]<-mean(observations[,i])
  scale_values[i]<-2*sd(observations[,i])
  observations[,i]<-(observations[,i]-mean_values[i])/scale_values[i]

}
scale_values[1]<-1000
observations[,1]<-observations[,1]/scale_values[1]

d<-ncol(observations)-1
datalength<-nrow(observations)

#####################################################################################
##Choose (h_t) where h_0=0 and h_t=c_h t^{-alpha} when h_t>0
#####################################################################################
alpha<-0.3
c_h<-3
h_seq<-c(c_h*(1:datalength)^(-alpha))
#####################################################################################
#Generate (t_p) 
#####################################################################################
scale1<-1
diff_min<-1
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
 


 


#####################################################################################
## Estimation 
#####################################################################################

Sigma_t<-diag(1,d)

#Sigma_t[2:d,2:d]<-cov.wt(observations[1:20000,3:(d+1)], cor=TRUE)$cor
#work<-lower.tri(Sigma_t, diag = TRUE)
#Sigma_t<-vec2symMat(c(Sigma_t[work==TRUE]), diag = TRUE)



tuning_parameters <- list(N=50000, K=10, J=1, h2_t=h_seq^2,  T_p=pert_times , student=stud_times, Sigma=Sigma_t, 
                          nu=50,  ESS_bound=0.7, resampling=SSP_Resampling) 
  
 
theta_star=0
T_end<-50000
#set.seed(30485) 
start_time <- Sys.time()
theta   <- online(tuning_parameters ,target, observations[1:T_end,], theta_star, cl=numCores-2, PI_BAR=FALSE,  t_bar=1 )
end_time <- Sys.time()

T_star<-T_end+1
T_end2<-10^6
#set.seed(30485) 
theta2   <- online(tuning_parameters ,target, observations[T_star:T_end2,], theta_star, cl=numCores-2, 
                     particles_in=theta$TILDE[,1:d], W_in=theta$TILDE[,(d+1)], t_in=T_end )




val<-seq(800000,900000,1000)
res_val<-rbind(theta$MEAN,theta2$MEAN)
ess_val<-rbind(theta$ESS,theta2$ESS)


plot(val, ess_val[val]/tuning_parameters$N,type='l')


plot(val, res_val[val,k],type='l')




#constant, married,age,ed1,ed2,ed3,weight_gain, natal2,natal3,novisit,cig1,cig2,cig3,sex, white,black, age^2)

plot(T_star:T_end2, theta2$ESS/tuning_parameters$N,type='l')





 
plot(600000:900000 , theta2$MEAN[600000:900000,k],type='l')
 












