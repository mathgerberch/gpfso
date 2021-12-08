

library(gtools)
library(ggplot2)
library(scales)
library(reshape2)
library(gridExtra)
library(latex2exp)


###############################################################################################
#Table 1
###############################################################################################
g<-function(x) x/(1+x)
pi.tilde<-function(c, alpha,y,sigma0=5,mu0=0, gamma=1){
   n<-length(y)
   h_t<-c*c(0,(1:n)^{-alpha})
   h2_t<-h_t^2
   var<-rep(0,n+1)
   mean<-rep(0,n+1)
   var.post<-rep(0,n+1)
   mean.post<-rep(0,n+1)
   mmd_vec<-rep(0,n)
   var[1]<-sigma0^2
   mean[1]<-mu0^2
   var.post[1]<-sigma0^2
   mean.post[1]<-mu0^2
   for(t in 2:(n+1)){
      var[t]<-g(var[t-1]+h2_t[t-1])
      mean[t]<-var[t]*y[t-1]+(1-var[t])*mean[t-1] 
      var.post[t]<-g(var.post[t-1])
      mean.post[t]<-var.post[t]*y[t-1]+(1-var.post[t])*mean.post[t-1]     
   }
   return(list(MEAN=mean[2:(n+1)], VAR=var[2:(n+1)], MEAN.P=mean.post[2:(n+1)], VAR.P=var.post[2:(n+1)]))
}
 

mu_star<-0
sigma_star<-1

set.seed(9585)
n<-10^7
y<-rnorm(n,mu_star,sigma_star)

#run code below for alpha\in\{0.1,0.3,0.5,0.7,1\}
res<-pi.tilde(c=1, alpha=1,y)

start<-10^5
error<- log(abs(res$MEAN-mu_star))
lm(error[start:n] ~log(start:n))$coef


thetaB<-cumsum(res$MEAN[start:n])/(start:n)
errorB<- log(abs(thetaB-mu_star))
lm(errorB ~log(start:n))$coef


###############################################################################################
#CQR model
###############################################################################################

#####figure 1(a)#################
S1<-read.table("Data/CQR/ADA_BoxplotS_d5_q50.txt")$V1
S2<-read.table("Data/CQR/ADA_BoxplotS_d5_q99.txt")$V1


all<-data.frame(dataset=c(rep('1',length(S1))   ,rep('2',length(S2))),value=c(S1,S2))

p1<-ggplot(data = all, aes(x=dataset, y=value)) + geom_boxplot(aes(fill=dataset))+theme_bw()+
    labs(y = " ") + theme(legend.position = 'none')+xlab(" ")+
    scale_fill_manual(values=c("white","white","white","white"))+
    theme(axis.text=element_text(size=30, colour="black"),
        axis.title=element_text(size=35))+theme(legend.text=element_text(size=25))+
    #scale_y_log10(limits = c(0.04,0.16), breaks = c( 10^{-1.5},10^{-1.25}, 10^{-1}),labels = trans_format("log10", math_format(10^.x)))+
    ylab("estimation error") +scale_x_discrete(labels=c('1' = expression(tau==0.5), '2'=expression(tau==0.99)))
    
#figure 1(a)
p1

##########figure 1(b)##################################
errorE<-read.table("Data/CQR/errorE_d5_q5_M20.txt")$V1  
errorE_jit<-read.table("Data/CQR/errorE_d5_q5_jit_M20.txt")$V1  

tr<-read.table("Data/CQR/tr_d5.txt")$V1 #time instants at which above quantities are computed
T_end<-length(tr) 


error<-errorE
error_jit<-errorE_jit

beta<-0.5
work<-(tr)^{-beta}
c<-error[T_end]/work[T_end]
work<-c*work

df<- data.frame(Method= factor(c(rep("1", T_end),rep("2", T_end),rep("3", T_end) )),
		axis	 	= c(tr),
                res 		= c(error,error_jit, work )
		)


p1<-ggplot(data=df,  aes(x=axis, y=res, group=Method, linetype=Method)) + geom_line(size=1.5)+
	scale_x_log10(limits = c(1,max(tr)),breaks = c(10^1,10^{3},10^5,10^7  ), labels = trans_format("log10", math_format(10^.x)))+
	scale_y_log10(limits = c(10^{-2.5}, 10^1.5), breaks = c(10^{-2},10^{-1},10^0,10^1   ),labels = trans_format("log10", math_format(10^.x)))+
	xlab("iteration t \n (true observations)") +ylab("estimation error") +theme_bw()+
        theme(legend.title=element_blank())+theme(legend.position="bottom")+
	theme(axis.text=element_text(size=30, colour="black"),
        axis.title=element_text(size=35))+theme(legend.text=element_text(size=30))+
	scale_linetype_manual(values=c(1, 2,3))+
        scale_size_manual(values=c(0.1,2))+
	theme(legend.key.size = unit(1.2, "cm"))+ theme(legend.position = 'none')+
        annotate("text", x = 10^{0.7}, y = 10^{-0.2},  label = expression(t^{-0.5}), size=15)


#figure 1(b)
p1


##########figure 1(c)##################################
errorE<-read.table("Data/CQR/errorE_d5_q99_M20.txt")$V1 #Euclidian norm
errorE3<-read.table("Data/CQR/errorE_d5_q99_M20_a3.txt")$V1 #Euclidian norm


tr<-read.table("Data/CQR/tr_d5.txt")$V1 #time instants at which above quantities are computed
T_end<-length(tr) 


error<-errorE
error3<-errorE3

beta<-0.5
work<-(tr)^{-beta}
c<-error[T_end]/work[T_end]
work<-c*work


beta<-0.3
work3<-(tr)^{-beta}
c<-error3[T_end]/work3[T_end]
work3<-c*work3



df<- data.frame(Method= factor(c(rep("1", T_end),rep("2", T_end),rep("3", T_end),rep("4", T_end) )),
		axis	 	= c(tr),
                res 		= c(error,errorE3, work,work3 )
		)



p1<-ggplot(data=df,  aes(x=axis, y=res, group=Method, linetype=Method)) + geom_line(size=1.5)+
	scale_x_log10(limits = c(1,max(tr)),breaks = c(10^1,10^{3},10^5,10^7  ), labels = trans_format("log10", math_format(10^.x)))+
	scale_y_log10(limits = c(10^{-2.5}, 10^1.5), breaks = c(10^{-2}, 10^{-2},10^{-1},10^{0},10^1),labels = trans_format("log10", math_format(10^.x)))+
	xlab("iteration t \n (true observations)") +ylab("estimation error") +theme_bw()+
        theme(legend.title=element_blank())+theme(legend.position="bottom")+
	theme(axis.text=element_text(size=30, colour="black"),
        axis.title=element_text(size=35))+theme(legend.text=element_text(size=30))+
	scale_linetype_manual(values=c(1, 2,3,3))+
        scale_size_manual(values=c(0.1,2))+
	theme(legend.key.size = unit(1.2, "cm"))+ theme(legend.position = 'none')+
         annotate("text", x = 10^{4.3}, y = 10^{-1.3},  label = expression(t^{-0.5}), size=15)+
         annotate("text", x = 10^{0.5}, y = 10^{0.4},  label = expression(t^{-0.3}), size=15)


#figure 1(c)
p1


###############################################################################################
#Toy Multimodal Example
###############################################################################################

###########figure 2(a)##################################

S1<-read.table("Data/Multimodal/OPSMC_boxplotS_d20_M100_T10p5.txt")$V1
S2<-read.table("Data/Multimodal/boxplotS_d20_M100_T10p5.txt")$V1
S3<-read.table("Data/Multimodal/boxplotS_d20_M100_T10p5_var3.txt")$V1
S4<-read.table("Data/Multimodal/boxplotS_d20_M100_T10p5_var1.txt")$V1
S5<-read.table("Data/Multimodal/boxplotS_d20_M100_T10p5_var1_stud.txt")$V1
 
 
all<-data.frame(dataset=c(rep('1',length(S1)) ,rep('2',length(S2)),rep('3',length(S3))  ,rep('4',length(S4)) ,rep('5',length(S5))),value=c(S1,S2,S3,S4,S5))

p1<-ggplot(data = all, aes(x=dataset, y=value)) + geom_boxplot(aes(fill=dataset))+theme_bw()+
    labs(y = " ") + theme(legend.position = 'none')+xlab(" ")+
    scale_fill_manual(values=c("white","white","white","white","white"))+
    theme(axis.text=element_text(size=30, colour="black"),
        axis.title=element_text(size=35))+theme(legend.text=element_text(size=25))+
    scale_y_log10(limits = c(0.03,5), breaks = c(10^{-2}, 10^{-1},10^{0}, 10),labels = trans_format("log10", math_format(10^.x)))+
    ylab("estimation error")  +scale_x_discrete(labels=c('1' = expression(theta["K,T'"]^N), 
          '2'= expression(bar(theta)["T'"]^N), '3'=expression(bar(theta)["T'"]^N), '4'=expression(bar(theta)["T'"]^N), '5'=expression(bar(theta)["T'"]^N)))
   
 
#figure 2(a)
p1


##########figure 2(b)##################################

E1<-read.table("Data/Multimodal/traj_theta14_d20.txt")$V1
T_end<-30000
points<-15000:T_end

df<- data.frame(Method= factor(c(rep("1", length(points)))),
		axis	 	= c(points),
                res 		= c(E1[points])
		)



p1<-ggplot(data=df,  aes(x=axis, y=res, group=Method, linetype=Method)) + geom_line(size=0.2)+
	xlab("iteration t (in thousands) \n (true observations)") +ylab(expression(theta[14])) +theme_bw()+#scale_x_continuous(limits=c(1000,T_end))+
	scale_x_continuous(breaks = c(15000,20000, 25000, 30000), labels=c("15","20","25", "30"))+
        theme(legend.title=element_blank())+theme(legend.position="bottom")+
	theme(axis.text=element_text(size=30, colour="black"),
        axis.title=element_text(size=35))+theme(legend.text=element_text(size=30))+
	scale_linetype_manual(values=c(1, 3))+
        scale_size_manual(values=c(0.1,2))+
	theme(legend.key.size = unit(1.2, "cm"))+ theme(legend.position = 'none')+geom_hline(yintercept=-1)
     
 
#figure 2(b)
p1


###########figure 2(c)##################################
errorE<-read.table("Data/Multimodal/errorE_d20_M10.txt")$V1 #Euclidian norm
tr<-read.table("Data/Multimodal/tr_d20.txt")$V1 #time instants at which above quantities are computed
T_end<-length(tr) 
error<-errorE
error2<-errorE2
error3<-errorE3


beta<-0.5
work<-(tr)^{-beta}
c<-error[T_end]/work[T_end]
work<-c*work

 

df<- data.frame(Method= factor(c(rep("1", T_end),rep("3", T_end) )),
		axis	 	= c(tr),
                res 		= c(error, work)
		)



p1<-ggplot(data=df,  aes(x=axis, y=res, group=Method, linetype=Method)) + geom_line(size=1.5)+
	scale_x_log10(limits = c(1,max(tr)),breaks = c(10^1,10^{3},10^5 ), labels = trans_format("log10", math_format(10^.x)))+
	scale_y_log10(limits = c(0.08, 10^2), breaks = c( 10^{-1},10^{0},10^1, 10^2   ),labels = trans_format("log10", math_format(10^.x)))+
	xlab("iteration t \n (true observations)") +ylab("estimation error") +theme_bw()+
        theme(legend.title=element_blank())+theme(legend.position="bottom")+
	theme(axis.text=element_text(size=30, colour="black"),
        axis.title=element_text(size=35))+theme(legend.text=element_text(size=30))+
	scale_linetype_manual(values=c(1, 3,2,5))+
        scale_size_manual(values=c(0.1,2))+
	theme(legend.key.size = unit(1.2, "cm"))+ theme(legend.position = 'none')+
         annotate("text", x = 10^{0.7}, y = 10^{1},  label = expression(t^{-0.5}), size=15)
 
#figure 2(c)
p1


###############################################################################################
#SAGM model
###############################################################################################

###########figure 3(a)##################################

S1<-read.table("Data/Mixture/BoxplotS_M100.txt")$V1
tS1<-read.table("Data/Mixture/tBoxplotS_M100.txt")$V1
S2<-read.table("Data/Mixture/BoxplotS_M100_T10p5.txt")$V1
tS2<-read.table("Data/Mixture/tBoxplotS_M100_T10p5.txt")$V1
tS3<-read.table("Data/Mixture/tBoxplotS_M100_st.txt")$V1 
tS4<-read.table("Data/Mixture/tBoxplotS_M100_T10p5_st.txt")$V1


 
all<-data.frame(dataset=c(rep('1',length(S1)),rep('2',length(tS1)),rep('3',length(S2)), rep('4',length(tS2)),
  rep('5',length(tS3)),  rep('6',length(tS4))),value=c(S1,tS1,S2,tS2,tS3,tS4))


p1<-ggplot(data = all, aes(x=dataset, y=value)) + geom_boxplot(aes(fill=dataset))+theme_bw()+
    labs(y = " ") + theme(legend.position = 'none')+xlab(" ")+
    scale_fill_manual(values=c("white","white","white","white","grey","grey"))+theme( axis.text.x = element_text(size = 30))+
    theme(axis.text=element_text(size=30, colour="black"),
        axis.title=element_text(size=35))+theme(legend.text=element_text(size=25))+
    scale_y_log10(limits = c(0.02,15), breaks = c( 10^{-1.5},10^{-0.5},   10^0.5),labels = trans_format("log10", math_format(10^.x)))+
    ylab("estimation error")   +scale_x_discrete(labels=c('1'= expression(bar(theta)[T[1]]^N), '2'=expression(tilde(theta)[T[1]]^N), 
    '3'=expression(bar(theta)[T[2]]^N),'4'=expression(tilde(theta)[T[2]]^N),'5'= expression(tilde(theta)[T[1]]^N),
           '6'= expression(tilde(theta)[T[2]]^N)))
    
#figure 3(a)
p1


M<-length(tS1)
m1<-rep(0,M)
m2<-rep(0,M)
m1[tS1>2.16]<-1
m2[tS2<0.23]<-1
sum(m1*m2) 


######figure 3(b)#################
theta_star<-read.table("Data/Mixture/theta_star.txt")$V1
E1<-read.table("Data/Mixture/traj_6.txt")$V1
k<-6
E1<-E1[1:100000]
T_end<-length(E1)
points<-1:T_end

df<- data.frame(Method= factor(c( rep("1", length(points)))),
		axis	 	= c(points),
                res 		= c(E1[points])
		)



p1<-ggplot(data=df,  aes(x=axis, y=res, group=Method, linetype=Method)) + geom_line(size=1)+
	xlab("iteration t (in thousands) \n (true observations)") +ylab(expression(theta[6])) +theme_bw()+
	scale_x_continuous(breaks = c(0, 20000,40000,60000,80000, 100000), labels=c("0","20","40", "60", "80","100"))+
        theme(legend.title=element_blank())+theme(legend.position="bottom")+
	theme(axis.text=element_text(size=30, colour="black"),
        axis.title=element_text(size=35))+theme(legend.text=element_text(size=30))+
	scale_linetype_manual(values=c(1, 2))+geom_hline(yintercept=theta_star[k])+
        scale_size_manual(values=c(0.1,2))+
	theme(legend.key.size = unit(1.2, "cm"))+ theme(legend.position = 'none')
         

#figure 3(b)
p1


######figure 3(c)#################

error<-read.table("Data/Mixture/errorE_M5.txt")$V1
tr<-read.table("Data/Mixture/tr.txt")$V1
 
start<-5000
error<-error[start:length(tr)]
tr<-tr[start:length(tr)]


beta<-0.5
work<-(tr)^{-beta}
c<-error[length(tr)]/work[length(error)]
work<-c*work

 
 
df<- data.frame(Method= factor(c(rep("1", length(tr)),rep("3", length(tr)))),
		axis	 	= c(tr),
                res 		= c(error, work)
		)

p1<-ggplot(data=df,  aes(x=axis, y=res, group=Method, linetype=Method)) + geom_line(size=1.5)+
	scale_x_log10(limits = c(min(tr),max(tr)),breaks = c(10^5, 10^5.5,  10^6 ), labels = trans_format("log10", math_format(10^.x)))+
	scale_y_log10(limits = c(0.02, 0.13), breaks = c(10^{-1.75}, 10^{-1.5} ,10^{-1.25},  10^{-1}, 10^{-0.75}),labels = trans_format("log10", math_format(10^.x)))+
	xlab("iteration t \n (true observations)") +ylab("estimation error") +theme_bw()+
        theme(legend.title=element_blank())+theme(legend.position="bottom")+
	theme(axis.text=element_text(size=30, colour="black"),
        axis.title=element_text(size=35))+theme(legend.text=element_text(size=30))+
	scale_linetype_manual(values=c(1, 3,2,5))+
        scale_size_manual(values=c(0.1,2))+
	theme(legend.key.size = unit(1.2, "cm"))+ theme(legend.position = 'none')+
         annotate("text", x = 10^{5.16}, y = 10^{-1.3},  label = expression(t^{-0.5}), size=15)


#figure 3(c)
p1
 
###############################################################################################
#g and k distribution
###############################################################################################
######figure 4(a)#################

library(winference)
library(numDeriv)
data2<-read.table("Data/g_and_k/exchange_rate.txt") 
x1<-data2$V4
x2<-data2$V5
observations<-cbind(diff(log(x1)),diff(log(x2)))
datalength<-nrow(observations)
observations<-100*observations
rm(data2,x1,x2)


theta_star<-read.table("Data/g_and_k/mle.txt")$V1
         
                   
cdf1_ <- function(z)sapply(z, FUN = function(v)gandkcdf(v, theta_star[1:4]))
pdf1 <- function(x) {return(grad(cdf1_, x))}           
cdf2_ <- function(z)sapply(z, FUN = function(v)gandkcdf(v, theta_star[5:8]))
pdf2 <- function(x) {return(grad(cdf2_, x))}


                  

move1<- -2
move2<-6
df <- data.frame(data1 = move1+observations[,1], data2 = move2+observations[,2],   curve1=pdf1(observations[,1]),curve2=pdf2(observations[,2]))


p1 <- ggplot(df) +  theme_bw()+labs(y = "density") + theme(legend.position = 'none')+xlab(" ")+#ylim(0,4)+
    geom_histogram(aes(x = data1, y = ..density..),
                   bins=1000, fill = "white", color = "gray60")+
                   geom_line(data = df, aes(x = move1+observations[,1], y = curve1))+
       geom_histogram(aes(x = data2, y = ..density..),
                   bins=1000, fill = "white", color = "gray60")+
                   geom_line(data = df, aes(x = move2+observations[,2], y = curve2))+
      theme(axis.text=element_text(size=30, colour="black"),
        axis.title=element_text(size=35))+theme(legend.text=element_text(size=25))+
        theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+coord_cartesian(xlim = c(-5,10), ylim=c(0,6))



theta_val<-read.table("Data/g_and_k/local_maxima.txt") 
theta_val<-theta_val$V2        
 
                   

cdf3_ <- function(z)sapply(z, FUN = function(v)gandkcdf(v, theta_val[1:4]))
pdf3 <- function(x) {return(grad(cdf3_, x))}
cdf4_ <- function(z)sapply(z, FUN = function(v)gandkcdf(v, theta_val[5:8]))
pdf4 <- function(x) {return(grad(cdf4_, x))}

move3<- -2
move4<-6
df2 <- data.frame(data3 = move3+observations[,1], data4 = move4+observations[,2],   curve3=pdf3(observations[,1]),curve4=pdf4(observations[,2]))


p2 <- ggplot(df2) +  theme_bw()+labs(y = "density") + theme(legend.position = 'none')+xlab(" ")+#ylim(0,4)+
    geom_histogram(aes(x = data3, y = ..density..),
                   bins=1000, fill = "white", color = "gray60")+
                   geom_line(data = df2, aes(x = move3+observations[,1], y = curve3))+
       geom_histogram(aes(x = data4, y = ..density..),
                   bins=1000, fill = "white", color = "gray60")+
                   geom_line(data = df2, aes(x = move4+observations[,2], y = curve4))+
      theme(axis.text=element_text(size=30, colour="black"),
        axis.title=element_text(size=35))+theme(legend.text=element_text(size=25))+
        theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +coord_cartesian(xlim = c(-5,10), ylim=c(0,6)) 
        
                   
#figure 4(a)                    
grid.arrange(p1,p2, nrow=1)
                         
######figure 4(b)#################
L1<-read.table("Data/g_and_k/errorS_Boxplot_nit10000.txt",)$V1
L2<-read.table("Data/g_and_k/errorS_tBoxplot_nit10000.txt")$V1
L3<-read.table("Data/g_and_k/errorS_Boxplot.txt",)$V1
L4<-read.table("Data/g_and_k/errorS_tBoxplot.txt",)$V1

L5<-read.table("Data/g_and_k/errorS_Boxplot_alpha05.txt",)$V1
L6<-read.table("Data/g_and_k/errorS_Boxplot_alpha05_c1_N500.txt")$V1
      
all<-data.frame(dataset=c(rep('1',length(L1)) ,rep('2',length(L2)), rep('3',length(L3)), rep('4',length(L4))   , 
                         rep('5',length(L5))   ,rep('6',length(L6))),value=c(L1,L2,L3,L4,L5,L6))


p1<-ggplot(data = all, aes(x=dataset, y=value)) + geom_boxplot(aes(fill=dataset))+theme_bw()+
    labs(y = " ") + theme(legend.position = 'none')+xlab(" ")+
    scale_fill_manual(values=c("white","white","white","white","black","grey"))+theme( axis.text.x = element_text(size = 30))+
    theme(axis.text=element_text(size=30, colour="black"),
        axis.title=element_text(size=35))+theme(legend.text=element_text(size=25))+
    scale_y_log10(limits = c(0.01,10), breaks = c(10^{-2}, 10^{-1},10^{0},   10^1),labels = trans_format("log10", math_format(10^.x)))+
    ylab("estimation error")   +scale_x_discrete(labels=c('1'= expression(bar(theta)[T[1]]^N), '2'=expression(tilde(theta)[T[1]]^N), 
    '3'=expression(bar(theta)[T[2]]^N),'4'=expression(tilde(theta)[T[2]]^N), '5'= expression(bar(theta)[T[2]]^N), '6'= expression(bar(theta)[T[2]]^N)))
        
#figure 4(b)
p1


######figure 4(c)#################

errorE<-read.table("Data/g_and_k/errorE_c1_500.txt")$V1 #Euclidian norm
 
tr<-read.table("Data/g_and_k/tr.txt")$V1 #time instants at which above quantities are computed
T_end<-length(tr) 

error<-errorE


beta<-0.5
work<-(tr)^{-beta}
c<-error[T_end]/work[T_end]
work<-c*work



df<- data.frame(Method= factor(c(rep("1", T_end),rep("2", T_end) )),
		axis	 	= c(tr),
                res 		= c(error, work )
		)



p1<-ggplot(data=df,  aes(x=axis, y=res, group=Method, linetype=Method)) + geom_line(size=1.5)+
	scale_x_log10(limits = c(min(tr),max(tr)),breaks = c(10^{4},10^5,10^6  ), labels = trans_format("log10", math_format(10^.x)))+
	scale_y_log10(limits = c(10^{-1.5}, 0.3), breaks = c(10^{-1.5}, 10^{-1.25}, 10^{-1},10^{-0.75}, 10^{-0.5}  ),labels = trans_format("log10", math_format(10^.x, format=force)))+
	xlab("iteration t \n (pseudo-observations)") +ylab("estimation error") +theme_bw()+
        theme(legend.title=element_blank())+theme(legend.position="bottom")+
	theme(axis.text=element_text(size=30, colour="black"),
        axis.title=element_text(size=35))+theme(legend.text=element_text(size=30))+
	scale_linetype_manual(values=c( 1,3))+
        scale_size_manual(values=c(0.1,2))+
	theme(legend.key.size = unit(1.2, "cm"))+ theme(legend.position = 'none')+ 
        annotate("text", x = 10^{5.2}, y = 10^{-0.72},  label = expression(t^{-0.5}), size=15)
 
#figure 4(c)
p1
       
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   




