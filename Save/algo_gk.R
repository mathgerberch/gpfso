











IBIS<-function(tuning_parameters, target, observations, cl=1){
  #//Input parameters
  N<- tuning_parameters$N
  ESS_bound<- N*tuning_parameters$ESS_bound   
  acc_rate<-{}              
  d<-target$dimension
  datalength <- nrow(observations)
  data_col<-ncol(observations)

  #//Initialization
  w<-rep(0,N)
  particles_0<-as.matrix(target$rprior(N, target$parameters)) ##sample from prior
  particles<-particles_0
  w<-rep(0,N)
  W<-rep(1/N,N)
  W_0<-rep(1/N,N)
  ESS<-1/sum(W^2)
  lik0<-rep(0,N)
  particles_bar<-matrix(0,N,d)
  mu_vec<-0
  Cov_mat<-0
  proposals<-matrix(0,N,d)
  w1<-rep(0,N)
  lik1<-0
  kt<-0
  scale<-(2.38*d^{-0.5})^2
  mean_vec<-matrix(0,datalength,d)
  mean2_vec<-matrix(0,datalength,d)
  ESS_vec<-rep(0,datalength)
  ESS_times<-{}

  for(t in 1:datalength){
    ESS_vec[t]<-ESS
    if(ESS<ESS_bound && t>1){
       ESS_times<-c(ESS_times,t-1)
       A<-tuning_parameters$resampling(W) 
       Cov_mat<-scale*cov.wt(particles, W)$cov 
       Cov_mat[lower.tri(Cov_mat)] = t(Cov_mat)[lower.tri(Cov_mat)]
       if(min(eigen(Cov_mat)$values)<=10^{-7})  Cov_mat<-diag(diag(Cov_mat),d)
       proposals<-as.matrix(particles[A,])+rmnorm(N, rep(0,d),Cov_mat)
       lik1<-target$loglikelihood(proposals, matrix(observations[(kt+1):(t-1),], ncol=data_col),target$parameters, cl)
       lik1[is.nan(lik1)]<--Inf
       if(mean(abs(lik0))==0){
              lik0<-target$loglikelihood(particles, matrix(observations[(kt+1):(t-1),], ncol=data_col), target$parameters, cl)
              lik0[is.nan(lik0)]<--Inf
       }
       alpha<-lik1-lik0[A]+target$dprior(proposals,target$parameters)-target$dprior(as.matrix(particles[A,]),target$parameters)
       U<-log(runif(N))
       work<- U<=alpha
       acc_rate<-c(acc_rate,sum(work)/N)
       particles<-as.matrix(particles[A,])
       particles[work,]<-proposals[work,] 
       lik0<-lik0[A]
       lik0[work]<-lik1[work]
       w<-target$loglikelihood(particles,  matrix(observations[t,], ncol=data_col),target$parameters, cl)
       w[is.nan(w)]<--Inf
       w1<- exp(w - max(w))
       W<- w1 / sum(w1)
       lik0<-lik0+w
       ESS<-1/sum(W^2)
    }else{
       work<-target$loglikelihood(particles, matrix(observations[t,], ncol=data_col),target$parameters, cl)
       work[is.nan(work)]<--Inf
       lik0<-lik0+work
       w<-w+ work 
       w1<- exp(w - max(w))
       W<- w1 / sum(w1)
       ESS<-1/sum(W^2)
    }
    mean_vec[t,]<-apply(W*particles,2,sum)
    mean2_vec[t,]<-apply(W*particles^2,2,sum)
  }


}

online_gk<-function(tuning_parameters, target, observations, mmd_times, cl=1, PI_BAR=FALSE,  t_bar=1){
  #//Input parameters
  N<- tuning_parameters$N
  J<- tuning_parameters$J
  K<- tuning_parameters$K
  ESS_bound<- N*tuning_parameters$ESS_bound   
  acc_rate<-{}      
  lik_eval<-{}        
  T_p <- c(target$parameters$T_p, datalength+1)
  d<-target$dimension
  datalength <- nrow(observations)
  data_col<-ncol(observations)

  #//Initialization
  w<-rep(0,N)
  particles_0<-as.matrix(target$rprior(N, target$parameters)) ##sample from prior
  particles<-particles_0
  w<-rep(0,N)
  W<-rep(1/N,N)
  W_0<-rep(1/N,N)
  ESS<-1/sum(W^2)
  lik0<-rep(0,N)
  particles_bar<-matrix(0,N,d)
  mu_vec<-0
  Cov_mat<-0
  proposals<-matrix(0,N,d)
  w1<-rep(0,N)
  lik1<-0
  kt<-0
  work<-1:length(T_p)
  p<-min(work[T_p>0])
  scale<-(2.38*d^{-0.5})^2
  mean_vec<-matrix(0,datalength,d)
  mean2_vec<-matrix(0,datalength,d)
  ESS_vec<-rep(0,datalength)
  
  tilde_mmd<-array(0,dim=c(length(mmd_times),N,d+1))
  bar_mmd<-{}
  if(PI_BAR==TRUE) bar_mmd<-array(0,dim=c(length(mmd_times),K,d))
  count_mmd<-1

  for(t in 1:datalength){
    ESS_vec[t]<-ESS
    if(ESS<ESS_bound && t>1){
       A<-tuning_parameters$resampling(W) 
       if(t-1==T_p[p]){  ##h_{t-1}>0
         particles_0<-target$rmutation(as.matrix(particles[A,]), t-1,  target$parameters)
         w<-target$loglikelihood(particles_0, matrix(observations[t,], ncol=data_col), target$parameters, cl)
         w[is.nan(w)]<--Inf
         w1<- exp(w - max(w))
    	 W_0<- w1 / sum(w1)
         particles<-particles_0
         W<-W_0 
         kt<-T_p[p]
         lik0<-rep(0,N)
         p<-p+1
       }else{ ##h_{t-1}=0 
         Cov_mat<-scale*cov.wt(particles, W)$cov 
         Cov_mat[lower.tri(Cov_mat)] = t(Cov_mat)[lower.tri(Cov_mat)]
         if(min(eigen(Cov_mat)$values)<=10^{-7})  Cov_mat<-diag(diag(Cov_mat),d)
         proposals<-as.matrix(particles[A,])+rmnorm(N, rep(0,d),Cov_mat)
         lik1<-target$loglikelihood(proposals, matrix(observations[(kt+1):(t-1),], ncol=data_col),target$parameters, cl)
         lik1[is.nan(lik1)]<--Inf
         if(mean(abs(lik0))==0){
                   lik0<-target$loglikelihood(particles, matrix(observations[(kt+1):(t-1),], ncol=data_col), target$parameters, cl)
                   lik0[is.nan(lik0)]<--Inf
         }
         if(kt==0){
             alpha<-lik1-lik0[A]+target$dprior(proposals,target$parameters)-target$dprior(as.matrix(particles[A,]),target$parameters)
         }else{
             if(J>1){
                A0<-matrix(sample(1:N,N*J,replace=T, prob=W_0),N,byrow=T)
                A1<-matrix(sample(1:N,N*J,replace=T, prob=W_0),N,byrow=T)
                work1<-0
                work2<-0
                for(j in 1:J){
                  work1<-work1+exp(target$dmutation(proposals-particles_0[A0[,j],], kt, target$parameters))
                  work2<-work2+exp(target$dmutation(as.matrix(particles[A,])-particles_0[A1[,j],], kt, target$parameters))
                }
                alpha<-lik1-lik0[A]+log(work1)-log(work2) 
                alpha[is.nan(alpha)==TRUE]<-log(0) 
             }else{ 
                A0<-sample(1:N,N,replace=T, prob=W_0)
                A1<-sample(1:N,N,replace=T, prob=W_0)
                alpha<-lik1-lik0[A]+
                target$dmutation(proposals-particles_0[A0,], kt, target$parameters)-
                target$dmutation(as.matrix(particles[A,])-particles_0[A1,], kt, target$parameters)
             }
          }
          U<-log(runif(N))
          work<- U<=alpha
          acc_rate<-c(acc_rate,sum(work)/N)
          lik_eval<-c(lik_eval,t-1)
          particles<-as.matrix(particles[A,])
          particles[work,]<-proposals[work,] 
          lik0<-lik0[A]
          lik0[work]<-lik1[work]
          w<-target$loglikelihood(particles,  matrix(observations[t,], ncol=data_col),target$parameters, cl)
          w[is.nan(w)]<--Inf
          w1<- exp(w - max(w))
    	  W<- w1 / sum(w1)
          lik0<-lik0+w
          #cat("ESS equals", ESS)
          #cat("time is", t)
       }
       ESS<-1/sum(W^2)
    }else{
       if(t-1==T_p[p]){
          particles_0<-target$rmutation(particles,t-1, target$parameters)
          work<-target$loglikelihood(particles_0,  matrix(observations[t,], ncol=data_col),target$parameters, cl)
          work[is.nan(work)]<--Inf
          w<-w+work
          w1<- exp(w - max(w))
    	  W_0<- w1 / sum(w1)
          particles<-particles_0
          W<-W_0
          kt<-T_p[p]
          lik0<-rep(0,N)
          p<-p+1
       }else{
          work<-target$loglikelihood(particles, matrix(observations[t,], ncol=data_col),target$parameters, cl)
          work[is.nan(work)]<--Inf
          lik0<-lik0+work
          w<-w+ work 
          w1<- exp(w - max(w))
    	  W<- w1 / sum(w1)
       }
       ESS<-1/sum(W^2)
    }
    mean_vec[t,]<-apply(W*particles,2,sum)
    mean2_vec[t,]<-apply(W*particles^2,2,sum)
    if(PI_BAR && t==t_bar){
        particles_bar<-particles[sample(1:N,K,replace=T, prob=W),]
    }
    if(PI_BAR && t>t_bar){
        V<-rbinom(K, 1, prob=1/(t-t_bar+1))
        if(sum(V)>0){
          work<-sample(1:N,sum(V),replace=T, prob=W) 
          particles_bar[V==1,]<-particles[work,]
        }
    } 
    if(t==mmd_times[count_mmd]){
         tilde_mmd[count_mmd,,]<-cbind(particles,W)
         if(PI_BAR){
             bar_mmd[count_mmd,,]<-particles_bar
         }
         count_mmd<-count_mmd+1
    }  
  }
  return(list(MEAN=mean_vec,MEAN2=mean2_vec, BAR=particles_bar, TILDE=cbind(particles,W), ESS=ESS_vec, ACC=acc_rate,TMMD=tilde_mmd, BMMD=bar_mmd, NLIK=lik_eval))


}





temp_smc<-function(tuning_parameters, target, observations){
  #//Input parameters
  N<- tuning_parameters$N
  ESS_bound<- N*tuning_parameters$ESS_bound   
  acc_rate<-{}          
  d<-target$dimension
  datalength <- nrow(observations)
  data_col<-ncol(observations)
  scale<-(2.38*d^{-0.5})^2
  #//Initialization
  particles<-as.matrix(target$rprior(N, target$parameters)) ##sample from prior 
  W<-rep(1/N,N)
  delta0<-0
  delta_seq<-delta0   
  
  lik0<-target$loglikelihood(particles, matrix(observations, ncol=data_col),target$parameters, cl) ##compute likelihood
  prior0<-target$dprior(as.matrix(particles),target$parameters)
  delta1<-delta_update(N,lik0, delta0, ESS_bound)       #compute delta1
  w<-(delta1-delta0)*lik0
  w[is.nan(w)]<--Inf
  w1<- exp(w - max(w))
  W<- w1 / sum(w1)   #estimate pi_delta1
  while(delta1<1){
    A<-tuning_parameters$resampling(W) 
    Cov_mat<-scale*cov.wt(particles, W)$cov 
    Cov_mat[lower.tri(Cov_mat)] = t(Cov_mat)[lower.tri(Cov_mat)]
    if(min(eigen(Cov_mat)$values)<=10^{-7})  Cov_mat<-diag(diag(Cov_mat),d)
    proposals<-as.matrix(particles[A,])+rmnorm(N, rep(0,d),Cov_mat)
    lik1<-target$loglikelihood(proposals, matrix(observations, ncol=data_col),target$parameters, cl)
    lik1[is.nan(lik1)]<--Inf
    prior1<-target$dprior(proposals,target$parameters)
    alpha<-(prior1-prior0[A])+delta1*(lik1-lik0[A])
    U<-log(runif(N))
    work<- U<=alpha
    acc_rate<-c(acc_rate,sum(work)/N)
    particles<-as.matrix(particles[A,])
    particles[work,]<-proposals[work,] 
    lik0<-lik0[A]
    lik0[work]<-lik1[work]
    prior0<-prior0[A]
    prior0[work]<-prior1[work] #resampled estimate of pi_delta1
    #update detla1
    delta0<-delta1
    delta1<-delta_update(N,lik0, delta0, ESS_bound)
    delta_seq<-c(delta_seq,delta0) 
    
    w<-(delta1-delta0)*lik0
    w[is.nan(w)]<--Inf
    w1<- exp(w - max(w))
    W<- w1 / sum(w1)  #estimate pi_delta1
  
    print(delta1)
    print(1/sum(W^2))
  }
  return(list(TILDE=cbind(particles,W), DELTA=delta_seq, ACC=acc_rate))
}


















