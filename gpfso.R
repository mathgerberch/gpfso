



##compute (h_t) and (t_p)
h_Tp<-function(alpha,datalength, A=1, B=1, varrho=0.1, t_0=5){
   h_seq<-(1:datalength)^(-alpha)
   stud_YN<-rep(0,datalength)
   t_1<-1
   while(t_1 <= datalength){
     t_1<-t_0+ceiling(max(A*t_0^varrho*log(t_0),B))
     if(t_1 <= datalength) stud_YN[t_1]<-1
     t_0<-t_1
  }
  return(list(H=h_seq,Tp=stud_YN))
}


###########################################################################################################################
# G_PFSO ALGORITHM  
###########################################################################################################################


rmutation_pi<- function(particles,t, target_parameters){
    if(target_parameters$student[t]==1){
       return(particles+rmt(nrow(particles), rep(0,ncol(particles)), target_parameters$h2_t[t]*target_parameters$Sigma, df=target_parameters$nu))
    }else{ 
       return(particles+rmnorm(nrow(particles),rep(0,ncol(particles)),  target_parameters$h2_t[t]*target_parameters$Sigma))
    }
}

G_PFSO<-function(N, target, observations, cl=1, t_in=1, c_ess=0.7){
  #//Input parameters
  ESS_bound<- N*c_ess             
  d<-target$dimension
  datalength <- nrow(observations)
  data_col<-ncol(observations)
  #//Initialization
  w<-rep(0,N)
  particles<-as.matrix(target$rprior(N, target$parameters)) ##sample from prior
   #To store useful quantities
  mean_vec<-matrix(0,(datalength-t_in+1),d)
  ESS_vec<-rep(0,(datalength-t_in+1))
 
  t<-t_in
  work<-target$loglikelihood(particles,  matrix(observations[t,], ncol=data_col),target$parameters, cl)
  work[is.nan(work)]<--Inf
  w<-work
  w1<- exp(w - max(w))
  W<- w1 / sum(w1)
  ESS<-1/sum(W^2)
  mean_vec[t-t_in+1,]<-apply(W*particles,2,sum)

  for(t in  (t_in+1):datalength){
    ESS_vec[t-t_in]<-ESS
    if(ESS<ESS_bound){
       A<-SSP_Resampling(W) 
       particles<-rmutation_pi(as.matrix(particles[A,]), t-1,  target$parameters)
       w<-target$loglikelihood(particles, matrix(observations[t,], ncol=data_col), target$parameters, cl)
       w[is.nan(w)]<--Inf
       w1<- exp(w - max(w))
       W<- w1 / sum(w1)
       ESS<-1/sum(W^2)
    }else{
       particles<-rmutation_pi(particles,t-1, target$parameters)
       work<-target$loglikelihood(particles,  matrix(observations[t,], ncol=data_col),target$parameters, cl)
       work[is.nan(work)]<--Inf
       w<-w+work
       w1<- exp(w - max(w))
       W<- w1 / sum(w1)
       ESS<-1/sum(W^2)
    }
    mean_vec[t-t_in+1,]<-apply(W*particles,2,sum)
  }
  return(list(MEAN=mean_vec,  ESS=ESS_vec))


}


###########################################################################################################################
#OPSMC
###########################################################################################################################



OPSMC<-function(tuning_parameters,  target, observations,scale_var=1, alpha=-1, cl=1){
  N<- tuning_parameters$N
  d<-target$dimension
  ESS_bound<- N*tuning_parameters$ESS_bound   
  if(alpha< 0){
      bN <- (4/((d+2)*N))^{2/(d+4)}
  }else{
      bN<- alpha
  }      
  rho<-sqrt(1-bN)
  d<-target$dimension
  datalength <- nrow(observations)
  data_col<-ncol(observations)
  #//Initialization
  w<-rep(0,N)
  particles<-as.matrix(target$rprior(N, target$parameters)) ##sample from prior
  w<-rep(0,N)
  W<-rep(1/N,N)
  
  ESS<-1/sum(W^2)  
  mean_vec<-matrix(0,datalength,d)
  ess_vec<-rep(0,datalength)
 
  for(t in 1:datalength){
    ess_vec[t]<-ESS
    if(ESS<ESS_bound && t>1){
       A<-tuning_parameters$resampling(W) 
       if(ESS==1){
         Cov_mat<-diag(scale_var,d)
       }else{
         Cov_mat<- cov.wt(particles, W)$cov 
         Cov_mat[lower.tri(Cov_mat)] = t(Cov_mat)[lower.tri(Cov_mat)]
         if(min(eigen(Cov_mat)$values)<=10^{-7})  Cov_mat<-diag(diag(Cov_mat),d)
      } 
      Cov_mat<-bN*Cov_mat
       
      particles<- t(t(rho*as.matrix(particles[A,]))+(1-rho)*mean_vec[t-1,])+rmnorm(N, rep(0,d),Cov_mat)
      lik<-target$loglikelihood(particles, matrix(observations[t,], ncol=data_col),target$parameters, cl)
      lik[is.nan(lik)]<--Inf
      w<-lik
      w1<- exp(w - max(w))
      W<- w1 / sum(w1)
      ESS<-1/sum(W^2)
    }else{
       work<-target$loglikelihood(particles, matrix(observations[t,], ncol=data_col),target$parameters, cl)
       work[is.nan(work)]<--Inf
       w<-w+ work 
       w1<- exp(w - max(w))
       W<- w1 / sum(w1)
       ESS<-1/sum(W^2)
    }
    mean_vec[t,]<-apply(W*particles,2,sum)

  }
  return(list(MEAN=mean_vec, ESS=ess_vec))


}















