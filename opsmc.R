


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
       A<-SSP_Resampling(W) 
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















