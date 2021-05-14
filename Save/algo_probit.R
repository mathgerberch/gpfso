




init_probit<-function(tuning_parameters, target, observations, t_in){
    N<- tuning_parameters$N
    Omega<-target$parameters$prior_covariance 
    d<-ncol(Omega)
    omega<-diag(sqrt(diag(Omega)),d)
    Omega_bar<-solve(omega)%*%Omega%*%solve(omega)
    particles<-matrix(rep(target$parameters$prior_mean,N),nrow=N, byrow=TRUE)
    val_probit<-const_probit(particles, matrix(observations[1:t_in,], ncol=(d+1)),target$parameters, omega, Omega, Omega_bar)
    return(rmutation_probit(1:N, particles, matrix(observations[1:t_in,], ncol=(d+1)), val_probit, omega))
}



online_probit<-function(tuning_parameters, target, observations, cl=1, PI_BAR=FALSE,  t_bar=1, particles_in,  t_in){
  #//Input parameters
  N<- tuning_parameters$N
  K<- tuning_parameters$K
  B<- tuning_parameters$B
  Omega<-target$parameters$Sigma
  d<-ncol(Omega)
  omega<-diag(sqrt(diag(Omega)),d)
  Omega_bar<-solve(omega)%*%Omega%*%solve(omega)
 
  
  T_p <- c(target$parameters$T_p, datalength+1)
  work<-1:length(T_p)
  p<-min(work[T_p>t_in])
  datalength <- nrow(observations)
  data_col<-ncol(observations)

  mean_vec<-matrix(0,datalength,d) 
  mean2_vec<-matrix(0,datalength,d)
  ESS_vec<-rep(0,datalength)
  particles_bar<-matrix(0,K,d)
  w<-rep(0,N)

  particles<-particles_in
  t_0<-t_in+1
  count<-0
  data_B<-floor(datalength)/B
  for(t in (t_in+1):data_B){
    
       if(target$parameters$student[t-1]==0 ||(t-1==t_in) ){  ##h_{t-1}>0
           val_probit<-const_probit(particles, matrix(observations[t_0:(t_0+B-1),], ncol=data_col),target$parameters, omega, Omega, Omega_bar)
           w<-w+val_probit$W
           w[is.nan(w)]<--Inf
           w1<- exp(w - max(w))
           W<- w1 / sum(w1)
           count<-count+1
           mean_vec[count,]<-apply(c(W)*particles,2,sum)
           mean2_vec[count,]<-apply(c(W)*particles^2,2,sum)
           ESS_vec[count]<-1/sum(W^2)
           A<-tuning_parameters$resampling(W) 
           particles<-rmutation_probit(A, particles, matrix(observations[t_0:(t_0+B-1),], ncol=data_col), val_probit, omega)
           p<-p+1
           t_0<-t_0+B
           w<-rep(0,N)
           if(PI_BAR && t==t_bar){
               particles_bar<-particles[sample(1:N,K,replace=T),]
           }
           if(PI_BAR && t>t_bar){
               V<-rbinom(K, 1, prob=1/(t-t_bar+1))
               if(sum(V)>0){
                   work<-sample(1:N,sum(V),replace=T) 
                   particles_bar[V==1,]<-particles[work,]
               }
           }   
       }else{
           particles<-particles+rmt(N, rep(0,d), target$parameters$h2_t[t-1]*target$parameters$Sigma, df=target$parameters$nu)
           w<-observations[t_0,1]*pnorm(particles%*%as.matrix(observations[t_0,2:data_col]), log=TRUE)+
               (1-observations[t_0,1])*pnorm(-particles%*%as.matrix(observations[t_0,2:data_col]), log=TRUE)
       
       }
      
       
  }
  return(list(MEAN=mean_vec[1:count,],MEAN2=mean2_vec[1:count,], BAR=particles_bar, TILDE= particles , ESS=ESS_vec[1:count]))


}
























