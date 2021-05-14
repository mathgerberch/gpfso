









online_time<-function(tuning_parameters, time_max, target, observations, cl=1, PI_BAR=FALSE,  t_bar=1, particles_in=0, W_in=0, t_in=0){
  #//Input parameters
  N<- tuning_parameters$N
  J<- tuning_parameters$J 
  K<- tuning_parameters$K
  ESS_bound<- N*tuning_parameters$ESS_bound         
  acc_rate<-{}           
  T_p <- c(target$parameters$T_p, datalength+1)
  d<-target$dimension
  datalength <- nrow(observations)
  data_col<-ncol(observations)

  #//Initialization
  w<-rep(0,N)
  if(t_in==0){ ##sample from prior
    N<- tuning_parameters$N
    particles_0<-as.matrix(target$rprior(N, target$parameters))
    particles<-particles_0
    w<-rep(0,N)
    W<-rep(1/N,N)
    W_0<-rep(1/N,N)
  }else{
    N<-length(W_in)
    particles_0<-particles_in
    particles<-particles_in
    W<-W_in
    W_0<-W_in
  }
  ESS<-1/sum(W^2)
  lik0<-rep(0,N)
  particles_bar<-matrix(0,N,d)

  mu_vec<-0
  Cov_mat<-0

  proposals<-matrix(0,N,d)
  w1<-rep(0,N)
  lik1<-0
  kt<-t_in
  work<-1:length(T_p)
  p<-min(work[T_p>t_in])

  scale<-(2.38*d^{-0.5})^2

  mean_vec<-matrix(0,datalength,d)
  mean2_vec<-matrix(0,datalength,d)
  ESS_vec<-rep(0,datalength)

  time_start <- as.numeric(Sys.time())
  time_diff<-0
  t<-t_in+1
  while( (t <= datalength) && time_diff< time_max){
    ESS_vec[(t-t_in)]<-ESS
    if(ESS<ESS_bound && t>1){
       A<-tuning_parameters$resampling(W) 
       if(t-1==T_p[p] ||(t-1==t_in && t_in>0) ){  ##h_{t-1}>0
         particles_0<-target$rmutation(as.matrix(particles[A,]), t-1, target$parameters)
         w<-target$loglikelihood(particles_0, matrix(observations[t,], ncol=data_col), target$parameters, cl)
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
         #work<-lower.tri(Cov_mat, diag = TRUE)
         #Cov_mat<-vec2symMat(c(Cov_mat[work==TRUE]), diag = TRUE)
         if(min(eigen(Cov_mat)$values)<=10^{-7})  Cov_mat<-diag(diag(Cov_mat),d)
         proposals<-as.matrix(particles[A,])+rmnorm(N, rep(0,d),Cov_mat)
         lik1<-target$loglikelihood(proposals, matrix(observations[(kt+1):(t-1),], ncol=data_col),target$parameters, cl)
         if(mean(abs(lik0))==0) lik0<-target$loglikelihood(particles, matrix(observations[(kt+1):(t-1),], ncol=data_col), target$parameters, cl)
         if(kt==0){
             alpha<-lik1-lik0[A]+target$dprior(proposals,target$parameters)-target$dprior(as.matrix(particles[A,]),target$parameters)
         }else{
             #if(J>1){
             #   A0<-matrix(sample(1:N,N*J,replace=T, prob=W_0),N,byrow=T)
             #   A1<-matrix(sample(1:N,N*J,replace=T, prob=W_0),N,byrow=T)
             #   work1<-0
             #   work2<-0
             #   for(j in 1:J){
             #     work1<-work1+exp(target$dmutation(proposals-particles_0[A0[,j],], kt,target$parameters))
             #     work2<-work2+exp(target$dmutation(as.matrix(particles[A,])-particles_0[A1[,j],], kt,target$parameters))
             #   }
             #   alpha<-lik1-lik0[A]+log(work1)-log(work2) 
             #   alpha[is.nan(alpha)==TRUE]<-log(0) 
             #}else{ 
                A0<-sample(1:N,N,replace=T, prob=W_0)
                A1<-sample(1:N,N,replace=T, prob=W_0)
                alpha<-lik1-lik0[A]+
                    target$dmutation(proposals-particles_0[A0,], kt,target$parameters)-
                    target$dmutation(as.matrix(particles[A,])-particles_0[A1,], kt,target$parameters)
             #}
          }
          U<-log(runif(N))
          work<- U<=alpha
          acc_rate<-c(acc_rate,sum(work)/N)
          particles<-as.matrix(particles[A,])
          particles[work,]<-proposals[work,] 
          lik0<-lik0[A]
          lik0[work]<-lik1[work]
          w<-target$loglikelihood(particles,  matrix(observations[t,], ncol=data_col),target$parameters, cl)
          w1<- exp(w - max(w))
    	  W<- w1 / sum(w1)
          lik0<-lik0+w
          #cat("time is", t)
       }
       ESS<-1/sum(W^2)
       #cat("ESS equals", ESS)
       
    }else{
       if(t-1==T_p[p]){
          particles_0<-target$rmutation(particles,t-1,target$parameters)
          w<-w+target$loglikelihood(particles_0,  matrix(observations[t,], ncol=data_col),target$parameters, cl)
          w1<- exp(w - max(w))
    	  W_0<- w1 / sum(w1)
          particles<-particles_0
          W<-W_0
          kt<-T_p[p]
          lik0<-rep(0,N)
          p<-p+1
       }else{
          work<-target$loglikelihood(particles, matrix(observations[t,], ncol=data_col),target$parameters, cl)
          lik0<-lik0+work
          w<-w+ work 
          w1<- exp(w - max(w))
    	  W<- w1 / sum(w1)
       }
       ESS<-1/sum(W^2)
    }
    mean_vec[(t-t_in),]<-apply(W*particles,2,sum)
    mean2_vec[(t-t_in),]<-apply(W*particles^2,2,sum)
    #ESS<-1/sum(W^2)
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
    t<-t+1
    time_diff<-as.numeric(Sys.time())-time_start
    #cat("ESS equals", ESS)
    #cat("Bias is ", max(abs(mean_vec[t,2:d]-theta_star[2:d])))
  }
  return(list(MEAN=mean_vec[1:(t-1),], MEAN2=mean2_vec[1:(t-1),], BAR=particles_bar, TILDE=cbind(particles,W), ESS=ESS_vec[1:(t-1)],ACC=acc_rate))


}
























