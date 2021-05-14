









online<-function(tuning_parameters, target, observations, cl=1, PI_BAR=FALSE,  t_bar=1, particles_in=0, W_in=0, t_in=0){
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

  for(t in (t_in+1):datalength){
    ESS_vec[(t-t_in)]<-ESS
    if(ESS<ESS_bound && t>1){
       A<-tuning_parameters$resampling(W) 
       if(t-1==T_p[p] ||(t-1==t_in && t_in>0) ){  ##h_{t-1}>0
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
         #work<-lower.tri(Cov_mat, diag = TRUE)
         #Cov_mat<-vec2symMat(c(Cov_mat[work==TRUE]), diag = TRUE)
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
    mean_vec[(t-t_in),]<-apply(W*particles,2,sum)
    mean2_vec[(t-t_in),]<-apply(W*particles^2,2,sum)
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
  }
  return(list(MEAN=mean_vec,MEAN2=mean2_vec, BAR=particles_bar, TILDE=cbind(particles,W), ESS=ESS_vec, ACC=acc_rate))


}

######################################################



online_Sigma<-function(tuning_parameters, target_in, observations, cl=1, PI_BAR=FALSE, Sig_par=c(1,0.01),  t_bar=1, particles_in=0, W_in=0, t_in=0){
  #//Input parameters
  target<-target_in
  N<- tuning_parameters$N
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
  Sigma0<-target$target_parameters$Sigma
  Sigma1<-target$target_parameters$Sigma
  Sigma2<-target$target_parameters$Sigma
  t_count0<-Sig_par[1]

  for(t in (t_in+1):datalength){
    ESS_vec[(t-t_in)]<-ESS
    if(ESS<ESS_bound && t>1){
       A<-tuning_parameters$resampling(W) 
       if(t-1==T_p[p] ||(t-1==t_in && t_in>0) ){  ##h_{t-1}>0
         if(t>T_p[Sig_par[1]]){
            target$target_parameters$Sigma<-Sigma2
         }
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
         Sigma2<-Sigma1
         Sigma1<-Sigma0
         Sigma0<-cov.wt(particles, W, cor=TRUE)$cor 
         Sigma0[lower.tri(Sigma0)] = t(Sigma0)[lower.tri(Sigma0)]
       }else{ ##h_{t-1}=0 
         Cov_mat<-scale*cov.wt(particles, W)$cov 
         Cov_mat[lower.tri(Cov_mat)] = t(Cov_mat)[lower.tri(Cov_mat)]
         #work<-lower.tri(Cov_mat, diag = TRUE)
         #Cov_mat<-vec2symMat(c(Cov_mat[work==TRUE]), diag = TRUE)
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
                A0<-sample(1:N,N,replace=T, prob=W_0)
                A1<-sample(1:N,N,replace=T, prob=W_0)
                alpha<-lik1-lik0[A]+
                    target$dmutation(proposals-particles_0[A0,], kt, target$parameters)-
                    target$dmutation(as.matrix(particles[A,])-particles_0[A1,], kt, target$parameters)
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
   
     


    mean_vec[(t-t_in),]<-apply(W*particles,2,sum)
    mean2_vec[(t-t_in),]<-apply(W*particles^2,2,sum)
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
  }
  return(list(MEAN=mean_vec,MEAN2=mean2_vec, BAR=particles_bar, TILDE=cbind(particles,W), ESS=ESS_vec, ACC=acc_rate))


}

######################################################




online_qmc<-function(tuning_parameters, target, observations, cl=1, PI_BAR=FALSE){
  #//Input parameters
  N<- tuning_parameters$N
  K<- tuning_parameters$K
  d<-target$dimension
  datalength <- nrow(observations)
  data_col<-ncol(observations)

  #//Initialization
  particles_bar<-matrix(0,N,d)
  mean_vec<-matrix(0,datalength,d)
  mean2_vec<-matrix(0,datalength,d)
  ESS_vec<-rep(0,datalength)
  
  #sample prior 
  points<-Sobol(N, d , scrambling=TRUE, seed=sample( 1:100000,1) )
  particles<-rmnorm_qmc(points, d, target$parameters$prior_mean,target$parameters$prior_covariance) 
  
  #time 1:
  w<-target$loglikelihood(particles,  matrix(observations[1,], ncol=data_col),target$parameters, cl)
  w[is.nan(w)]<--Inf
  w1<- exp(w - max(w))
  W<- w1 / sum(w1)
  ESS_vec[1]<-1/sum(W^2)
  particles_bar<-particles[sample(1:N,K,replace=T, prob=W),]

  for(t in 2:datalength){
       points<-Sobol(N,d+1 , scrambling=TRUE)	
       A<-Hilbert_Resampling(points[,1], particles, W)
       if(target$parameters$student[t-1]==0){
            particles<-as.matrix(particles[A,])+rmnorm_qmc(as.matrix(points[,2:(d+1)]), d, rep(0,d), target$parameters$h2_t[t-1]*target$parameters$Sigma)
       }else{
            particles<-as.matrix(particles[A,])+rmvt_qmc(as.matrix(points[,2:(d+1)]), d, rep(0,d), target$parameters$h2_t[t-1]*target$parameters$Sigma, target$parameters$nu) 
       }
       w<-target$loglikelihood(particles, matrix(observations[t,], ncol=data_col), target$parameters, cl)
       w[is.nan(w)]<--Inf
       w1<- exp(w - max(w))
       W<- w1 / sum(w1)
       ESS_vec[t]<-1/sum(W^2)
       mean_vec[t,]<-apply(W*particles,2,sum)
       mean2_vec[t,]<-apply(W*particles^2,2,sum)
       
       if(PI_BAR){
          V<-rbinom(K, 1, prob=1/t)
          if(sum(V)>0){
            work<-sample(1:N,sum(V),replace=T, prob=W) 
            particles_bar[V==1,]<-particles[work,]
          }
       }   
  }
  free_Sobol()
  return(list(MEAN=mean_vec,MEAN2=mean2_vec, BAR=particles_bar, TILDE=cbind(particles,W), ESS=ESS_vec))


}



##############################################





mmd2<-function(theta,theta_post,  gamma){
    Kval<-length(theta$BAR[,1]) 
    d<-ncol(theta$BAR)

    w1<-c(outer(theta$BAR[,1],theta$BAR[,1], FUN='-'))^2
    for(i in 2:d){
        w1<-w1+c(outer(theta$BAR[,i],theta$BAR[,i], FUN='-'))^2
    }
    w1<-sqrt(w1)

    part1<-mean(exp(-w1/gamma))
    
    w1<-c(outer(theta_post$TILDE[,1],theta_post$TILDE[,1], FUN='-'))^2
    for(i in 2:d){
       w1<-w1+c(outer(theta_post$TILDE[,i],theta_post$TILDE[,i], FUN='-'))^2
    }
    w1<-sqrt(w1)
    weight1<-c(outer(theta_post$TILDE[,d+1],theta_post$TILDE[, d+1], FUN='*'))
    part2<-sum(exp(-w1/gamma)*weight1)

    w1<-c(outer(theta$BAR[,1],theta_post$TILDE[,1], FUN='-'))^2
    for(i in 2:d){
       w1<-w1+c(outer(theta$BAR[,i],theta_post$TILDE[,i], FUN='-'))^2
    }
    w1<-sqrt(w1)
    weight1<-c(outer(rep(1/Kval,Kval),theta_post$TILDE[, d+1], FUN='*'))
    part3<-sum(exp(-w1/gamma)*weight1)
    return(part1+part2-2*part3)
}



mmd2_bis<-function(theta,theta_post,  gamma){
    Kval<-length(theta$BAR[,1]) 
    d<-ncol(theta$BAR)

    w1<-c(outer(theta$TILDE[,1],theta$TILDE[,1], FUN='-'))^2
    for(i in 2:d){
       w1<-w1+c(outer(theta$TILDE[,i],theta$TILDE[,i], FUN='-'))^2
    }
    w1<-sqrt(w1)
    weight1<-c(outer(theta$TILDE[,d+1],theta$TILDE[, d+1], FUN='*'))
    part1<-mean(exp(-w1/gamma))
    
    w1<-c(outer(theta_post$TILDE[,1],theta_post$TILDE[,1], FUN='-'))^2
    for(i in 2:d){
       w1<-w1+c(outer(theta_post$TILDE[,i],theta_post$TILDE[,i], FUN='-'))^2
    }
    w1<-sqrt(w1)
    weight1<-c(outer(theta_post$TILDE[,d+1],theta_post$TILDE[, d+1], FUN='*'))
    part2<-sum(exp(-w1/gamma)*weight1)

    w1<-c(outer(theta$TILDE[,1],theta_post$TILDE[,1], FUN='-'))^2
    for(i in 2:d){
       w1<-w1+c(outer(theta$TILDE[,i],theta_post$TILDE[,i], FUN='-'))^2
    }
    w1<-sqrt(w1)
    weight1<-c(outer(theta$TILDE[,d+1],theta_post$TILDE[, d+1], FUN='*'))
    part3<-sum(exp(-w1/gamma)*weight1)
    return(part1+part2-2*part3)
}


##############################################

compute_param<-function(alpha, Delta, datalength, c_h=1){
    h_seq<-c(c_h*(1:datalength)^(-alpha))
    scale1<-10
    diff_min<-10
    f_m<-function(x){return (max(0.5,x^{alpha*0.8}))}
    stud_YN<-rep(0,datalength)
    t_0<-10
    t_1<-1
    while(t_1 <= datalength){
        t_1<-t_0+max(diff_min,scale1*floor(f_m(t_0)*log(1/h_seq[t_0])))
        if(t_1 <= datalength) stud_YN[t_1]<-1
        t_0<-t_1
   }
   if(Delta==1){
     start<-1
     pert_times<-start:datalength
     stud2_YN<-stud_YN
   }else{
     gap<-Delta
     start<-Delta
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
  }
  return(list(H=h_seq, T_seq=pert_times, TP=stud2_YN))

}

















