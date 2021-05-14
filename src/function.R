	
##################################################
# SSP RESAMPLING
##################################################
SSP_Resampling<-function(weights){
	N<-length(weights)
	return(1+.C('SSPResampler_R', as.double(runif(N)), as.integer(N), as.double(weights), res=double(N))$res)
}
##################################################
# ADAGRAD FOR CQR MODEL
##################################################

ADA_cqr<-function(observations, quantile, c, theta0, epsilon=0.000001){
        d<-ncol(observations)-1
        T_end<-nrow(observations)
        if(length(theta0)==d){
           M<-1
        }else{
           M<-nrow(theta0)
        }
        val<-.C("ADA_censored", as.double(quantile), as.double(epsilon), as.double(c(t(theta0))), as.integer(M), 
                 as.integer(d), as.integer(T_end), as.double(c), as.double(c(t(observations))), res=double(2*M*d))
        est_mat<-matrix(val$res, ncol=2*d, nrow=M, byrow=TRUE)
        return(list(TILDE=est_mat[,1:d], BAR=est_mat[,(d+1):(2*d)]))
}
##################################################
# GPFSO FOR CQR MODEL
##################################################

GPFSO_cqr<-function(tuning_parameters, target , Data, cl=1, seed){
   N<-tuning_parameters$N
   target_parameters<-target$parameters
   c_ess<-tuning_parameters$c_ess
   d<-target$dimension
   T_end<-nrow(Data)
   val<-.C('censored_online', as.integer(cl), as.double(target_parameters$quantile), as.double(c_ess*N), as.double(target_parameters$prior_mean), 
         as.double(sqrt(target_parameters$prior_var)), as.integer(seed),  as.integer(N),  as.integer(d), as.double(target_parameters$nu), as.integer(nrow(Data)), 
                   as.double(target_parameters$learning_rate), as.integer(target_parameters$student), as.double(c(t(Data))), res=double(d*nrow(Data)))
    theta_mat<-matrix(val$res,ncol=d,byrow=TRUE)
    return(list(MEAN=theta_mat))
}
##################################################
# JITTERING FOR CQR MODEL
##################################################
jittering_cqr<-function(tuning_parameters, target, iota, Data, cl=1, seed){
   N<-tuning_parameters$N
   target_parameters<-target$parameters
   c_ess<-tuning_parameters$ESS_bound
   d<-target$dimension
   val<-.C('censored_online_jitter', as.integer(cl), as.double(target_parameters$quantile), as.double(c_ess*N), as.double(target_parameters$prior_mean), 
          as.double(sqrt(target_parameters$prior_var)), as.integer(seed), as.integer(N),  as.integer(d), as.integer(nrow(Data)), 
                   as.double(iota), as.double(c(t(Data))), res=double(d*nrow(Data)))
    theta_mat<-matrix(val$res,ncol=d,byrow=TRUE)
    return(list(MEAN=theta_mat))
}

##################################################
# GPFSO FOR toy multimodal MODEL
##################################################

GPFSO_multi<-function(tuning_parameters, target , Data, cl=1, seed){
   N<-tuning_parameters$N
   target_parameters<-target$parameters
   c_ess<-tuning_parameters$ESS_bound
   d<-target$dimension
   T_end<-nrow(Data)
   val<-.C('GPFSO_multimodal', as.integer(cl),  as.double(c_ess*N),  as.integer(seed), 
         as.integer(N),  as.integer(d), as.double(target_parameters$nu), as.integer(nrow(Data)), 
                   as.double(target_parameters$learning_rate), as.integer(target_parameters$student), as.double(c(t(Data))), res=double(d*nrow(Data)))
    theta_mat<-matrix(val$res,ncol=d,byrow=TRUE)
    return(list(MEAN=theta_mat))
}





loglikelihood_multimodal<-function(theta, Data, parameters, cl){
	val<-.C('Lik_multimodal', as.integer(cl), as.integer(nrow(theta)), as.integer(ncol(theta)), as.integer(nrow(Data)),  
                  as.double(c(t(theta))), as.double(c(t(Data))), res=double(nrow(theta)) )
        return( val$res)
}

##################################################
# LIKELIHOOD FUNCTION FOR G-AND-K DISTRIBUTION
##################################################



fast_gandk<- function(theta, y, param=NULL, cl=1, maxsteps = 1000,  tolerance = 1e-10, lower = 1e-20, upper = 1-1e-20){
     val<-.C('Lik_gandk', as.integer(cl), as.integer(nrow(theta)),  as.double(c(t(theta))), as.double(y), 
        as.integer(maxsteps), as.double(tolerance), as.double(lower), as.double(upper), res=double(nrow(theta)))$res
     return(val)
}



fast_gandk_multi_obs<- function(theta, y, cl=1, maxsteps = 1000,  tolerance = 1e-10, lower = 1e-20, upper = 1-1e-20){
     if(length(c(y))==2){
        T_end<1
     }else{
        T_end<-nrow(y)
     }
     val<-.C('Lik_gandk_multi_obs', as.integer(cl),   as.double(c(t(theta))), as.integer(T_end), as.double(c(t(y))), 
        as.integer(maxsteps), as.double(tolerance), as.double(lower), as.double(upper), res=double(1))$res
     return(val)
}


fast_gandk_multi_obsN<- function(theta, y, parameters=NULL, cl=1, maxsteps = 1000,  tolerance = 1e-10, lower = 1e-20, upper = 1-1e-20){
     if(length(y)==2){
        Data<-matrix(y,nrow=1,ncol=2)
     }else{
        Data<-y
     }
     T_end<-nrow(Data)
     N<-nrow(theta)
     val<-.C('Lik_gandk_multi_obs_N', as.integer(cl),  as.integer(nrow(theta)), as.double(c(t(theta))), as.integer(nrow(Data)), 
            as.double(c(t(y))),  as.integer(maxsteps), as.double(tolerance), as.double(lower), as.double(upper), res=double(N))$res
     return(val)
}




##################################################
# GPFSO_ FOR SAGM MODEL
##################################################
GPFSO_mixture<-function(N, target , Data, cl=1, seed=-1, c_ess=0.7){
    T_end<-nrow(Data)
    if(seed== -1){
        seed_use<-sample(1:10^5,1)
    }else{
        seed_use=seed
    }
    d=target$d
    val<-.C('GPFSO_mixture_C',  as.integer(cl), as.integer(ncol(Data)), as.double(c_ess*N), as.double(target$parameters$mu), as.double(target$parameters$sd), 
          as.integer(target$parameters$prob),as.integer(target$parameters$mu1), as.integer(target$parameters$mu2), as.integer(target$parameters$sig1), 
           as.integer(target$parameters$sig2), as.integer(c(length(target$parameters$prob), length(target$parameters$mu1),length(target$parameters$mu2),
              length(target$parameters$sig1),length(target$parameters$sig2))), as.integer(seed_use), as.integer(N), as.integer(d), as.double(target$parameters$nu), as.integer(nrow(Data)), 
                   as.double(target$parameters$learning_rate), as.integer(target$parameters$student), as.double(c(t(Data))),
                   res=double(d*nrow(Data)))
   theta_mat<-matrix(val$res,ncol=d,byrow=TRUE)
    return(list(MEAN=theta_mat))
}


##################################################
# JITTERING FOR SAGM MODEL
##################################################

jittering_mixture<-function(N, target, Data, iota, cl=1, seed=-1, c_ess=0.7){
    T_end<-nrow(Data)
    if(seed== -1){
        seed_use<-sample(1:10^5,1)
    }else{
        seed_use=seed
    }
    d=target$d
    val<-.C('PSMCO_mixture_C',  as.integer(cl), as.integer(ncol(Data)), as.double(c_ess*N), as.double(target$parameters$mu), as.double(target$parameters$sd), 
          as.integer(target$parameters$prob),as.integer(target$parameters$mu1), as.integer(target$parameters$mu2), as.integer(target$parameters$sig1), 
           as.integer(target$parameters$sig2), as.integer(c(length(target$parameters$prob), length(target$parameters$mu1),length(target$parameters$mu2),
              length(target$parameters$sig1),length(target$parameters$sig2))), as.integer(seed_use), as.integer(N), as.integer(d), as.double(target$parameters$nu), as.integer(nrow(Data)), 
                   as.double(iota), as.integer(target$parameters$student), as.double(c(t(Data))),
                   res=double(d*nrow(Data)))
   theta_mat<-matrix(val$res,ncol=d,byrow=TRUE)
    return(list(MEAN=theta_mat))
}
##################################################################

 

















