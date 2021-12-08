

#include <stdio.h>
#include <string.h>
#include <omp.h>
#include "../include/functions.h"	





void SSPResampler(double* points, int* N,  double *W, int *J)
{
	int i,j,k, index1, index2;
	double s,p1,p2;

	j=0;
	for(i=0;i<*N;i++)
	{
	       for(k=0;k<floor(*N*W[i]); k++)
	       {
               	J[j]=i;
			j++;
		}
	}

        if(j<*N){
	   s=0;
	   index1=0;
	   index2=1;
	   p1=*N*W[index1]-floor(*N*W[index1]);
	   for(i=0;i<*N;i++)
	   {
	        p2=*N*W[index2]-floor(*N*W[index2]);
		s=p1+p2;
		if(s<1.0)
		{
			if(points[i]<p1/s)
			{
				p1=s;
				index2++;
			}
			else
			{
			        index1=index2;
				 p1=s;
				 index2++;
			}
		}
		else if(s==1)
		{
			if(points[i]<p1)
			{
				J[j]=index1;
				index1=index2+1;
				index2=index2+2;
				p1=*N*W[index1]-floor(*N*W[index1]);
				j++;
			}
			else
			{
				J[j]=index2;
				index1=index2+1;
				index2=index2+2;
				p1=*N*W[index1]-floor(*N*W[index1]);
				j++;
			}
		}
		else
		{
			if(points[i]<(1.0-p1)/(2.0-s))
			{
				p1=s-1.0;
				J[j]=index2;
				index2++;
				j++;
			}
			else
			{
			        J[j]=index1;
			        index1=index2;
				p1=s-1.0;
				index2++;
				j++;
			}
		}
		if(index2==*N)
		{
			break;
		}
	   }
	   if(j==*N-1)
	   {
		  J[*N-1]=index1;
	   }
	}
}


void SSPResampler_R(double* points, int* N, double *W, double *J)
{
	int i,j,k, index1, index2;
	double s,p1,p2;

	j=0;
	for(i=0;i<*N;i++)
	{
	       for(k=0;k<floor(*N*W[i]); k++)
	       {
               	J[j]=i;
			j++;
		}
	}

        if(j<*N){
	   s=0;
	   index1=0;
	   index2=1;
	   p1=*N*W[index1]-floor(*N*W[index1]);
	   for(i=0;i<*N;i++)
	   {
	        p2=*N*W[index2]-floor(*N*W[index2]);
		s=p1+p2;
		if(s<1.0)
		{
			if(points[i]<p1/s)
			{
				p1=s;
				index2++;
			}
			else
			{
			        index1=index2;
				 p1=s;
				 index2++;
			}
		}
		else if(s==1)
		{
			if(points[i]<p1)
			{
				J[j]=index1;
				index1=index2+1;
				index2=index2+2;
				p1=*N*W[index1]-floor(*N*W[index1]);
				j++;
			}
			else
			{
				J[j]=index2;
				index1=index2+1;
				index2=index2+2;
				p1=*N*W[index1]-floor(*N*W[index1]);
				j++;
			}
		}
		else
		{
			if(points[i]<(1.0-p1)/(2.0-s))
			{
				p1=s-1.0;
				J[j]=index2;
				index2++;
				j++;
			}
			else
			{
			        J[j]=index1;
			        index1=index2;
				p1=s-1.0;
				index2++;
				j++;
			}
		}
		if(index2==*N)
		{
			break;
		}
	   }
	   if(j==*N-1)
	   {
		  J[*N-1]=index1;
	   }
	}
}


//////////////////////////////////////////////////////////////////////////////////
                   
                   
void censored_online(int *cl, double *q, double *ess_bound, double *mu, double *sigma, int *seed, int *N, int *d, double *nu, int *datalength, 
                   double *alpha, int *tp,  double *observations, double *res){

    omp_set_dynamic(0);     		
    omp_set_num_threads(*cl);

    int t, i, k,  T_max;
    double ess_val, h_t, work;
    double sqrt_nu=pow(*nu,0.5);
  
    T_max=*datalength;
    ess_val=*N;


    double *theta=(double*)malloc(sizeof(double)*(*d*(*N)));
    double *thetah=(double*)malloc(sizeof(double)*(*d*(*N)));
    double *w=(double*)malloc(sizeof(double)*(*N));
    double *W=(double*)malloc(sizeof(double)*(*N));
    int *J=(int*)malloc(sizeof(int)*(*N));

    //MTRand r = seedRand(*seed);

    gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set(r, *seed);

    //sample from tilde{pi}_0
    for(i=0; i<*N; i++){
          for(k=0; k<*d; k++){
               theta[*d*i+k]=mu[k]+ sigma[k]*gsl_ran_gaussian(r,1.0);
          }
    }
    //Process observations 1
    t=0;
    //2.1 Compute the weights
    #pragma omp parallel for private(work,k)
    for(i=0; i<*N; i++){
        work=0;
        for(k=0; k<*d; k++){
             work+= observations[(*d+1)*t+(k+1)]*theta[*d*i+k];
        }
        if(work<=0){
            work=0;
        }
        w[i]=-(fabs(observations[(*d+1)*t]-work)+(*q*2.0-1.0)*(observations[(*d+1)*t]-work))*0.5;   
    }
    work=weight(w,  W, *N);
    //Compute the ESS
    ess_val=0;
    for(i=0;i<*N;i++){
            ess_val+=W[i]*W[i];
    }
    ess_val=1.0/ess_val;

    //Compute the theta_tilde
    for(k=0; k<*d;k++){
         res[*d*t+k]=0;
         for(i=0;i<*N;i++){
              res[*d*t+k]+=W[i]*theta[*d*i+k];
         }
    }
   
    //Process other observations
    for(t=1; t<T_max; t++)
    {
             h_t=alpha[0]*pow(t,-alpha[1]);  //compute h_t
             if(ess_val> *ess_bound){ //no resampling
               if(tp[t-1]==1){  //Student move
                 for(i=0; i<*N; i++){
                       work=pow(gsl_ran_chisq(r, *nu),-0.5);
                       for(k=0; k<*d;k++){
                           theta[*d*i+k]+=h_t*gsl_ran_gaussian(r,1.0)*sqrt_nu*work;
                       }
                 }
              }
              else{ //Gaussian move
                 for(i=0; i<*N; i++){
                       for(k=0; k<*d;k++){
                           theta[*d*i+k]+=h_t*gsl_ran_gaussian(r,1.0);
                       }
                 }
              }
              #pragma omp parallel for private(work,k)
              for(i=0; i<*N; i++){
                 work=0; 
                 for(k=0; k<*d; k++){
                      work+= observations[(*d+1)*t+(k+1)]*theta[*d*i+k];
                 }
                 if(work<=0){
                     work=0;
                 }
                 w[i]+=-(fabs(observations[(*d+1)*t]-work)+(*q*2.0-1.0)*(observations[(*d+1)*t]-work))*0.5;      
              }
              work=weight(w,  W, *N);
           }
           else{//resampling
                for(i=0; i<*N; i++){
                      w[i]=gsl_rng_uniform(r);
                }
                
                SSPResampler(w, N, W, J);

               // #pragma omp parallel for private(k)
                for(i=0; i<*N; i++){
                      for(k=0; k<*d;k++){
                          thetah[*d*i+k]=theta[*d*J[i]+k];
                      }
                }
                  
                if(tp[t-1]==1){
                     for(i=0; i<*N; i++){
                         work=pow(gsl_ran_chisq(r, *nu),-0.5);
                         for(k=0; k<*d;k++){
                             theta[*d*i+k]=thetah[*d*i+k]+h_t*gsl_ran_gaussian(r,1.0)*sqrt_nu*work;
                         }
                     }
                }
                else{
                    for(i=0; i<*N; i++){
                        for(k=0; k<*d;k++){
                           theta[*d*i+k]=thetah[*d*i+k]+h_t*gsl_ran_gaussian(r,1.0);
                        }
                    }
                }
                #pragma omp parallel for private(work,k)
                for(i=0; i<*N; i++){
                  work=0; 
                 for(k=0; k<*d; k++){
                      work+= observations[(*d+1)*t+(k+1)]*theta[*d*i+k];
                 }
                 if(work<=0){
                     work=0;
                 }
                 w[i]=-(fabs(observations[(*d+1)*t]-work)+(*q*2.0-1.0)*(observations[(*d+1)*t]-work))*0.5;   
               }
              work=weight(w,  W, *N);
            }

            ess_val=0;
            for(i=0;i<*N;i++){
                 ess_val+=W[i]*W[i];
            }
            ess_val=1.0/ess_val;

            for(k=0; k<*d;k++){
                  res[*d*t+k]=0;
                   for(i=0;i<*N;i++){
                       res[*d*t+k]+=W[i]*theta[*d*i+k];
                   }
            }                   
   }
    free(theta);
    theta=NULL;
    free(thetah);
    thetah=NULL;
    free(w);
    w=NULL;
    free(W);
    W=NULL;
    free(J);
    J=NULL;
    gsl_rng_free(r);

}

//////////////////////////////////////////////////////////////

                 
void ADA_censored(double *quantile, double *epsilon, double *theta_in, int *M, int *d, int *datalength,  double *alpha,   double *observations, double *res)
{
    int t,k, m;
    double grad, work1, work, step;
    double *grad_vec=(double*)malloc(sizeof(double)*(*d));

     
    for(m=0;m<*M;m++){
    
       for(k=0; k<*d;k++){
           res[*d*m*2+*d+k]=0;
       }
       for(k=0; k<*d;k++){
           res[*d*m*2+k]=theta_in[*d*m+k];
       }
       for(k=0; k<*d;k++){
          grad_vec[k]=0;
       }

       for(t=0; t< *datalength;t++){
           work1=0; 
           for(k=0; k<*d; k++){
	          work1+= observations[(*d+1)*t+(k+1)]*res[*d*m*2+k];
           }
           if(work1>0){
              work= observations[(*d+1)*t]-work1;
              if(work>0){
                 for(k=0;k<*d;k++){
                    grad= -*quantile*observations[(*d+1)*t+(k+1)];
                    grad_vec[k]+= pow(grad,2.0);
                    res[*d*m*2+k]+= -*alpha*grad/pow(*epsilon+grad_vec[k],0.5);
                 }
             }
             else
             {
                for(k=0;k<*d;k++){
                    grad=(1.0-*quantile)*observations[(*d+1)*t+(k+1)];
                    grad_vec[k]+=pow(grad,2.0);
                    res[*d*m*2+k]+=-*alpha*grad/pow(*epsilon+grad_vec[k],0.5);
                }
             }
          }
          for(k=0;k<*d;k++){
                 res[*d*m*2+*d+k]=(t/(t+1.0))*res[*d*m*2+*d+k]+res[*d*m*2+k]/(t+1.0);
          }
     }
    }
    free(grad_vec);
    grad_vec=NULL;
}

//////////////////////////////////////////////////////////////



void censored_online_jitter(int *cl, double *q, double *ess_bound, double *mu, double *sigma, int *seed, int *N, int *d, int *datalength, 
                   double *alpha,  double *observations, double *res){

    omp_set_dynamic(0);     		
    omp_set_num_threads(*cl);

    int t, i, k;
    double ess_val, h_t, work;
    double proba=1.0/pow(*N,0.5);
    ess_val=*N;


    double *theta=(double*)malloc(sizeof(double)*(*d*(*N)));
    double *thetah=(double*)malloc(sizeof(double)*(*d*(*N)));
    double *w=(double*)malloc(sizeof(double)*(*N));
    double *W=(double*)malloc(sizeof(double)*(*N));
    int *J=(int*)malloc(sizeof(int)*(*N));


    gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set(r, *seed);

    //sample from tilde{pi}_0
    for(i=0; i<*N; i++){
          for(k=0; k<*d; k++){
               theta[*d*i+k]=mu[k]+ sigma[k]*gsl_ran_gaussian(r,1.0);
          }
    }
    //Process observations 1
    t=0;
    //2.1 Compute the weights
    #pragma omp parallel for private(work,k)
    for(i=0; i<*N; i++){
        work=0;
        for(k=0; k<*d; k++){
             work+= observations[(*d+1)*t+(k+1)]*theta[*d*i+k];
        }
        if(work<=0){
            work=0;
        }
        w[i]=-(fabs(observations[(*d+1)*t]-work)+(*q*2.0-1.0)*(observations[(*d+1)*t]-work))*0.5;   
    }
    work=weight(w,  W, *N);
    //Compute the ESS
    ess_val=0;
    for(i=0;i<*N;i++){
            ess_val+=W[i]*W[i];
    }
    ess_val=1.0/ess_val;

    //Compute the theta_tilde
    for(k=0; k<*d;k++){
         res[*d*t+k]=0;
         for(i=0;i<*N;i++){
              res[*d*t+k]+=W[i]*theta[*d*i+k];
         }
    }
   
    //Process other observations
    for(t=1; t<*datalength; t++)
    {
             if(ess_val> *ess_bound){ //no resampling
                 for(i=0; i<*N; i++){
                     if(gsl_rng_uniform(r)<proba){
                       for(k=0; k<*d;k++){
                           theta[*d*i+k]+=alpha[0]*gsl_ran_gaussian(r,1.0);
                       }
                    }
                 }       
                 #pragma omp parallel for private(work,k)
                 for(i=0; i<*N; i++){
                   work=0; 
                   for(k=0; k<*d; k++){
                      work+= observations[(*d+1)*t+(k+1)]*theta[*d*i+k];
                   }
                   if(work<=0){
                      work=0;
                   }
                   w[i]+=-(fabs(observations[(*d+1)*t]-work)+(*q*2.0-1.0)*(observations[(*d+1)*t]-work))*0.5;      
                 }
                 work=weight(w,  W, *N);
           }
           else{//resampling
                for(i=0; i<*N; i++){
                      w[i]=gsl_rng_uniform(r);
                }
                SSPResampler(w, N, W, J);
                for(i=0; i<*N; i++){
                      for(k=0; k<*d;k++){
                          thetah[*d*i+k]=theta[*d*J[i]+k];
                      }
                }      
                for(i=0; i<*N; i++){
                     if(gsl_rng_uniform(r)<proba){
                        for(k=0; k<*d;k++){
                           theta[*d*i+k]=thetah[*d*i+k]+alpha[0]*gsl_ran_gaussian(r,1.0);
                        }
                    }
                    else{
                        for(k=0; k<*d;k++){
                           theta[*d*i+k]=thetah[*d*i+k];
                        }   
                   }
                }
                #pragma omp parallel for private(work,k)
                for(i=0; i<*N; i++){
                  work=0; 
                 for(k=0; k<*d; k++){
                      work+= observations[(*d+1)*t+(k+1)]*theta[*d*i+k];
                 }
                 if(work<=0){
                     work=0;
                 }
                 w[i]=-(fabs(observations[(*d+1)*t]-work)+(*q*2.0-1.0)*(observations[(*d+1)*t]-work))*0.5;   
               }
              work=weight(w,  W, *N);
            }

            ess_val=0;
            for(i=0;i<*N;i++){
                 ess_val+=W[i]*W[i];
            }
            ess_val=1.0/ess_val;

            for(k=0; k<*d;k++){
                  res[*d*t+k]=0;
                   for(i=0;i<*N;i++){
                       res[*d*t+k]+=W[i]*theta[*d*i+k];
                   }
            }                   
   }
    free(theta);
    theta=NULL;
    free(thetah);
    thetah=NULL;
    free(w);
    w=NULL;
    free(W);
    W=NULL;
    free(J);
    J=NULL;
    gsl_rng_free(r);

}

/////////////////////////////////////////////////////////////////////////////////


void GPFSO_multimodal(int *cl,  double *ess_bound, int *seed, int *N, int *d, double *nu, int *datalength, 
                   double *alpha, int *tp,  double *observations, double *res){

    omp_set_dynamic(0);     		
    omp_set_num_threads(*cl);

    int t, i, k,  T_max;
    double ess_val, h_t, work;
    double sqrt_nu=pow(*nu,0.5);
   
    double bound1=-21, bound2=19;
    T_max=*datalength; 
    ess_val=*N;


    double *theta=(double*)malloc(sizeof(double)*(*d*(*N)));
    double *thetah=(double*)malloc(sizeof(double)*(*d*(*N)));
    double *w=(double*)malloc(sizeof(double)*(*N));
    double *W=(double*)malloc(sizeof(double)*(*N));
    int *J=(int*)malloc(sizeof(int)*(*N));

   
    gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set(r, *seed);

    //sample from tilde{pi}_0
    for(i=0; i<*N; i++){
          for(k=0; k<*d; k++){
               theta[*d*i+k]=bound1+  (bound2-bound1)*gsl_rng_uniform(r);
          }
    }
    //Process observations 1
    t=0;
    //2.1 Compute the weights
    #pragma omp parallel for private(work,k)
    for(i=0; i<*N; i++){
        work=0;
        for(k=0;k<*d;k++){
            if(theta[*d*i+k]<bound1 || theta[*d*i+k]>bound2){
               work=1.0;
            }
        }
        if(work==1.0){
           w[i]= -INFINITY;
        }
        else{
          work=observations[(*d+1)*t];
          for(k=0; k<*d; k++){
	        work+= -exp(-observations[(*d+1)*t+(k+1)]*pow(theta[*d*i+k], 2.0))-observations[(*d+1)*t+(k+1)]*theta[*d*i+*d-k-1];
          }
         w[i]=-0.5*fabs(work);   
        }
    }
    work=weight(w,  W, *N);
    //Compute the ESS
    ess_val=0;
    for(i=0;i<*N;i++){
            ess_val+=W[i]*W[i];
    }
    ess_val=1.0/ess_val;
    //Compute the theta_tilde
    for(k=0; k<*d;k++){
         res[*d*t+k]=0;
         for(i=0;i<*N;i++){
              res[*d*t+k]+=W[i]*theta[*d*i+k];
         }

    }
   //Process other observations
    for(t=1; t<T_max; t++)
    {
             h_t=alpha[0]*pow(t,-alpha[1]);  //compute h_t
             if(ess_val> *ess_bound){ //no resampling
               if(tp[t-1]==1){  //Student move
                 for(i=0; i<*N; i++){
                       work=pow(gsl_ran_chisq(r, *nu),-0.5);
                       for(k=0; k<*d;k++){
                           theta[*d*i+k]+=h_t*gsl_ran_gaussian(r,1.0)*sqrt_nu*work;
                       }
                 }
              }
              else{ //Gaussian move
                 for(i=0; i<*N; i++){
                       for(k=0; k<*d;k++){
                           theta[*d*i+k]+=h_t*gsl_ran_gaussian(r,1.0);
                       }
                 }
              }
              #pragma omp parallel for private(work,k)
              for(i=0; i<*N; i++){
                 work=0;
                 for(k=0;k<*d;k++){
                    if(theta[*d*i+k]<bound1 || theta[*d*i+k]>bound2){
                       work=1.0;
                    }
                 }
                 if(work==1.0){
                    w[i]= -INFINITY;
                 }
                 else{
                   work=observations[(*d+1)*t];
                   for(k=0; k<*d; k++){
	               work+= -exp(-observations[(*d+1)*t+(k+1)]*pow(theta[*d*i+k], 2.0))-observations[(*d+1)*t+(k+1)]*theta[*d*i+*d-k-1];
                   }
                   w[i]+=-0.5*fabs(work);  
                 }      
              }
              work=weight(w,  W, *N);
           }
           else{//resampling
                for(i=0; i<*N; i++){
                      w[i]=gsl_rng_uniform(r);
                }
                
                SSPResampler(w,  N, W, J);

               // #pragma omp parallel for private(k)
                for(i=0; i<*N; i++){
                      for(k=0; k<*d;k++){
                          thetah[*d*i+k]=theta[*d*J[i]+k];
                      }
                }
                  
                if(tp[t-1]==1){
                     for(i=0; i<*N; i++){
                         work=pow(gsl_ran_chisq(r, *nu),-0.5);
                         for(k=0; k<*d;k++){
                             theta[*d*i+k]=thetah[*d*i+k]+h_t*gsl_ran_gaussian(r,1.0)*sqrt_nu*work;
                         }
                     }
                }
                else{
                    for(i=0; i<*N; i++){
                        for(k=0; k<*d;k++){
                           theta[*d*i+k]=thetah[*d*i+k]+h_t*gsl_ran_gaussian(r,1.0);
                        }
                    }
                }
                #pragma omp parallel for private(work,k)
                for(i=0; i<*N; i++){
                  work=0;
                  for(k=0;k<*d;k++){
                    if(theta[*d*i+k]<bound1 || theta[*d*i+k]>bound2){
                       work=1.0;
                    }
                  }
                  if(work==1.0){
                     w[i]= -INFINITY;
                  }
                else{
                   work=observations[(*d+1)*t];
                   for(k=0; k<*d; k++){
	               work+= -exp(-observations[(*d+1)*t+(k+1)]*pow(theta[*d*i+k], 2.0))-observations[(*d+1)*t+(k+1)]*theta[*d*i+*d-k-1];
                   }
                    w[i]=-0.5*fabs(work);   
                  }  
                }
                work=weight(w,  W, *N);
            }

            ess_val=0;
            for(i=0;i<*N;i++){
                 ess_val+=W[i]*W[i];
            }
            ess_val=1.0/ess_val;

            for(k=0; k<*d;k++){
                  res[*d*t+k]=0;
                   for(i=0;i<*N;i++){
                       res[*d*t+k]+=W[i]*theta[*d*i+k];
                   }
                     
            }                   
   }
    free(theta);
    theta=NULL;
    free(thetah);
    thetah=NULL;
    free(w);
    w=NULL;
    free(W);
    W=NULL;
    free(J);
    J=NULL;
    gsl_rng_free(r);

}



//Normalize the importance weights and return the normaliyzing constant

double weight(double *wei, double *W, int N)
{
	double s,acc,nc;
	int i;
	
	s = maxmin(wei,N,1);

	nc = 0.0;
	
	for(i=0;i<N;i++)
	{
		nc += exponen(wei[i]-s);
	}
	
	nc = s + log(nc); 		   
	
	for(i=0;i<N;i++)
	{
		W[i] = exp(wei[i] - nc);  
	}
	return(nc);

}



//Computation of the max/min of a "a" dimensional vector "vec"

double maxmin(double *vec,int a,int b)
{
	double val;
	int i;

	val = vec[0];

	if(b==0)	//minimum
	{
		for(i=1;i<a;i++)
		{
			if(vec[i]<val)
			{
				val = vec[i];
			}	
		}
	}
	else      //maximum
	{
		for(i=1;i<a;i++)
		{
			if(vec[i]>val)
			{
				val = vec[i];
			}	
		}
	}
	
	return(val);
}



double exponen(double x)
{
	double val;
		
	if(x>302.0)
	{
		val = 1e300;
	}
	else if(x<-302.0)
	{
		val = 0.0;
	}
	else
	{
		val = exp(x);
	}

	return(val);
}











double sum2(double *vect,int *N)
{
	int i;
	double s=0;
	for(i=0;i<*N;i++)
	{
           s=s+pow(vect[i],2.0);
	}
    return(s);
}



void Lik_gandk_multi_obs(int *cl,  double *theta, int *datalength, double *observations, int *maxsteps, double *tolerance, double *lower, double *upper, double *res)
 {
    omp_set_dynamic(0);     		
    omp_set_num_threads(*cl);

    int i, k,t;
    int d=9;
    double work1, work2, theta1[4], theta2[4];
    double Fy1, Fy2, x1,x2; 
    double delta=0.000;
    double h = 1e-05;
    double bound=5.5;
    
    double *W=(double*)malloc(sizeof(double)*(*datalength));
  
//       if((theta[1]< delta) ||  (theta[3]< delta) || (theta[5]< delta) || (theta[7]< delta)|| (theta[8]< -1.0+delta) || (theta[8]> 1.0-delta) ){
if((theta[1]< delta) || (fabs(theta[2])> bound-delta)|| (theta[3]< -0.045-0.01*pow(theta[2], 2.0)) || (theta[5]< delta) ||
                     (fabs(theta[6])> bound-delta)||  (theta[7]< -0.045-0.01*pow(theta[6], 2.0))|| (fabs(theta[8])> 1.0-delta)  ){
             res[0]= -INFINITY;
       }else
       {
          for(k=0;k<4;k++){
               theta1[k]=theta[k];
               theta2[k]=theta[4+k];
          }
          #pragma omp parallel for private(Fy1,Fy2,x1,x2, work1, work2)
          for(t=0; t<*datalength; t++){    
            Fy1=gandkcdf_(observations[2*t], theta1, maxsteps, tolerance, lower,  upper);
            if(Fy1> 1-1e-15){
               Fy1= 1-1e-15;
            }
            x1=gsl_cdf_ugaussian_Pinv(Fy1);
            Fy2=gandkcdf_(observations[2*t+1], theta2, maxsteps, tolerance, lower,  upper);
            if(Fy2> 1-1e-15){
               Fy2= 1-1e-15;
            }
            x2=gsl_cdf_ugaussian_Pinv(Fy2); 
            work1= 1.0/(1.0-pow(theta[8],2.0));
            work2= -theta[8]*work1;
            W[t]= -0.5*(work1*pow(x1,2.0)+work1*pow(x2,2.0)+2.0*x1*x2*work2)+0.5*log(work1);
            W[t]+= 0.5*pow(x1,2.0)+0.5*pow(x2,2.0);
            Fy1=gandkcdf_(observations[2*t]-h, theta1, maxsteps, tolerance, lower,  upper);
            Fy2=gandkcdf_(observations[2*t]+h, theta1, maxsteps, tolerance, &Fy1,  upper);
            work1= (Fy2-Fy1)/(2.0*h);
            if(work1< 1e-300){
               work1= 1e-300;
            }
            W[t]+=log(work1);
            Fy1=gandkcdf_(observations[2*t+1]-h, theta2, maxsteps, tolerance, lower,  upper);
            Fy2=gandkcdf_(observations[2*t+1]+h, theta2, maxsteps, tolerance, &Fy1,  upper);
            work2= (Fy2-Fy1)/(2.0*h);
            if(work2< 1e-300){
               work2= 1e-300;
            }
            W[t]+=log(work2);
          }
        }
        
        res[0]=0;
        for(t=0;t<*datalength; t++){
           res[0]+=W[t];
        }
  
    free(W);
    W=NULL;
}


void Lik_gandk_multi_obs_N(int *cl, int *N,    double *theta, int *datalength, double *observations, int *maxsteps, double *tolerance, double *lower, double *upper, double *res)
 {
    omp_set_dynamic(0);     		
    omp_set_num_threads(*cl);

    int i, k,t;
    int d=9;
    double work1, work2, theta1[4], theta2[4], W;
    double Fy1, Fy2, x1,x2; 
    double delta=0.000;
    double h = 1e-05;
    double bound=5.5;

   #pragma omp parallel for private(k, theta1,theta2, W, Fy1,Fy2,x1,x2, work1, t,work2)
   for(i=0; i<*N; i++){
       //if((theta[d*i+1]< delta) ||  (theta[d*i+3]< delta) || (theta[d*i+5]< delta) || (theta[d*i+7]< delta)|| (theta[d*i+8]< -1.0+delta) || (theta[d*i+8]> 1.0-delta) ){
       if((theta[d*i+1]< delta) || (fabs(theta[d*i+2])> bound-delta)|| (theta[d*i+3]< -0.045-0.01*pow(theta[d*i+2], 2.0)) || (theta[d*i+5]< delta) ||
                     (fabs(theta[d*i+6])> bound-delta)||  (theta[d*i+7]< -0.045-0.01*pow(theta[d*i+6], 2.0))|| (fabs(theta[d*i+8])> 1.0-delta)  ){
             res[i]= -INFINITY;
       }else{

           for(k=0;k<4;k++){
              theta1[k]=theta[d*i+k];
              theta2[k]=theta[d*i+4+k];
           }
           W=0;
           for(t=0; t<*datalength;t++){
             Fy1=gandkcdf_(observations[2*t], theta1, maxsteps, tolerance, lower,  upper);
             if(Fy1> 1-1e-15){
                Fy1= 1-1e-15;
             }
             x1=gsl_cdf_ugaussian_Pinv(Fy1);
             Fy2=gandkcdf_(observations[2*t+1], theta2, maxsteps, tolerance, lower,  upper);
             if(Fy2> 1-1e-15){
                Fy2= 1-1e-15;
             }
             x2=gsl_cdf_ugaussian_Pinv(Fy2); 
             work1= 1.0/(1.0-pow(theta[d*i+8],2.0));
             work2= -theta[d*i+8]*work1;
             W+= -0.5*(work1*pow(x1,2.0)+work1*pow(x2,2.0)+2.0*x1*x2*work2)+0.5*log(work1);
             W+= 0.5*pow(x1,2.0)+0.5*pow(x2,2.0);
             Fy1=gandkcdf_(observations[2*t]-h, theta1, maxsteps, tolerance, lower,  upper);
             Fy2=gandkcdf_(observations[2*t]+h, theta1, maxsteps, tolerance, &Fy1,  upper);
             work1= (Fy2-Fy1)/(2.0*h);
             if(work1< 1e-300){
                work1= 1e-300;
             }
             W+=log(work1);
             Fy1=gandkcdf_(observations[2*t+1]-h, theta2, maxsteps, tolerance, lower,  upper);
             Fy2=gandkcdf_(observations[2*t+1]+h, theta2, maxsteps, tolerance, &Fy1,  upper);
             work2= (Fy2-Fy1)/(2.0*h);
             if(work2< 1e-300){
                work2= 1e-300;
             }
             W+=log(work2);
           }
           res[i]=W;
       }
    }
}
    

  
  
void Lik_gandk(int *cl, int *N,    double *theta, double *observations, int *maxsteps, double *tolerance, double *lower, double *upper, double *res)
 {
    omp_set_dynamic(0);     		
    omp_set_num_threads(*cl);

    int i, k;
    int d=9;
    double work1, work2, theta1[4], theta2[4];
    double Fy1, Fy2, x1,x2; 
    double delta=0.000;
    double h = 1e-05;
    double bound=5.5;
   #pragma omp parallel for private(k, theta1,theta2, Fy1,Fy2,x1,x2, work1, work2)
   for(i=0; i<*N; i++){
//       if((theta[d*i+1]< delta) ||  (theta[d*i+3]< delta) || (theta[d*i+5]< delta) || (theta[d*i+7]< delta)|| (theta[d*i+8]< -1.0+delta) || (theta[d*i+8]> 1.0-delta) ){
if((theta[d*i+1]< delta) || (fabs(theta[d*i+2])> bound-delta)|| (theta[d*i+3]< -0.045-0.01*pow(theta[d*i+2], 2.0)) || (theta[d*i+5]< delta) ||
                     (fabs(theta[d*i+6])> bound-delta)||  (theta[d*i+7]< -0.045-0.01*pow(theta[d*i+6], 2.0))|| (fabs(theta[d*i+8])> 1.0-delta)  ){
             res[i]= -INFINITY;
       }else{
           for(k=0;k<4;k++){
              theta1[k]=theta[d*i+k];
              theta2[k]=theta[d*i+4+k];
           }
           Fy1=gandkcdf_(observations[0], theta1, maxsteps, tolerance, lower,  upper);
           if(Fy1> 1-1e-15){
              Fy1= 1-1e-15;
           }
           x1=gsl_cdf_ugaussian_Pinv(Fy1);
           Fy2=gandkcdf_(observations[1], theta2, maxsteps, tolerance, lower,  upper);
           if(Fy2> 1-1e-15){
              Fy2= 1-1e-15;
           }
           x2=gsl_cdf_ugaussian_Pinv(Fy2); 
           work1= 1.0/(1.0-pow(theta[d*i+8],2.0));
           work2= -theta[d*i+8]*work1;
           res[i]= -0.5*(work1*pow(x1,2.0)+work1*pow(x2,2.0)+2.0*x1*x2*work2)+0.5*log(work1);
           res[i]+= 0.5*pow(x1,2.0)+0.5*pow(x2,2.0);
           Fy1=gandkcdf_(observations[0]-h, theta1, maxsteps, tolerance, lower,  upper);
           Fy2=gandkcdf_(observations[0]+h, theta1, maxsteps, tolerance, &Fy1,  upper);
           work1= (Fy2-Fy1)/(2.0*h);
           if(work1< 1e-300){
              work1= 1e-300;
           }
           res[i]+=log(work1);
           Fy1=gandkcdf_(observations[1]-h, theta2, maxsteps, tolerance, lower,  upper);
           Fy2=gandkcdf_(observations[1]+h, theta2, maxsteps, tolerance, &Fy1,  upper);
           work2= (Fy2-Fy1)/(2.0*h);
           if(work2< 1e-300){
              work2= 1e-300;
           }
           res[i]+=log(work2);
       }
    }
}




  
  
   

double gandkcdf_(double y, double *theta, int *maxsteps, double *tolerance, double *lower, double *upper){
  double A = theta[0];
  double B = theta[1];
  double c = 0.8;
  double g = theta[2];
  double k = theta[3];

  int istep = 0;
  double current_try = 0.5*(*upper + *lower);
  double current_size = 0.25*(*upper - *lower);
  double dd= current_try;
  double z= gsl_cdf_ugaussian_Pinv(dd);
  double fattempt = A + B * (1.0 + c * (1.0 - exp(- g * z)) / (1.0 + exp(- g * z))) * pow((1.0 + pow(z, 2.0)), k) * z;
  while (!(fattempt > y-*tolerance && fattempt < y+*tolerance) && (istep < *maxsteps)){
    istep++;
    if (fattempt > y-*tolerance){
      current_try+=  - current_size;
      dd = current_try;
      z = gsl_cdf_ugaussian_Pinv(dd);
      fattempt = A + B * (1.0 + c * (1.0 - exp(- g * z)) / (1.0 + exp(- g * z))) * pow((1.0 + pow(z, 2.0)), k) * z;
      current_size*=0.5;
    } else {
      current_try+=  current_size;
      dd = current_try;
      z = gsl_cdf_ugaussian_Pinv(dd);
      fattempt = A + B * (1.0 + c * (1.0 - exp(- g * z)) / (1.0 + exp(- g * z))) * pow((1.0 + pow(z, 2.0)), k) * z;
      current_size*= 0.5;
    }
  }
  return(current_try);
}



void gandkcdf2_(double *y, double *theta, int *maxsteps, double *tolerance, double *lower, double *upper, double *res){
  double A = theta[0];
  double B = theta[1];
  double c = 0.8;
  double g = theta[2];
  double k = theta[3];

  int istep = 0;
  double current_try = 0.5*(*upper + *lower);
  double current_size = 0.25*(*upper - *lower);
  double dd= current_try;
  double z= gsl_cdf_ugaussian_Pinv(dd);
  double fattempt = A + B * (1.0 + c * (1.0 - exp(- g * z)) / (1.0 + exp(- g * z))) * pow((1.0 + pow(z, 2.0)), k) * z;
   while (!(fattempt > *y-*tolerance && fattempt < *y+*tolerance) && (istep < *maxsteps)){
    istep++;
    if (fattempt > *y-*tolerance){
      current_try+=  - current_size;
      dd = current_try;
      z = gsl_cdf_ugaussian_Pinv(dd);
      fattempt = A + B * (1.0 + c * (1.0 - exp(- g * z)) / (1.0 + exp(- g * z))) * pow((1.0 + pow(z, 2.0)), k) * z;
      current_size*=0.5;
    } else {
      current_try+=  current_size;
      dd = current_try;
      z = gsl_cdf_ugaussian_Pinv(dd);
      fattempt = A + B * (1.0 + c * (1.0 - exp(- g * z)) / (1.0 + exp(- g * z))) * pow((1.0 + pow(z, 2.0)), k) * z;
      current_size*= 0.5;
    }
  }
  
  res[0]=current_try;
}



void GPFSO_mixture_C(int *cl,  int *dobs, double *ess_bound, double *mu, double *sigma, int *a_vec, int *m1_vec,int *m2_vec,
                     int *sig1_vec, int *sig2_vec, int *vec_length, int *seed, int *N, int *d, double *nu, int *datalength, 
                     double *alpha, int *tp,  double *observations, double *res){

    omp_set_dynamic(0);     		
    omp_set_num_threads(*cl);

    int t, i, k, count;
    double ess_val, h_t, work, prob, m1,m2,sig1,sig2; 
    double sqrt_nu=pow(*nu,0.5);
    


    double *theta=(double*)malloc(sizeof(double)*(*d*(*N)));
    double *thetah=(double*)malloc(sizeof(double)*(*d*(*N)));
    double *w=(double*)malloc(sizeof(double)*(*N));
    double *W=(double*)malloc(sizeof(double)*(*N));
    int *J=(int*)malloc(sizeof(int)*(*N));

   
    gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set(r, *seed);

    //sample from tilde{pi}_0
    for(i=0; i<*N; i++){
          theta[*d*i]=gsl_ran_exponential(r,mu[0]);
          for(k=1; k<*d; k++){
               theta[*d*i+k]=mu[k]+sigma[k-1]*gsl_ran_gaussian(r,1.0);
          }
    }
    //Process observations 1
    t=0;
    //2.1 Compute the weights
    #pragma omp parallel for private(k,count,work,prob,m1,m2,sig1,sig2)
    for(i=0; i<*N; i++){
          work=0;
          for(k=0; k< vec_length[0]; k++){
	        work+= observations[*dobs*t+a_vec[k]]*theta[*d*i+k];
          }
          prob=exp(-work)/(1.0+exp(-work));
          count=vec_length[0];
          m1=0;
          for(k=0; k<vec_length[1]; k++){
	        m1+= observations[*dobs*t+m1_vec[k]]*theta[*d*i+count+k];
          }
          count+=vec_length[1];
          m2=0;
          for(k=0; k<vec_length[2]; k++){
	        m2+= observations[*dobs*t+m2_vec[k]]*theta[*d*i+count+k];
          }
          count+=vec_length[2];
          work=0;
          for(k=0; k<vec_length[3]; k++){
	        work+= observations[*dobs*t+sig1_vec[k]]*theta[*d*i+count+k];
          }
          sig1=exp(-work);
          count+=vec_length[3];
          work=0;
          for(k=0; k<vec_length[4]; k++){
	        work+= observations[*dobs*t+sig2_vec[k]]*theta[*d*i+count+k];
          }
          sig2=exp(-work);
          w[i]=log(prob*gsl_ran_gaussian_pdf(observations[*dobs*t]-m1,sig1)+(1.0-prob)*gsl_ran_gaussian_pdf(observations[*dobs*t]-m2,sig2));     
    }
    work=weight(w,  W, *N);
    //Compute the ESS
    ess_val=0;
    for(i=0;i<*N;i++){
            ess_val+=W[i]*W[i];
    }
    ess_val=1.0/ess_val;
    //Compute the theta_tilde
    for(k=0; k<*d;k++){
         res[*d*t+k]=0;
         for(i=0;i<*N;i++){
              res[*d*t+k]+=W[i]*theta[*d*i+k];
         }
    }
    //Process other observations
    for(t=1; t< *datalength; t++)
    {
             h_t=alpha[0]*pow(t,-alpha[1]);  //compute h_t
             if(ess_val> *ess_bound){ //no resampling
               if(tp[t-1]==1){  //Student move
                 for(i=0; i<*N; i++){
                       work=pow(gsl_ran_chisq(r, *nu),-0.5);
                       for(k=0; k<*d;k++){
                           theta[*d*i+k]+=h_t*gsl_ran_gaussian(r,1.0)*sqrt_nu*work;
                       }
                 }
              }
              else{ //Gaussian move
                 for(i=0; i<*N; i++){
                       for(k=0; k<*d;k++){
                           theta[*d*i+k]+=h_t*gsl_ran_gaussian(r,1.0);
                       }
                 }
              }
              #pragma omp parallel for private(k,count,work,prob,m1,m2,sig1,sig2)
              for(i=0; i<*N; i++){
                 if(theta[*d*i]<0){
                     w[i]=-INFINITY;
                 }
                 else{              
                    work=0;
                    for(k=0; k<vec_length[0]; k++){
	                 work+= observations[*dobs*t+a_vec[k]]*theta[*d*i+k];
                    }
                    prob=exp(-work)/(1.0+exp(-work));
                    count=vec_length[0];
                    m1=0;
                    for(k=0; k<vec_length[1]; k++){
	                 m1+= observations[*dobs*t+m1_vec[k]]*theta[*d*i+count+k];
                    }
                    count+=vec_length[1];
                    m2=0;
                    for(k=0; k<vec_length[2]; k++){
	                 m2+= observations[*dobs*t+m2_vec[k]]*theta[*d*i+count+k];
                    }
                    count+=vec_length[2];
                    work=0;
                    for(k=0; k<vec_length[3]; k++){
	                 work+= observations[*dobs*t+sig1_vec[k]]*theta[*d*i+count+k];
                    }
                    sig1=exp(-work);
                    count+=vec_length[3];
                    work=0;
                    for(k=0; k<vec_length[4]; k++){
	                work+= observations[*dobs*t+sig2_vec[k]]*theta[*d*i+count+k];
                    }
                    sig2=exp(-work);
                    w[i]+=log(prob*gsl_ran_gaussian_pdf(observations[*dobs*t]-m1,sig1)+(1.0-prob)*gsl_ran_gaussian_pdf(observations[*dobs*t]-m2,sig2));     
                 }
              }
           }
           else{//resampling
                for(i=0; i<*N; i++){
                      w[i]=gsl_rng_uniform(r);
                }
                
                SSPResampler(w,  N, W, J);

               // #pragma omp parallel for private(k)
                for(i=0; i<*N; i++){
                      for(k=0; k<*d;k++){
                          thetah[*d*i+k]=theta[*d*J[i]+k];
                      }
                }
                  
                if(tp[t-1]==1){
                     for(i=0; i<*N; i++){
                         work=pow(gsl_ran_chisq(r, *nu),-0.5);
                         for(k=0; k<*d;k++){
                             theta[*d*i+k]=thetah[*d*i+k]+h_t*gsl_ran_gaussian(r,1.0)*sqrt_nu*work;
                         }
                     }
                }
                else{
                    for(i=0; i<*N; i++){
                        for(k=0; k<*d;k++){
                           theta[*d*i+k]=thetah[*d*i+k]+h_t*gsl_ran_gaussian(r,1.0);
                        }
                    }
                }
                #pragma omp parallel for private(k,count,work,prob,m1,m2,sig1,sig2)
                for(i=0; i<*N; i++){
                  if(theta[*d*i]<0){
                      w[i]=-INFINITY;
                  }
                  else{              
                    work=0;
                    for(k=0; k<vec_length[0]; k++){
	                 work+= observations[*dobs*t+a_vec[k]]*theta[*d*i+k];
                    }
                    prob=exp(-work)/(1.0+exp(-work));
                    count=vec_length[0];
                    m1=0;
                    for(k=0; k<vec_length[1]; k++){
	                 m1+= observations[*dobs*t+m1_vec[k]]*theta[*d*i+count+k];
                    }
                    count+=vec_length[1];
                    m2=0;
                    for(k=0; k<vec_length[2]; k++){
	                 m2+= observations[*dobs*t+m2_vec[k]]*theta[*d*i+count+k];
                    }
                    count+=vec_length[2];
                    work=0;
                    for(k=0; k<vec_length[3]; k++){
	                 work+= observations[*dobs*t+sig1_vec[k]]*theta[*d*i+count+k];
                    }
                    sig1=exp(-work);
                    count+=vec_length[3];
                    work=0;
                    for(k=0; k<vec_length[4]; k++){
	                work+= observations[*dobs*t+sig2_vec[k]]*theta[*d*i+count+k];
                    }
                    sig2=exp(-work);
                    w[i]=log(prob*gsl_ran_gaussian_pdf(observations[*dobs*t]-m1,sig1)+(1.0-prob)*gsl_ran_gaussian_pdf(observations[*dobs*t]-m2,sig2));     
                 }
               }    
            }
            work=weight(w,  W, *N);
            ess_val=0;
            for(i=0;i<*N;i++){
                 ess_val+=W[i]*W[i];
            }
            ess_val=1.0/ess_val;

            for(k=0; k<*d;k++){
                  res[*d*t+k]=0;
                   for(i=0;i<*N;i++){
                       res[*d*t+k]+=W[i]*theta[*d*i+k];
                   }
            }                   
    }
    free(theta);
    theta=NULL;
    free(thetah);
    thetah=NULL;
    free(w);
    w=NULL;
    free(W);
    W=NULL;
    free(J);
    J=NULL;
    gsl_rng_free(r);

}





void PSMCO_mixture_C(int *cl,  int *dobs, double *ess_bound, double *mu, double *sigma, int *a_vec, int *m1_vec,int *m2_vec,
                     int *sig1_vec, int *sig2_vec, int *vec_length, int *seed, int *N, int *d, double *nu, int *datalength, 
                     double *alpha, int *tp,  double *observations, double *res){

    omp_set_dynamic(0);     		
    omp_set_num_threads(*cl);

    int t, i, k, count;
    double ess_val,  work, prob, m1,m2,sig1,sig2; 
    ess_val=*N;
    double proba=1.0/pow(*N,0.5); 

    double *theta=(double*)malloc(sizeof(double)*(*d*(*N)));
    double *thetah=(double*)malloc(sizeof(double)*(*d*(*N)));
    double *w=(double*)malloc(sizeof(double)*(*N));
    double *W=(double*)malloc(sizeof(double)*(*N));
    int *J=(int*)malloc(sizeof(int)*(*N));

   
    gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set(r, *seed);

    //sample from tilde{pi}_0
    for(i=0; i<*N; i++){
          theta[*d*i]=gsl_ran_exponential(r,mu[0]);
          for(k=1; k<*d; k++){
               theta[*d*i+k]=mu[k]+sigma[k-1]*gsl_ran_gaussian(r,1.0);
          }
    }
    //Process observations 1
    t=0;
    //2.1 Compute the weights
    #pragma omp parallel for private(k,count,work,prob,m1,m2,sig1,sig2)
    for(i=0; i<*N; i++){
          work=0;
          for(k=0; k< vec_length[0]; k++){
	        work+= observations[*dobs*t+a_vec[k]]*theta[*d*i+k];
          }
          prob=exp(-work)/(1.0+exp(-work));
          count=vec_length[0];
          m1=0;
          for(k=0; k<vec_length[1]; k++){
	        m1+= observations[*dobs*t+m1_vec[k]]*theta[*d*i+count+k];
          }
          count+=vec_length[1];
          m2=0;
          for(k=0; k<vec_length[2]; k++){
	        m2+= observations[*dobs*t+m2_vec[k]]*theta[*d*i+count+k];
          }
          count+=vec_length[2];
          work=0;
          for(k=0; k<vec_length[3]; k++){
	        work+= observations[*dobs*t+sig1_vec[k]]*theta[*d*i+count+k];
          }
          sig1=exp(-work);
          count+=vec_length[3];
          work=0;
          for(k=0; k<vec_length[4]; k++){
	        work+= observations[*dobs*t+sig2_vec[k]]*theta[*d*i+count+k];
          }
          sig2=exp(-work);
          w[i]=log(prob*gsl_ran_gaussian_pdf(observations[*dobs*t]-m1,sig1)+(1.0-prob)*gsl_ran_gaussian_pdf(observations[*dobs*t]-m2,sig2));     
    }
    work=weight(w,  W, *N);
    //Compute the ESS
    ess_val=0;
    for(i=0;i<*N;i++){
            ess_val+=W[i]*W[i];
    }
    ess_val=1.0/ess_val;
    //Compute the theta_tilde
    for(k=0; k<*d;k++){
         res[*d*t+k]=0;
         for(i=0;i<*N;i++){
              res[*d*t+k]+=W[i]*theta[*d*i+k];
         }
    }
   
    //Process other observations
    for(t=1; t< *datalength; t++)
    {
             if(ess_val> *ess_bound){ //no resampling
                for(i=0; i<*N; i++){
                     if(gsl_rng_uniform(r)<proba){
                       for(k=0; k<*d;k++){
                           theta[*d*i+k]+= *alpha*gsl_ran_gaussian(r,1.0);
                       }
                    }
                }
              #pragma omp parallel for private(k,count,work,prob,m1,m2,sig1,sig2)
              for(i=0; i<*N; i++){
                 if(theta[*d*i]<0){
                     w[i]=-INFINITY;
                 }
                 else{              
                    work=0;
                    for(k=0; k<vec_length[0]; k++){
	                 work+= observations[*dobs*t+a_vec[k]]*theta[*d*i+k];
                    }
                    prob=exp(-work)/(1.0+exp(-work));
                    count=vec_length[0];
                    m1=0;
                    for(k=0; k<vec_length[1]; k++){
	                 m1+= observations[*dobs*t+m1_vec[k]]*theta[*d*i+count+k];
                    }
                    count+=vec_length[1];
                    m2=0;
                    for(k=0; k<vec_length[2]; k++){
	                 m2+= observations[*dobs*t+m2_vec[k]]*theta[*d*i+count+k];
                    }
                    count+=vec_length[2];
                    work=0;
                    for(k=0; k<vec_length[3]; k++){
	                 work+= observations[*dobs*t+sig1_vec[k]]*theta[*d*i+count+k];
                    }
                    sig1=exp(-work);
                    count+=vec_length[3];
                    work=0;
                    for(k=0; k<vec_length[4]; k++){
	                work+= observations[*dobs*t+sig2_vec[k]]*theta[*d*i+count+k];
                    }
                    sig2=exp(-work);
                    w[i]+=log(prob*gsl_ran_gaussian_pdf(observations[*dobs*t]-m1,sig1)+(1.0-prob)*gsl_ran_gaussian_pdf(observations[*dobs*t]-m2,sig2));     
                 }
              }
           }
           else{//resampling
                for(i=0; i<*N; i++){
                      w[i]=gsl_rng_uniform(r);
                }
                
                SSPResampler(w,  N, W, J);

               // #pragma omp parallel for private(k)
                for(i=0; i<*N; i++){
                      for(k=0; k<*d;k++){
                          thetah[*d*i+k]=theta[*d*J[i]+k];
                      }
                }
                  
                for(i=0; i<*N; i++){
                     if(gsl_rng_uniform(r)<proba){
                        for(k=0; k<*d;k++){
                           theta[*d*i+k]=thetah[*d*i+k]+*alpha*gsl_ran_gaussian(r,1.0);
                        }
                     }
                     else{
                        for(k=0; k<*d;k++){
                           theta[*d*i+k]=thetah[*d*i+k];
                        }   
                    }
                 }
                #pragma omp parallel for private(k,count,work,prob,m1,m2,sig1,sig2)
                for(i=0; i<*N; i++){
                  if(theta[*d*i]<0){
                      w[i]=-INFINITY;
                  }
                  else{              
                    work=0;
                    for(k=0; k<vec_length[0]; k++){
	                 work+= observations[*dobs*t+a_vec[k]]*theta[*d*i+k];
                    }
                    prob=exp(-work)/(1.0+exp(-work));
                    count=vec_length[0];
                    m1=0;
                    for(k=0; k<vec_length[1]; k++){
	                 m1+= observations[*dobs*t+m1_vec[k]]*theta[*d*i+count+k];
                    }
                    count+=vec_length[1];
                    m2=0;
                    for(k=0; k<vec_length[2]; k++){
	                 m2+= observations[*dobs*t+m2_vec[k]]*theta[*d*i+count+k];
                    }
                    count+=vec_length[2];
                    work=0;
                    for(k=0; k<vec_length[3]; k++){
	                 work+= observations[*dobs*t+sig1_vec[k]]*theta[*d*i+count+k];
                    }
                    sig1=exp(-work);
                    count+=vec_length[3];
                    work=0;
                    for(k=0; k<vec_length[4]; k++){
	                work+= observations[*dobs*t+sig2_vec[k]]*theta[*d*i+count+k];
                    }
                    sig2=exp(-work);
                    w[i]=log(prob*gsl_ran_gaussian_pdf(observations[*dobs*t]-m1,sig1)+(1.0-prob)*gsl_ran_gaussian_pdf(observations[*dobs*t]-m2,sig2));     
                 }
               }    
            }
            work=weight(w,  W, *N);
            ess_val=0;
            for(i=0;i<*N;i++){
                 ess_val+=W[i]*W[i];
            }
            ess_val=1.0/ess_val;

            for(k=0; k<*d;k++){
                  res[*d*t+k]=0;
                   for(i=0;i<*N;i++){
                       res[*d*t+k]+=W[i]*theta[*d*i+k];
                   }
            }                   
    }
    free(theta);
    theta=NULL;
    free(thetah);
    thetah=NULL;
    free(w);
    w=NULL;
    free(W);
    W=NULL;
    free(J);
    J=NULL;
    gsl_rng_free(r);

}



void Lik_multimodal(int *cl, int *N, int *d , int *t,  double *theta, double *observations, double *res)
 {
   omp_set_dynamic(0);     		
   omp_set_num_threads(*cl);

   int s, i, k;
   double work, bound1=-21.0, bound2=19.0;
 

  #pragma omp parallel for private(s,work,k)
  for(i=0; i<*N; i++){
        work=0;
        for(k=0;k<*d;k++){
            if(theta[*d*i+k]<bound1 || theta[*d*i+k]>bound2){
               work=1.0;
            }
        }
        if(work==1.0){
           res[i]= -INFINITY;
        }
        else{
        res[i]=0;
          for(s=0; s<*t; s++){
             work=observations[(*d+1)*s];
          
             for(k=0; k<*d; k++){
	          work+= -exp(-observations[(*d+1)*s+(k+1)]*pow(theta[*d*i+k], 2.0))-observations[(*d+1)*s+(k+1)]*theta[*d*i+*d-k-1];
             }
             res[i]+=-0.5*fabs(work);
        }  
       }
  }
}










