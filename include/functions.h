


#include<stdlib.h>  		
#include<stdio.h>	
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include<gsl/gsl_eigen.h>

//SSP resampling
void SSPResampler(double* points, int* N,  double *W, int *J);
void SSPResampler_R(double* points, int* N, double *W, double *J);


//useful functions
double weight(double*, double*, int);
double maxmin(double*,int,int);
double exponen(double);
double sum2(double *vect,int *N);         
                 
             
// CQR model
void ADA_censored(double *quantile, double *, double *theta_in, int *M, int *d, int *datalength,  double *alpha,   double *observations, double *res);     
void censored_online_jitter(int *cl, double *q, double *ess_bound, double *mu, double *sigma, int *seed, int *N, int *d, int *datalength, 
                   double *alpha,  double *observations, double *res);
void censored_online(int *cl, double *q, double *ess_bound, double *mu, double *sigma, int *seed, int *N, int *d, double *nu, int *datalength, 
                   double *alpha, int *tp,  double *observations, double *res);
      

//Toy multimodal example
void GPFSO_multimodal(int *cl,  double *ess_bound, int *seed, int *N, int *d, double *nu, int *datalength, 
                   double *alpha, int *tp,  double *observations, double *res);
void SGD_multi(int *d, double *theta_in, int *n_val, int *datalength, double *alpha,   double *observations, double *res);                                      
void Lik_multimodal(int *cl, int *N, int *d , int *t,  double *theta, double *observations, double *res);       
//SAGM model
void GPFSO_mixture_C(int *cl,  int *dobs, double *ess_bound, double *mu, double *sigma, int *a_vec, int *m1_vec, int *m2_vec,
                     int *sig1_vec, int *sig2_vec, int *vec_length, int *seed, int *N, int *d, double *nu, int *datalength, 
                     double *alpha, int *tp,  double *observations, double *res);
                   
void PSMCO_mixture_C(int *cl,  int *dobs, double *ess_bound, double *mu, double *sigma, int *a_vec, int *m1_vec,int *m2_vec,
                     int *sig1_vec, int *sig2_vec, int *vec_length, int *seed, int *N, int *d, double *nu, int *datalength, 
                     double *alpha, int *tp,  double *observations, double *res);
                              
                   
// g-and-k distribution               
double  gandkcdf_(double y, double *theta, int *maxsteps, double *tolerance, double *lower, double *upper);
void Lik_gandk(int *cl, int *N,    double *theta, double *observations, int *maxsteps, double *tolerance, double *lower, double *upper, double *res);
void gandkcdf2_(double *y, double *theta, int *maxsteps, double *tolerance, double *lower, double *upper, double *res);
void Lik_gandk_multi_obs(int *cl,  double *theta, int *datalength, double *observations, int *maxsteps, double *tolerance, double *lower, double *upper, double *res);
void Lik_gandk_multi_obs_N(int *cl, int *N,    double *theta, int *datalength, double *observations, int *maxsteps, double *tolerance, double *lower, double *upper, double *res);




                   
                   
                   
                   
                   
                   
                   
                   
                   

