// [[Rcpp::depends(BH)]]

#include<Rcpp.h>
#include<R.h>
#include<algorithm>
#include<boost/random.hpp>
#include<chrono>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<fstream>
#include<iomanip>
#include<iostream>
#include<math.h>
#include<random>
#include<sstream>
#include<string>
#include<thread>
#include<vector>
#include<stdlib.h> 


#ifndef PDF_NORMAL_H
#define PDF_NORMAL_H
// [[Rcpp::export]]
double pdf_normal(double mu, double sigma, double x)
{
  double part1 = 1 / sqrt(2 * M_PI * sigma * sigma);
  double part2 = exp( (-1 * (x - mu) * (x - mu) ) / ( 2 * sigma * sigma ));
  double answer = part1 * part2;
  return answer;
}
#endif

#ifndef CDF_NORMAL_H
#define PDF_NORMAL_H
// [[Rcpp::export]]
double cdf_normal(double mu, double sigma, double x)
{
  double answer;
  answer = erf( (x - mu ) / ( sigma * sqrt(2) ) );
  answer ++;
  answer = answer / 2;
  return answer;
}
#endif


#ifndef PDF_F_T_H
#define PDF_F_T_H
// [[Rcpp::export]]
double f_t (double t, double v, double b, double z, double s)
{
  double outvar = pdf_normal (v, s, (b - z) / t) * (b - z) / t / t;
  return outvar; 
}
#endif

#ifndef CDF_F_T_H
#define CDF_F_T_H
// [[Rcpp:export]]
double F_t (double t, double v, double b, double z, double s)
{
  double outvar = cdf_normal (v, s, (b - z) / t);
  return outvar;
}
#endif


#ifndef LBAR_H
#define LBAR_H
// [[Rcpp::export]]
double lbar (double t, double v1, double v2, double z1, double z2, double b, double s1, double s2)
{
  double term1 = f_t (t, v1, b, z1, s1);
  // Rprintf("\n term1: %f", term1);
  double term2 = F_t (t, v2, b, z2, s2);
  // Rprintf("\n term2: %f", term2);
  return term1*term2;
}
#endif


#ifndef LOGDENS_LBAR_H
#define LOGDENS_LBAR_H
// [[Rcpp::export]]
double logdens_lbar (double rt, int resp, int stim, double v1, double v2, double z1, double z2, double b, double s1 = 1, double s2 = 1)
{
  double dens = 0.0;
  double lbartemp = 0.0;
  if (resp == stim) // MAY NEED TO CHANGE
  {
    lbartemp = lbar(rt, v1, v2, z1, z2, b, s1, s2);
  }else{
    lbartemp = lbar(rt, v2, v1, z2, z1, b, s2, s1);
  }
  if(lbartemp <= 0.0){
    // Rprintf("true");
    dens = -739.9974;
  }else{
    dens = std::log(lbartemp);
  }
  return(dens);
}
#endif

// ------------------------------------------------------------------------------------------------

#ifndef LOGDENS_BALBA_H
#define LOGDENS_BALBA_H
// [[Rcpp::export]]
double logdens_balba (std::vector < double > rt, 
                      std::vector < int > resp, 
                      std::vector < int > stim, 
                      double v1, 
                      double v2, 
                      double A, 
                      std::vector < double >  alpha, 
                      std::vector < double >  beta, 
                      double b, 
                      int mc_size = 100, 
                      int seed = 123)
{
  int N = rt.size();
  double z1;
  double z2;
  double dens = 0.0;
  
  typedef boost::mt19937 RNGType; 	// select a generator, MT good default
  RNGType rng(seed);			// instantiate and seed
  RNGType rng2(seed+1);			// instantiate and seed
  
  for(int i = 0; i < N; i++){ // loop over all RTs
    // determine starting point distribution for racer 1
    boost::random::beta_distribution<double> distributionz1(alpha[i], beta[i]); 
    boost::variate_generator< RNGType, boost::random::beta_distribution<> > rngz1(rng, distributionz1);
    
    // determine starting point distribution for racer 2
    boost::random::beta_distribution<double> distributionz2(beta[i], alpha[i]);
    boost::variate_generator< RNGType, boost::random::beta_distribution<> > rngz2(rng2, distributionz2);
    
    // sample likelihoods from starting point distributions
    for(int j = 0; j < mc_size; j++){
      z1 = rngz1() * A; // sample racer 1 starting point with upper limit A
      z2 = rngz2() * A; // sample racer 2 starting point with upper limit A
      // Rprintf("\n z1: %f", z1);
      // Rprintf("\n z2: %f", z2);
      dens += logdens_lbar(rt[i], resp[i], stim[i], v1, v2, z1, z2, b)/mc_size;
    }
  }
  return (dens);
}  
#endif


#ifndef RBALBA_H
#define RBALBA_H
// [[Rcpp::export]]
Rcpp::NumericVector rbalba (double v1, 
                            double v2, 
                            double A = .75, 
                            double alpha = 1, 
                            double beta = 1, 
                            double b = 1, 
                            double s1 = 1, 
                            double s2 = 1, 
                            int seed = 123)
{
  // initiate starting point and racer end times
  double z1;
  double z2;
  double t_1;
  double t_2;
  
  typedef boost::mt19937 RNGType; 	// select a generator, MT good default
  RNGType rng(seed);			// instantiate and seed
  RNGType rng2(seed+1);			// instantiate and seed
  
  // create random sampler for racer 1
  boost::normal_distribution<double> distributionV1(v1, s1);
  boost::variate_generator< RNGType, boost::normal_distribution<> > rngv1(rng, distributionV1);
  
  // create random sampler for racer 2
  boost::normal_distribution<double> distributionV2(v2, s2);
  boost::variate_generator< RNGType, boost::normal_distribution<> > rngv2(rng2, distributionV2);
  
  // create random sampler for starting point 1
  boost::random::beta_distribution<double> distributionz1(alpha, beta);
  boost::variate_generator< RNGType, boost::random::beta_distribution<> > rngz1(rng, distributionz1);
  
  // create random sampler for starting point 2
  boost::random::beta_distribution<double> distributionz2(beta, alpha);
  boost::variate_generator< RNGType, boost::random::beta_distribution<> > rngz2(rng2, distributionz2);
  
  double v1_i = -1;
  double v2_i = -1;
  // Sample positive racer 1
  while(v1_i < 0){
    v1_i = rngv1();
  }
  // Sample positive racer 2
  while(v2_i < 0){
    v2_i = rngv2();
  }
  z1 = rngz1() * A; // Sample starting point 1
  z2 = rngz2() * A; // Sample starting point 2
  
  // Rprintf("\n z1: %f", z1);
  // Rprintf("\n z2: %f", z2);
  // Rprintf("\n v1: %f", v1_i);
  // Rprintf("\n v2: %f", v2_i);
  
  t_1 = (b - z1) / v1_i; // determine racer 1 time
  t_2 = (b - z2) / v2_i; // determine racer 2 time
  if (t_1 < t_2)
  {
    // racer 1 wins
    Rcpp::NumericVector v = {t_1, 1};
    return v;
  }
  else
  { 
    // racer 2 wins
    Rcpp::NumericVector v = {t_2, 2};
    return v;
  }
}
#endif