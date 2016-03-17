// Compile for R with "R CMD SHLIB prodpval_stat_integrand.c"
// Compiler flags can be changed with "export MAKEFLAGS='CFLAGS=-O3\ -fno-strict-aliasing\ -Wall'"
#include <math.h>
#include <float.h>    // FLT_MIN
#include <limits.h>
#include <stdio.h>     // printf

// Calculate the product of p-values statistics:
// Given N independent p-values (distributed uniformly over [0,1]), p1,...,pN, denote the sorted values by p'1,...,p'N.
// Calculate the probabilty for the geometric mean of the K lowest p-values to be lower than some q
//
//   Prob( Prod_{k=1,K} p'k < q^K) = N!/(N-K-1)!/K! * Integral_{0,1} (1-p)^(N-K-1) p^K P(K,(q/p)^K) dp 
//
// where
//   P(K,(q/p)^K) = min{1, (q/p)^K Sum_{k=0,K-1} (-K*log(q/p))^k/k! }
//
// is the probability that the product of K p-values out of K will be lower than (q/p)^K and
// N!/(N-K-1)!/K! * (1-p)^(N-K-1) p^K dp is the probability that the (K+1)th-ranked P-value is between p and p+dp.

// This function uses the approximations N-K >> K, e.g. N-K >= 10*K, and Pval <<1.
// The last approximation is equivalent to q << K/(N-K-1)/e, since by chance q is expected to be around K/N/e:
// q^K ~ K!/N^K ~ (K/e)^K/N^K, hence q ~ K/N/e. 
// The Laplace approximation is used to approximate the integrand by a Gaussian function, which can be integrated
// analytically. For this purpose, the point p1 where the log integrand reaches its maximum is determined.
// The height of the Gaussian is then f(p1), where f(p) is the integrand
//
//    f(p) = (1-p)^(N-K-1) * p^K * min{ 1, (q/p)^K * sum_{k=0,K-1} (-K*log(q/p))^k/k! }
//
// The standard deviation of the Gaussian is given by 1/sigma^2 = -1 / ( d^2/dp^2 log f(p) ).
// The integral is then  f(p1) * sqrt(2*pi*sigma^2)
//
//

double ddlogf(const int N, const double Kd, const double q, const double y, const double p);
double logsum(const int K, const double y);
double prodpval_stat(const int N, const int *K, const double *q);

// Calculate the approximate second derivative d^2/dp^2 log f(p) of the log integrand 
//    f(p) = (q/p)^K * sum_{k=0,K-1} (-K*log(q/p))^k/k!  wrt p
inline double ddlogf(const int N, const double Kd, const double q, const double y, const double p)
{
  double zK = pow(y,-Kd);
  double zK1 = y*zK;
  return ((double)N-Kd-1)/(1-p)/(1-p) - Kd/p/p/y * (1-zK1)/(1-zK) + Kd*Kd/(Kd-1)/p/p * ( (Kd-1)*zK/y/(1-zK) - (1-zK1)/y/y/(1-zK)/(1-zK) )*(1+(Kd-1)*zK); 
}


// Calculate the sum 1 + (K-1)/y + (K-1)(K-2)/y^2 + ... +(K-1)!/y^(K-1)
inline double logsum(const int K, const double y)
{
  double sum=1.0;
  double term=1.0;
  int k;
  
  for (k=K-1; k>=1 && term>1E-4; k--)  // 1E-4????
    { 
      term *= (double)k/y; 
      sum += term;
    }
  return log(sum);
}

inline double prodpval_stat(const int N, const int *K, const double *q)
{
  double res;
  double Kd = (double)*K;
  if (*K>N)
    {
      res = -DBL_MAX; // Error: return minus largest float
    }
  else if (*K==N)
    {
      // Calculate log P-value of product of N independent P-values, which is the log of 
      //   q^K sum_{k=0,K-1} (-K log q)^k/k! = q^K y^(K-1)/(K-1)! * [ 1 + (K-1)/y + (K-1)(K-2)/y^2 + ... + (K-1)!/y^(K-1) ]
      // with y = -K log q.
      res = Kd*log(*q) + (Kd-1)*log(-Kd*log(*q)) - lgamma(Kd) + logsum(*K,-Kd*log(*q));
    }
  else if (*K==1) // res = log ( 1 - (1-q)^N )
    {
      if (*q>1E-9)
     res = -pow(1-*q,N);
      else
     res = -exp(-(*q)*N);
      if (res<-1E-9) res = log(1+res);
    }
  else  // 1<K<N (general case)
    {
      const double Kd = (double)*K;
      const double p0 = Kd/(double)(N - 1);
      if (log(p0/(*q)) < 1.1) {res = 0.0; return res;} // not significant

      double p1prev = 1000;
      double p1 = sqrt(p0*(*q)*2.71828);
      double y, zK, p0q, sigma2;
      int i = 0;

      // Estimate position p1 of maximum of integrand iteratively
      while (fabs(1-p1prev/p1) > 1E-3) 
	{
	  if (i>15) {
	    printf("Warning: p1 iteration in logpval didn't converge. i=%i p1prev=%10.2E p1=%10.2E N=%i K=%i q=%10.2E\n",i,p1prev,p1,N,*K,*q);
	    break;}
	  p1prev = p1;
	  y = log(p1/(*q)) * Kd/(Kd-1);
	  zK = pow(y,-Kd);
	  p0q = p0*(1-y*zK)/(1-zK)*Kd/(Kd-1);
	  p1 = p0q / (y + p0q);
	  if (++i == 2) p1 = sqrt(p1prev*p1);
	}
      
      // Estimate integral
      p1 = sqrt(p1prev*p1);
      y = log(p1/(*q)) * Kd;
      sigma2 = -1/ddlogf(N,Kd,*q,y/(Kd-1),p1);
      if (sigma2 < 0) res  = 0.0; // not significant
      else 
    	  res = ((double)N-Kd-1.0)*log(1-p1) + Kd*log(*q) + (Kd-1)*log(y) - lgamma(Kd) + logsum(*K,y) + 0.5*log(6.283185*sigma2) + lgamma(N+1) - lgamma(N-Kd) - lgamma(Kd+1);
	  }

  return res;
}


