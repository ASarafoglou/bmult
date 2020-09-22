
#include <Rcpp.h>
#include <mpfr.h>

using namespace Rcpp;

mpfr_t BigHi, BigLo, BigA, BigLoA, BigMinusZi, BigHighBound, BigUni0, BigUni1, BigTheta, BigBetaMinusOne, BigExponent, BigResult;

// [[Rcpp::export]]
void initBigNumbers(int bitsPrecision)
{   
    auto setPrec = [&](mpfr_t & setMyPrecision){ mpfr_init2(setMyPrecision, bitsPrecision); };
    
    setPrec(BigHi);
    setPrec(BigLo);
    setPrec(BigA);
    setPrec(BigLoA);
    setPrec(BigMinusZi);
    setPrec(BigHighBound);
    setPrec(BigUni0);
    setPrec(BigUni1);
    setPrec(BigTheta);
    setPrec(BigBetaMinusOne);
    setPrec(BigExponent);
    setPrec(BigResult);
    
}

double truncatedSamplingSubFunction()
{
  mpfr_pow(BigLoA, BigLo, BigA, MPFR_RNDD); //Lo^a[k]
  mpfr_pow(BigResult, BigHi, BigA, MPFR_RNDD); //Hi^a[k]
  
  mpfr_sub(BigResult, BigResult, BigLoA, MPFR_RNDD); //(Hi^a[k]) - (Lo^a[k])
  mpfr_mul(BigResult, BigResult, BigUni1, MPFR_RNDD); //runif(1) * ((Hi^a[k]) - (Lo^a[k]))
  mpfr_add(BigResult, BigResult, BigLoA, MPFR_RNDD); //(runif(1) * ((Hi^a[k]) - (Lo^a[k]))) + (Lo^a[k])
  
  mpfr_d_div(BigA, 1.0, BigA, MPFR_RNDD); // 1/a[k]
  mpfr_pow(BigResult, BigResult, BigA, MPFR_RNDD); //((runif(1) * ((Hi^a[k]) - (Lo^a[k]))) + (Lo^a[k])) ^ (1/a[k])
  
  double result = mpfr_get_d(BigResult, MPFR_RNDD);

return result;
}

// [[Rcpp::export]]
double truncatedSamplingSubiterationBinomialCDF(double uniformSample, double alpha, double Lo, double Hi) 
{
  mpfr_set_d(BigUni1,         uniformSample,  MPFR_RNDD);
  mpfr_set_d(BigA,            alpha,          MPFR_RNDD);
  mpfr_set_d(BigLo,           Lo,             MPFR_RNDD);
  mpfr_set_d(BigHi,           Hi,             MPFR_RNDD);
  
  return truncatedSamplingSubFunction();
}

// [[Rcpp::export]]
double truncatedSamplingSubiteration(double uniformSample0, double uniformSample1, double minusZi, double Lo, double ai, bool thereIsAHigherBound, double theHigherBound ) 
{
    mpfr_set_d(BigUni0,         uniformSample0, MPFR_RNDD);
    mpfr_set_d(BigUni1,         uniformSample1, MPFR_RNDD);
    mpfr_set_d(BigMinusZi,      minusZi,        MPFR_RNDD);
    mpfr_set_d(BigLo,           Lo,             MPFR_RNDD);
    mpfr_set_d(BigA,            ai,             MPFR_RNDD);
    mpfr_set_d(BigHighBound,    theHigherBound, MPFR_RNDD);

    
    mpfr_exp(BigHi,   BigMinusZi,          MPFR_RNDD);
    mpfr_mul(BigHi,   BigHi,      BigUni0, MPFR_RNDD);
    mpfr_log(BigHi,   BigHi,               MPFR_RNDD);
    mpfr_mul_d(BigHi, BigHi,      -1.0,    MPFR_RNDD);

    if(thereIsAHigherBound)
        mpfr_min(BigHi, BigHi, BigHighBound, MPFR_RNDD);

    return truncatedSamplingSubFunction();
}

// [[Rcpp::export]]
double truncatedSamplingSubiterationBinomialY(double uniformSample, double theta, double betaMinusOne) 
{
  mpfr_set_d(BigUni0,         uniformSample,  MPFR_RNDD);
  mpfr_set_d(BigTheta,        theta,          MPFR_RNDD);
  mpfr_set_d(BigBetaMinusOne, betaMinusOne,   MPFR_RNDD);
  
  //Rcout << "BigBetaMinusOne: " << mpfr_get_d(BigBetaMinusOne, MPFR_RNDD) << std::endl;
  
  mpfr_d_div(BigExponent, 1.0, BigBetaMinusOne, MPFR_RNDD); // 1/(beta - 1)
  mpfr_d_sub(BigResult,    1.0,       BigTheta,         MPFR_RNDD); //(1 - theta[k])
  mpfr_pow(BigResult,      BigResult, BigBetaMinusOne,  MPFR_RNDD); //(1 - theta[k])^(beta - 1)
  mpfr_mul(BigResult,      BigResult, BigUni0,          MPFR_RNDD); //runif(1) * (1 - theta[k])^(beta - 1)
  mpfr_pow(BigResult,      BigResult, BigExponent,      MPFR_RNDD); //(runif(1) * (1 - theta[k])^(beta - 1))^exponent
  mpfr_d_sub(BigResult,    1.0,       BigResult,         MPFR_RNDD); //1 - (runif(1) * (1 - theta[k])^(beta - 1))^exponent
  
  double result = mpfr_get_d(BigResult, MPFR_RNDD);
  
  return result;
}
