#include <iostream>
#include <cmath>
#include <random>
#include "MuS.h"
#include "Lib.h"
#include "gig.h"
#include "Rlob.h"
#include "global.h"
//extern "C" {
#include <R.h>
//};

// public methods -------------------------------------------------
void MuS::drawPost()
{
    mu=post_m;
}


double MuS::getLogILik(){

  double logres = 0.0;
  // Validate input parameters first
  if (std::isnan(rhk) || std::isnan(shk) || std::isinf(rhk) || std::isinf(shk) || -aa + rhk <= 0. || shk <= 0.) {
    //std::cerr << "Invalid value for logres" << std::endl;
    return -10000.0;  // Sentinel value for invalid inputs
  }

  double zValue1, zValue2, zValue3;

  try {
    zValue1 = computeZ(-aa + rhk, 2.0 * bb, 2.0 * shk);
    //std::cout << "zValue1 = " << zValue1 << std::endl;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return -10000.0;
  }

    zValue2 = computeZ(aa + rhk, 0, 2.0 * (bb + shk));
    zValue3 = computeZ(aa, 0, 2.0 * bb);
  // Check for non-positive values which are invalid for logarithm
  if (zValue1 <= 0 || zValue2 <= 0 || zValue3 <= 0) {
    std::cerr << "Invalid value for logarithm encountered in Z function." << std::endl;
    return -10000.0;  // Sentinel value for invalid calculation
  }

  try {
    double lognumerator = std::log(zValue1 + zValue2);
    double logdenominator = std::log(2.0 * zValue3);
    logres = lognumerator - logdenominator;
    //std::cout << "logres = " << logres << std::endl;

  } catch (const std::exception& e) {
    std::cerr << "Exception during log calculations: " << e.what() << std::endl;
    return -10000.0;  // Sentinel value for exception in logarithm calculation
  }

  // Check for non-finite result
  if (std::isinf(logres) || std::isnan(logres)) {
    std::cerr << "Non-finite result encountered in logres calculation." << std::endl;
    return -10000.0;  // Sentinel value for invalid result
  }

  return logres;
}






void MuS::updatepost() {
  int i;
  if(nob) {
    //ch=1.0;
    rhk=0.0;
    shk=0.0;
    for(i=1;i<=nob;i++){
      //ch *= wi[indices[i]]*pow(exp(y[indices[i]]),de[indices[i]]);
     rhk += de[indices[i]];
      shk += exp(y[indices[i]]) * v[indices[i]];
    }
    //std::cout << "rhk = " << rhk << std::endl;
    //std::cout << "shk = " << shk << std::endl;
    if (std::isnan(rhk) || std::isnan(shk) || shk <= 0. || ( 2*shk == 0. && rhk-aa >= 0.)) {
      post_m = 0.0;
   } else {
      double U;
      U = Bern(0.5);
      double W1;
      W1 = geGammaSample(aa+rhk, bb+shk);
      double W2;
      Random randomGenerator(0);
      W2= randomGenerator.gig( rhk-aa, 2*bb, 2*shk);
      post_m =std::log(U*W1+(1-U)*W2);
      }
  }else{
    post_m=0.0;
  }
}


void MuS::setData(int nob, double **x, double *y, double*de, double*wi, double*v,
                  int *indices, double *w)
{
  this->nob=nob;
  this->y=y;
  this->de=de;
  this->wi=wi;
  this->v=v;
  this->indices=indices;
  updatepost();
}



double* MuS::getParameterEstimate()
{
   double *rv = new double[2];
   rv[1] = post_m;
   return rv;
}
double* MuS::getFits(int np, double **xpred, int* indpred)
{
   double *rv = new double[np+1];
   for(int i=1;i<=np;i++) rv[i]=post_m;
   return rv;
}
void MuS::toScreen() const
{
   Rprintf("defunct MuS::toScreen called");
   /*
   std::cout << "MuS EndNodeModel:\n";
   std::cout << "  Normal mean given the standard deviation\n";
   std::cout << "mu: " << mu << " ,sigma: " << std::sqrt(sigma2) << "\n";
   std::cout << "prior sd: " << 1.0/std::sqrt(a) << "\n";
   std::cout << "nob: " << nob << std::endl;
   if(nob) {
      std::cout << "posterior mean and sd: " <<
          post_m << " , " << post_s << "\n";
      std::cout << "nob, ybar,s2: " << nob << ", " <<  ybar << ", " << s2 << std::endl;
   }
   */
}
