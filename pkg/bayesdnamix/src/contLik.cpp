//This script contains functions for calculating several conditional likelihoods.
//calculating likelihood of data with genotypes marginalised out
//cd ~/Dropbox/Forensic/MixtureProj/myDev/quantLR/bayesdnamix0
//R CMD SHLIB contLik.cpp -llapack -lblas

#include <cmath>
#include <stdio.h>
#include <RcppArmadillo.h> //require RcppArmadillopackage and Namespaced defined
//#include <armadillo> 
using namespace std;
using namespace arma;
const double PIVAL = std::acos(0.0)*2;
//#define DEBUG

 //PE - conditional liklihood value
 //theta - parameter in model
 //model - integer for selected data-model: {1=mixsep,2=strmix}
 //nC - number of contributors
 //nL - number of loci
 //nA - number of allele at each loci
 //nQi - number of genotype combinations at each loci
 //pG - matrix of genotype-probabilities for each combination and loci
 //CnQ - cumulative index for genotype-combintions 
 //allY - vector of peak heights
 //CnA - Cumulative index for peak heights
 //allX - matrix of genotype-contribution
 //sY - sum peak height at each loci

extern "C" {

//function for calculating likelihood of data with genotypes marginalised out
void contlikC(double *PE, double *theta, int *model, int *nC, int *nL, int *nA, int *nQi, double *pG, int *CnQ, double *allY, int *CnA, double *allX, int *cdfX, double *sY) {
 int i,j,k,cc;
 Col<double> *Yi; 
 Col<double> *omega; //mixture proportion 
 Mat<double> *Xij; //(ni x C) - matrix. Taking out columns to use for restriction
 omega = new Col<double>( theta, *nC-1, false); //insert mixture proportion: (C-1)-vector

 //temporary variables:
 Col<double> mui;
 Col<double> ri; //column vector
 Col<double> tmp; //column vector
 double Di; //residual sum squared
 double lik;
 double bigsum;
 double konstant;
 double tau; //variance parameter in model
 konstant = sqrt(2*PIVAL); 
 tau = theta[*nC-1]; //last element is variance parameter

 //SELECTION OF MODELS:
 //MIXSEP: 
 if(*model==1) {  //r ~ N(muij,tau)
  double sigma,sigmasq; //square root of variance parameter
  cc = 0; //counter for X-matrix and O-matrix
  for(i=0; i<*nL; i++) {
   sigmasq = tau*sY[i];
   sigma = sqrt(sigmasq);
   Yi =  new Col<double>( &allY[CnA[i]], nA[i], false); //init datavector
   bigsum = 0.0; //sum over all genotypes
   for(j=0; j<nQi[i]; j++) { //for all genotypes: Init space and insert relevant matrice directly
    Xij = new Mat<double>( &allX[cdfX[cc]] ,nA[i],*nC,false); //init. Xij-matrix
    mui = sY[i]/2*((Xij->cols(0,*nC-2))*(*omega) + Xij->col(*nC-1)); //mean is a vector
    ri = *Yi - mui; //residual is a vector
//    Di = as_scalar(trans(ri)*(ri))/sY[i]; //mahalonobis distance 
//    Di = as_scalar(sum(pow(ri,2)))/sY[i]; //mahalonobis distance 
    Di = dot(ri,ri)/sigmasq; //FASTEST!!
    lik = exp(-0.5*Di)/pow(konstant*sigma,nA[i]); 
    bigsum = bigsum + lik*pG[CnQ[i]+j]; //# add likval
    cc = cc + 1; //update counter
    delete Xij;
   }//end for each j: genotype combination
   *PE = (*PE)*bigsum;
   delete Yi;
  }//end for each loci i:
 } //end for model mixsep

 //Log-normal (logged mixsep)
 else if(*model==2) { //log(r) ~ N(log(muij),tau)
  double sigma,sigmasq; //square root of variance parameter
  double yprod;
  cc = 0; //counter for X-matrix and O-matrix
  for(i=0; i<*nL; i++) {
   sigmasq = tau; //no covariate weight
   sigma = sqrt(sigmasq);
   Yi =  new Col<double>( &allY[CnA[i]], nA[i], false); //init datavector
   yprod = prod(*Yi); //product of y
   bigsum = 0.0; //sum over all genotypes
   for(j=0; j<nQi[i]; j++) { //for all genotypes: Init space and insert relevant matrice directly
    Xij = new Mat<double>( &allX[cdfX[cc]] ,nA[i],*nC,false); //init. Xij-matrix
    mui = sY[i]/2*((Xij->cols(0,*nC-2))*(*omega) + Xij->col(*nC-1)); //mean is a vector
    ri = log(*Yi) - log(mui); //residual is vector
    Di = dot(ri,ri)/sigmasq; //FASTEST!!
    lik = exp(-0.5*Di); 
    bigsum = bigsum + lik*pG[CnQ[i]+j]; //# add likval
    cc = cc + 1; //update counter
    delete Xij;
   }//end for each j: genotype combination
   *PE = (*PE)*bigsum/pow(konstant*sigma,nA[i])/yprod; //add scaling-constant
   delete Yi;
  }//end for each loci i:
 } //end for model mixsep


 //strmix0
 else if(*model==3) { //r~N(muij,tau/muij)
  Col<double> sigmavec; //column vector used for variances
  double yprod;
  cc = 0; //counter for X-matrix and O-matrix
  for(i=0; i<*nL; i++) {
   Yi =  new Col<double>( &allY[CnA[i]], nA[i], false); //init datavector
   yprod = prod(*Yi); //product of y
   bigsum = 0.0; //sum over all genotypes
   for(j=0; j<nQi[i]; j++) { //for all genotypes: Init space and insert relevant matrice directly
    Xij = new Mat<double>( &allX[cdfX[cc]] ,nA[i],*nC,false); //init. Xij-matrix
    mui = sY[i]/2*( (Xij->cols(0,*nC-2))*(*omega) + Xij->col(*nC-1) ); //mean is vector
    ri = log(*Yi) - log(mui); //residual is vector
    sigmavec = sqrt(tau/mui); //variance-vector depending on genotype
    tmp = ri/sigmavec; //temporary variable
    lik = exp(-0.5*dot(tmp,tmp))/prod(sigmavec);
    bigsum = bigsum + lik*pG[CnQ[i]+j]; //# add likval
    cc = cc + 1; //update counter
    delete Xij;
   }//end for each j: genotype combination
   *PE = (*PE)*bigsum/pow(konstant,nA[i])/yprod; //add scaling-constant
   delete Yi;
  }//end for each loci i:
 } //end for model

 //COWELL: MAIES: DIRICHLET
 else if(*model==4) { //r ~ Dirichlet(tau*mui)
  Col<double> alphavec; //column vector used for variances
  konstant = std::tr1::lgamma(tau); //constant used in dirichlet
  cc = 0; //counter for X-matrix and O-matrix
  for(i=0; i<*nL; i++) {
   Yi =  new Col<double>( &allY[CnA[i]], nA[i], false); //needs to be scaled
   bigsum = 0.0; //sum over all genotypes
   for(j=0; j<nQi[i]; j++) { //for all genotypes: Init space and insert relevant matrice directly
    Xij = new Mat<double>( &allX[cdfX[cc]] ,nA[i],*nC,false); //init. Xij-matrix
    mui = ( (Xij->cols(0,*nC-2))*(*omega) + Xij->col(*nC-1) )/2; //mean is vector
    alphavec = tau*mui; //just multiply
    lik = konstant; //multiplicator for each j
    for(k=0;k<nA[i];k++) {
     lik = lik + (alphavec.at(k)-1)*log(Yi->at(k))-std::tr1::lgamma(alphavec.at(k)); //log-likelihood for responses
    }
    bigsum = bigsum + exp(lik)*pG[CnQ[i]+j]; //# add likval
    cc = cc + 1; //update counter
    delete Xij;
   }//end for each j: genotype combination
   *PE = (*PE)*bigsum;
   delete Yi;
  }//end for each loci i:
 } //end for model

 //mixsepPos
 else if(*model==5) { //r ~ Norm+(mui,tau)
  double sigma,sigmasq; //square root of variance parameter
  cc = 0; //counter for X-matrix and O-matrix
  for(i=0; i<*nL; i++) {
   sigmasq = tau*sY[i];
   sigma = sqrt(sigmasq);
   Yi =  new Col<double>( &allY[CnA[i]], nA[i], false); //init datavector
   bigsum = 0.0; //sum over all genotypes
   for(j=0; j<nQi[i]; j++) { //for all genotypes: Init space and insert relevant matrice directly
    Xij = new Mat<double>( &allX[cdfX[cc]] ,nA[i],*nC,false); //init. Xij-matrix
    mui = sY[i]/2*((Xij->cols(0,*nC-2))*(*omega) + Xij->col(*nC-1)); //mean is a vector
    ri = *Yi - mui; //residual is a vector
    Di = dot(ri,ri)/sigmasq; //FASTEST!!
    lik = exp(-0.5*Di)/pow(konstant*sigma,nA[i]); //log-likelihood for responses
    Di = 1;
    for(k=0;k<nA[i];k++) {
     Di = Di*std::tr1::erfc(-mui.at(k)/(sqrt(2)*sigma))/2; //scale for 0-truncation
    }
    bigsum = bigsum + lik/Di*pG[CnQ[i]+j]; //# add likval
    cc = cc + 1; //update counter
    delete Xij;
   }//end for each j: genotype combination
   *PE = (*PE)*bigsum;
   delete Yi;
  }//end for each loci i:
 } //end models
 //MAIES version normal distributed (require scaling)
// if(*model==7) { //selecting model strmix0
// double eta;
//  eta = theta[*nC-2]; //next last element is the second variance parameter
//  Col<double> sigmavec; //column vector used for variances
//  cc = 0; //counter for X-matrix and O-matrix
//  for(i=0; i<*nL; i++) {
//   Yi =  new Col<double>( &allY[CnA[i]], nA[i], false); //should be scaled
//   bigsum = 0.0; //sum over all genotypes
//   for(j=0; j<nQi[i]; j++) { //for all genotypes: Init space and insert relevant matrice directly
//    Xij = new Mat<double>( &allX[cdfX[cc]] ,nA[i],*nC,false); //init. Xij-matrix
//    mui = sY[i]/2*( (Xij->cols(0,*nC-2))*(*omega) + Xij->col(*nC-1) ); //mean is vector
//    ri = *Yi - mui; //r esidual is vector
//    sigmavec = sqrt(eta*mui + tau); //variance-vector depending on genotype
//    tmp = ri/sigmavec; //temporary variable
//    lik = exp(-0.5*dot(tmp,tmp))/prod(sigmavec);
//    bigsum = bigsum + lik*pG[CnQ[i]+j]; //# add likval
//    cc = cc + 1; //update counter
//    delete Xij;
//   }//end for each j: genotype combination
//   *PE = (*PE)*bigsum/pow(konstant,nA[i]); //add scaling-constant
//   delete Yi;
//  }//end for each loci i:
// } //end for model


//NB: NOT WORKING FOR >2 PERSONS
 //Relaxion model which increases variance for small peak heights.
// if(*model==6) { //selecting model strmix0
// int k,l;
//  double q; //power of relaxion
//  double rri; //residual one element
//  double sumtmp; //temperary sums
//  double sigmasq; //column vector used for variances
//  Row<double> muij; //expected heights for each contributor
//  Row<double> mvec; //vector of mixture proportions
//  mvec = Row<double>(*nC); //initialize vector 
//  mvec.subvec(0,*nC-2) = *omega;
//  mvec.at(*nC-1) = 1-sum(*omega);
//  q = theta[*nC-2]; //next last element is power of variance
//  cc = 0; //counter for X-matrix and O-matrix
//  for(i=0; i<*nL; i++) {
//   Yi =  new Col<double>( &allY[CnA[i]], nA[i], false); //should be scaled
//   bigsum = 0.0; //sum over all genotypes
//   for(j=0; j<nQi[i]; j++) { //for all genotypes: Init space and insert relevant matrice directly
//    Xij = new Mat<double>( &allX[cdfX[cc]] ,nA[i],*nC,false); //init. Xij-matrix
//    lik = 1;
//    for(k=0; k<nA[i]; k++) { //for all alllees
//     muij = sY[i]/2*( (Xij->row(k))%mvec ); //mean for each contributor
//     rri = Yi->at(k) - sum(muij); //r esidual is vector
//     sumtmp = 0;
//     for(l=0; l<(*nC); l++) { //for each contributor
//      if(muij.at(l)>0) { //only if contributed!
//       sumtmp = sumtmp + 1/pow(muij.at(l),q); 
//      }
//     }
//     sigmasq = tau*sumtmp; //variance in parameter
//     lik = lik*exp(-0.5*rri*rri/sigmasq)/sqrt(sigmasq); //each element are multiplied
//    } //end k: for each alleles
//    bigsum = bigsum + lik*pG[CnQ[i]+j]; //# add likval
//    cc = cc + 1; //update counter
//    delete Xij;
//   }//end for each j: genotype combination
//   *PE = (*PE)*bigsum/pow(konstant,nA[i]); //add scaling-constant
//   delete Yi;
//  }//end for each loci i:
// } //end for model

 delete omega;
}


} 
