//This script contains functions for calculating several conditional likelihoods.
//calculating likelihood of data with genotypes marginalised out
//cd ~/Dropbox/Forensic/MixtureProj/myDev/
//R CMD SHLIB contLikDrop.cpp -llapack -lblas

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

 //allY - vector of peak heights
 //allA - vector of allele (must be of type uvec?)
 //CnA - Cumulative number of alleles in each loci
 //sY - sum peak height at each loci

 //Gset - integers (uvec?) of alleles. Vectorized matrix (sum(Ni) x 2). 
 //nG - number of genotypes in each loci
 //CnG - Cumulative number of genotypes in each loci
 //pG - vector of genotype-probabilities in Gset
 //pDvec - matrix of drop-out probability for each genotypes in Gset (row). Each contributors (col) may have different drop-out prob.
 //condRef - integers for conditioned on Gset-indizes for each contributors. Vectorized matrix (nL x nC)
 //t0 - peak height imputed for dropped alleles.

 //We let allele-names be integers from 0:(N-1) where N are possible alleles in population (in popfreq). 
 //This correspond to index of placement in vectors. Makes assignment of X,Z,Y, very fast.
 //Y will be a N height vector (static), X is a (NxC) contribution matrix (dynamic), Z is a N vector which colSums X (this is required to quickly know what alleles are relevant in calculation).

class recurseClass { //recurse-class for each loci
 public: 
  //input:
  int *model; //number of contr
  int *nC; //number of contr
  int *nGi; //Number of genotypes in populaiton
  int *nAi; //number of alleles in locus
  int *nAalli; //Number of alleles in population
  Col<uword> *Ai; //alleles
  Col<double> *Hi; //heights (small vector)
  Col<double> *Yi; //heights (full vector)
  Col<double> Ytmp; //temporary used
  Mat<uword> *Gset; //(Nx2)-Matrix of allele-genotypes (full vector)
  Col<double> *pGvec; //probability of genotypes
  Col<double> *pAvec; //probability of alleles
  double *pDvec; //probability of dropout for each contributor
  double *prC; //probability of dropin for each allele

  Col<double> *mvec;
  Mat<int> *Xij; //(N x C) - matrix. Taking out columns to use for restriction
  Col<int> *Zi; //N-vector. used as colsum of Xij
  Row<int> condVec; //Define row in condmatrix
  Mat<int> XijSub; //(tilde(N) x C) - matrix. Subset of X for non-dropout: X[E,]


  //parameters
  double sigma;
  double sigmasq;
  double konstant;
  double sYi;

  //temporary variables:
  Col<double> mui;
  Col<double> ri; //column vector
  double Di; //residual sum squared
  double lik;
  double pGprod; //genotype probability product
  double pDprod; //dropout/nondropout probability product
  double pDprodTmp; //dropout/nondropout probability product
  double pAprod; //allele probability product (used for dropin)
  double bigsum; //summed marginal
  Col<uword> psiR; //row-indice of Y,X in calculation: contributing alleles

  void recUpdateXM1(int k) { 
   int j; //used to traverse genotype probabilites
   for(j=0; j<*nGi; j++) { //for each loci we need to run recursion on each combination
    if( (condVec.at(k))>=0 ) { j = condVec.at(k); } //Known genotype: insert static combination 
    Zi->elem((Gset->row(j))) += 1; //update values in Zi. Updates twice if homozygote.
    if(k==(*nC-1) && *prC==0) { //if in last contributor and dropin-probability is zero -> no necessary evaluation
     if(sum((Zi->elem(*Ai))==0) >0) { //still missing allele contributor in evidence: not valid combinations 
      Zi->elem((Gset->row(j))) -= 1; //update values in Zi
      continue; //continoue loop with j++
     }
    } 
    Xij->elem((Gset->row(j) + k*(*nAalli))) += 1; //update elements in Xij-matrix
    if( (condVec.at(k))<0 ) { pGprod = pGprod*(pGvec->at(j));  } //If unknown: multiply with unknown genotype product

    if(k==(*nC-1)) { //IF IN LAST CONTRIBUTOR
     int l,m; //used to traverse alleles in X
     int nondrop; //if non-droped
     psiR = find((*Zi)>0); //Indices for contributing alleles. checks those non-zero
     //Step 1) calculate dropout/dropin probabilities using info in X
     pDprod = 1; //must multiply drop-out probability for each allele (both dropped and non-dropped)
     XijSub = Xij->rows(psiR); //get relevant alleles for X
     //1a) find info of contributing alleles (B1,B2) dropout-sets
     for(l=0; l<psiR.n_elem; l++) { 
      nondrop = sum((*Ai)==psiR[l]); //check whether allele are non-dropped (0,1)
      pDprodTmp = 1;
      for(m=0; m<*nC; m++) { //for each contributor
        pDprodTmp = pDprodTmp*pow(pDvec[m],XijSub.at(l,m)); //powered with number of contributing alleles
      }  //end for each contributors
      if(nondrop) { 
       pDprod = pDprod*(1-pDprodTmp);  //if non-droped allele (get probability for no dropout)
      }  else {
       pDprod = pDprod*pDprodTmp;  //if non-droped allele (get probability for no dropout) 
      } 
     } //end for each allele-elements
     //1b) find info of non-contributing alleles (C) drop-in
     //drop-in heights are considered as artifacts (not peak height modeled)
     if(*prC>0) { //only if drop-in probability is >0. 
      pAprod = 1;
      for(l=0; l<*nAi; l++) { //check each allele
       nondrop = sum(psiR==Ai->at(l)); //check whether allele are drop-in (dropin=0,no dropin=1)
       if(!nondrop) { //if dropin
        pAprod = pAprod*(*prC)*(pAvec->at( Ai->at(l) )); //multiply with allele probability
       }
      }
      if(pAprod==1) { //if no dropin found
       pDprod = pDprod*(1-*prC); //scale with probability of not dropping in
      } else {  
       pDprod = pDprod*pAprod; //multiply with dropin-probability of droped in alleles
      }
     } //end if drop-in probability given
     //Step 2) calculating likelihood of peak heights
     lik = 1; //binary model is default
     if(*model>0) {
      mui = sYi/2*( (Xij->rows(psiR))*(*mvec) ); //mean peak height
      Ytmp = (Yi->elem(psiR));
      ri = Ytmp/sum(Ytmp) - mui; //residual is a vector
      Di = dot(ri,ri)/sigmasq; //FASTEST!! 
      lik = exp(-0.5*Di)/pow(konstant,psiR.n_elem); //likelihood of model
     }
     bigsum = bigsum + lik*pGprod*pDprod; //# multiply with dropout-prob and genotype prob
    } else { //IF NOT IN LAST CONTRIBUTOR
     recUpdateXM1(k+1); //recurse to next contributor if not in last
    }

    //reversing:
    Xij->elem((Gset->row(j) + k*(*nAalli))) -= 1; //update elements in Xij-matrix
    Zi->elem((Gset->row(j))) -= 1; //update values in Zi. Updates twice if homozygote.
    if(condVec.at(k)>=0) { 
     j = *nGi;  //end loop if known combination
    } else {
     pGprod = pGprod/(pGvec->at(j)); //reverse genotype probability if unknown combination
    }
   } //end for each combination
   return;
  } //end recursive function

  recurseClass(int *mod, double *pC, double *pD, double *pG, double *pA, Row<int> condV, uword *A, double *H, uword *Gvec, int *C, int *nA, int *nAall, int *nG, Col<double> *omega, double *sY,double *tau ,double *t0) {
   nC = C; //copy pointer to number of contributors
   nAi = nA; //copy pointer to number of alleles in evidence
   nAalli = nAall; //copy pointer to number of alleles in population
   nGi = nG; //copy pointer to number of genotypes
   mvec = omega;  //copy pointer to mix-prop vector
   sYi = *sY;
   sigmasq = *tau;  //copy pointer to mix-prop vector
   condVec = condV; //assign value where condV points 
   pDvec = pD; //copy pointer to drop-out probabilities
   prC = pC; //copy pointer to drop-in probability
   model = mod; //copy pointer to model choice

   sigma = sqrt(*tau);
   konstant = sqrt(2*PIVAL)*sigma;
   Ai = new Col<uword>( A, *nAi, false); //insert allele evidence (allele indices)
   Hi = new Col<double>( H, *nAi, false); //insert allele evidence
   Gset = new Mat<uword>( Gvec, *nGi, 2,false); //insert genotype-combinations (allele indices)
   pGvec = new Col<double>( pG, *nGi, false); //insert genotype probabilities
   pAvec = new Col<double>( pA, *nAalli, false); //insert allele probabilities
   
   Yi =  new Col<double>(*nAalli); //init datavector
   Zi =  new Col<int>(*nAalli); //init datavector
   Xij = new Mat<int>(*nAalli, *nC); //init. Xij-matrix with zeros
   Yi->fill(*t0); //insert imputed peak heights
   Xij->zeros();
   Zi->zeros();
   Yi->elem(*Ai) = *Hi; //insert heights
   bigsum = 0.0; //init big sum over all genotypes
   pGprod = 1.0; //init genotype-prob product

   //start recursion
   recUpdateXM1(0); 
   //delete when finished
   delete Xij;
   delete Yi;
   delete Zi;
  } //end constructor
}; //end recursiveClass


extern "C" {

//function for calculating likelihood of data with genotypes marginalised out
void contlikdropC(double *PE, double *theta, int *model, int *nC, int *nL, int *nA, double *allY, uword *allA , int *CnA, double *sY, int *nAall, int *CnAall, uword *Gvec, int *nG, int *CnG,int *CnG2, double *pG,  double *pA, double *pDvec,  double *pC, int* condRef, double *t0) {

 //pC  -  drop-in probability
 //pDvec  -  drop-out probability vector (each contributor)
 //pG  -  genotype-probability vector
 //pA  -  allele-probability vector
 //CnG  - cumulative position of genotype-vector
 //CnG2  - cumulative position of genotype-matrix
 //nA - number of alleles on each loci (in evidence)
 //CnA - Cumulative number of alleles on each loci (in evidence)
 //nAall - number of alleles on each loci (all in population)
 //CnAall - Cumulative number of alleles on each loci (all in population)
 //model - {0 = binary, 1=mixsep}

 int i;
 Col<double> *omega; //mixture proportion 
 Mat<int> *condMatrix; //conditional matrix for each contributor (values equal Gset-indices)
 omega = new Col<double>( theta, *nC-1, false); //insert mixture proportion: (C-1)-vector
 condMatrix = new Mat<int>(condRef, *nL, *nC,false); //insert condRef-matrix directly
 Col<double> mvec; //vector of mixture proportions. Update in X easier
 mvec = Col<double>(*nC); //initialize vector 
 mvec.subvec(0,*nC-2) = *omega;
 mvec.at(*nC-1) = 1-sum(*omega); 

 //for each loci we need to run the recursion:
 for(i=0; i<*nL; i++) { 
   recurseClass *rec = new recurseClass(model, pC, pDvec, &pG[CnG[i]],&pA[CnAall[i]] ,condMatrix->row(i),&allA[CnA[i]],&allY[CnA[i]], &Gvec[CnG2[i]], nC, &nA[i], &nAall[i],&nG[i],  &mvec, &sY[i], &theta[*nC-1]  , t0); //create object of recursion
   *PE = (*PE)*(rec->bigsum);
   delete rec;
 }//end for each loci i:

 delete omega;
 delete condMatrix;
} //end function



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