//This script contains functions for calculating several conditional likelihoods.
//calculating likelihood of data with genotypes marginalised out
//In Makevars: Change compiler to -std=c++11
//cd ~/Dropbox/Forensic/MixtureProj/myDev/quantLR/bayesdnamix0
//R CMD SHLIB loglikgamma.cpp -O2 -larmadillo -llapack -lblas
//Other:
//R CMD SHLIB loglikgamma.cpp -llapack -lblas

#include <cmath> //includes special functions
#include <stdio.h> //indludes printf-function
#include <RcppArmadillo.h> //require RcppArmadillopackage and Namespaced defined
//#include <armadillo> 
#include <boost/math/special_functions/gamma.hpp> 
using namespace std;
using namespace arma;
using namespace boost::math;
const double PIVAL = std::acos(0.0)*2;
//#define DEBUG

 //About: We let allele-names be integers from 0:(N-1) where N are possible alleles in population (in popfreq). 
 //This correspond to index of placement in vectors. Makes assignment of X,Z,Y, very fast.
 //Y will be a N height vector (static), X is a (NxC) contribution matrix (dynamic), Z is a N vector which colSums X (this is required to quickly know what alleles are relevant in calculation).

//dropout is modeled through model (depending only on threshold)
//dropin is modeled with prior probability pC and dropped in peak heights are shifted exponential modeled with lambda-param:  


class recurseClassStutter { //recurse-class for each loci
 public: 
  //input:
  int *nC; //number of contr
  int *nGi; //Number of genotypes in populaiton
  int *nAi; //number of alleles in locus
  int *nAalli; //Number of alleles in population
  Col<uword> *Ai; //alleles
  Col<double> *Hi; //heights (small vector)
  Col<double> *Yi; //heights (full vector)
  Mat<uword> *Gset; //(Nx2)-Matrix of allele-genotypes (full vector)
  Col<double> *pGvec; //probability of genotypes
  Col<double> *pAvec; //probability of alleles
  double *prC; //probability of dropin for each allele
  double *t0; //threshold for detecting peak height
  Col<uword> *Abpind; //Indices of alleles
  double *lambda; //parameter for dropin-peak height exponential model

  //fst-correction:
  double *fst; //theta-correction
  Col<int> *mkvec; //count-vector for each alleles
  int *nkval; //counts number of sampled alleles
  Col <double> *pGenotmp; //pointer for storing genotypeproduct if fst>0

  Col<double> *mvec;
  Mat<int> *Xij; //(N x C) - matrix. Taking out columns to use for restriction
  Col<int> *Zi; //N-vector. used as colsum of Xij
  Row<int> condVec; //Define row in condmatrix
  Mat<int> XijSub; //(tilde(N) x C) - matrix. Subset of X for non-dropout: X[E,]
  Mat<int> Xijtmp; //used for stutter
  Col<double> Ytmp; //temporary used

  //parameters
  double rho;  //1. parameter part of theta
  double tau;  //2. parameter part of theta
  double xi;  //3. parameter part of theta
  double konstant;
  double konstant2;

  //temporary variables:
  Col<double> mui;
  Col<double> mutmp;
  Col<double> ri; //column vector
  double Di; //residual sum squared
  double lik;
  double pGprod; //genotype probability product
  double pDprod; //dropout/in/nondrop probability product
  double pAprod; //allele probability product (used for dropin)
  double bigsum; //summed marginal
  Col<uword> psiR; //row-indice of Y,X in calculation: contributing alleles (finds where Z>0)
  Col<uword> psiDO; //row-indice of Y,X in calculation: dropped-out alleles (finds where Y=0 AND Z>0)
  Col<uword> psiDI; //row-indice of Y,X in calculation: dropped-in alleles (finds where Y>0 AND mu=0)
  Col<uword> psiYmu; //row-indice of Y,X in calculation: non-dropped out alleles (finds where Y>0 AND mu>0)
  Col<uword> psiBP; //indices relevant to stutter (bp>0)
  Col<uword> psiRtoS; //row-indice of Y,X in calculation: indices beeing Stutter to
  Col<uword> psiRfromS; //row-indice of Y,X in calculation: indices beeing Stutter from
  Col<uword> psitmp; //temporary variable
 
  void recUpdateXM1(int k) { 
   int j,l; //used to traverse genotype probabilites
   uword l2; //used when comparing psi_elem
   for(j=0; j<*nGi; j++) { //for each loci we need to run recursion on each combination
    if( (condVec.at(k))>=0 ) { j = condVec.at(k); } //Known genotype: insert static combination 
    Zi->elem((Gset->row(j))) += 1; //update values in Zi. Updates twice if homozygote.
    if(k==(*nC-1) && *prC==0 && xi==0) { //if in last contributor, dropin-probability and stutterratio is zero -> no eval
     if(sum((Zi->elem(*Ai))==0) >0) { //still missing allele contributor in evidence: not valid combinations 
      Zi->elem((Gset->row(j))) -= 1; //update values in Zi
      continue; //continoue loop with j++
     }
    } 
    Xij->elem((Gset->row(j) + k*(*nAalli))) += 1; //update elements in Xij-matrix
    if( (condVec.at(k))<0 ) { //only if no restriction
     if(*fst>0) { //if theta-correction given
      pGenotmp->at(k) = 1; //important to reset element first!
      for(l=0; l<2; l++) { //for each allele
       pGenotmp->at(k) *= ( mkvec->at(Gset->at(j,l))*(*fst) + (1-*fst)*pAvec->at(Gset->at(j,l)) ) / ( 1+(*nkval-1)*(*fst));
       mkvec->at(Gset->at(j,l))+=1; //update counter of sampled allele
       (*nkval)+=1; //update with sampled 
      }
      if(Gset->at(j,0)!=Gset->at(j,1)) pGenotmp->at(k)*=2; 
      pGprod *=  pGenotmp->at(k);
     } else { //use genotypes directly (Hardy Weinberg)
      pGprod *= (pGvec->at(j));   //If unknown: multiply with unknown genotype product
     }
    }
    if(k==(*nC-1)) { //IF IN LAST CONTRIBUTOR
     pDprod = 1;
     if(xi>0) { //incorporate stutter ratio
      psiR = find(*Zi>0); //Indices for contributing alleles. checks those non-zero
      psitmp = Abpind->elem(psiR); //take out bp of contributing alleles.
      psiBP = find(psitmp>0); //indices relevant to stutter
      if(psiBP.n_elem>0) {
       psiRtoS = psitmp.elem(psiBP)-1; //find indices which may be stuttered to. Adjust index (starting from 0)
       psiRfromS = psiR.elem(psiBP); //find indices which may be stuttered from
       //fix stutter-contributor-X:
       Xijtmp = *Xij;
       Xijtmp.zeros();
       Xijtmp.rows(psiRtoS) = Xij->rows(psiRfromS); //insert relevant stutter-contributors
       mui = ( (1-xi)*((*Xij)*(*mvec)) + xi*(Xijtmp*(*mvec))  ); //mean peak height of model (contributed means only)
      } else { //if none in stutter-position
       mui = ( (1-xi)*((*Xij)*(*mvec)) ); //mean peak height of model (contributed means only)
      }
     } else { //no stutter ratio
      mui = ( (*Xij)*(*mvec) ); //mean peak height of model
     }
     mui = mui*rho; //this is now alpha-parameter!
     psiYmu = find( ((*Yi>0) + (mui>0))==2 ); //Indices for modelled alleles
     psiDO = find( ((*Yi==0) + (mui>0))==2 ); //Indices for dropped out alleles
     psiDI = find( ((*Yi>0) + (mui==0))==2 ); //Indices for dropped in alleles

     //calculate dropin:
     if(*prC>0) { //only if drop-in probability is >0. 
      if(psiDI.n_elem>0) {        
       pDprod = prod( (*prC)*(pAvec->elem(psiDI))) ; //multiply with allele probability
       if(*lambda>0) { //weighting drop-out probability with peak height. Lambda=0 reduces to standard procedure
        pDprod *= exp( psiDI.n_elem*( log(*lambda) + (*lambda)*(*t0) ) - (*lambda)*sum(Yi->elem(psiDI)) );
       }
      } else { //if no dropin found
       pDprod = (1-*prC); //scale with probability of not dropping in
      }
     } //end if drop-in probability given
     if(*prC>0 || psiDI.n_elem==0) { //if possible for dropin or no drop-in given
      //Step 2) calculating log-likelihood of models
      lik = 0; //summarize log-likelihood
       Ytmp = Yi->elem(psiYmu); //take out relevant peak heights
       mutmp = mui.elem(psiYmu); //this is alpha-parameter!
       lik = dot(mutmp,log(Ytmp))- konstant*sum(mutmp)  - sum(log(Ytmp)) - sum(Ytmp)*konstant2; 
       for(l2=0;l2<psiYmu.n_elem;l2++) { //for each alleles in non-dropped out allees
        lik = lik - std::tr1::lgamma(mui.at(psiYmu.at(l2))); //add last expression in sum 
       }
       if(psiDO.n_elem>0) { //there are drop.out elements (i.e. contributing genos gives peak 0)
        for(l2=0;l2<psiDO.n_elem;l2++) { //for each dropped out alleles (erf only takes elements)
         Di = mui.at(psiDO.at(l2)); //this is alpha-parameter!
         Di = gamma_p(Di,(*t0)*konstant2); //see formula for cdf-gamma
         lik = lik + log(Di);// - std::lgamma(mui.at(psiDO.at(l))*konstant2); //add log-probability
        }
       } //end dropout
      bigsum += exp(lik)*pGprod*pDprod; //# multiply with dropout-prob and genotype prob
     } //end possible combination
    } else { //IF NOT IN LAST CONTRIBUTOR
     recUpdateXM1(k+1); //recurse to next contributor if not in last
    }
    //reversing:
    Xij->elem((Gset->row(j) + k*(*nAalli))) -= 1; //update elements in Xij-matrix
    Zi->elem((Gset->row(j))) -= 1; //update values in Zi. Updates twice if homozygote.
    if(condVec.at(k)>=0) { 
     j = *nGi;  //end loop if known combination
    } else {
     if(*fst>0) { //if theta-correction given
      pGprod /=  pGenotmp->at(k);
      mkvec->elem(Gset->row(j))-=1; //update counter of sampled allele
      (*nkval)-=2; //update with sampled 
     } else { //use genotypes directly (Hardy Weinberg)
      pGprod /= (pGvec->at(j));   //If unknown: multiply with unknown genotype product
     }
    }
   } //end for each combination j
   //return;
  } //end recursive function

  recurseClassStutter(double *pC, double *pG, double *pA, Row<int> condV, uword *A, double *H, uword *Gvec, int *C, int *nA, int *nAall, int *nG, Col<double> *omega, uword *allAbpind,double *th ,double *t0in, double *fstin, int *mkvector, int *nkvalue, double *lam) {
   nC = C; //copy pointer to number of contributors
   nAi = nA; //copy pointer to number of alleles in evidence
   nAalli = nAall; //copy pointer to number of alleles in population
   nGi = nG; //copy pointer to number of genotypes
   mvec = omega;  //copy pointer to mix-prop vector
   condVec = condV; //assign value where condV points 
   prC = pC; //copy pointer to drop-in probability
   fst = fstin; //copy pointer theta-correction
   nkval = nkvalue; //copy pointer to sample-counter
   t0 = t0in;//copy pointer to threshold
   lambda = lam; //copy pointer to parameter of exponential drop-in model
//reparameterization: (mu,sd) -> (rho,tau)
   tau = th[1]*th[1]/th[0];  
   rho = th[0]/tau;  //get value of model parameter (distr-param)
//   rho = th[0];  //get value of model parameter (distr-param)
//   tau = th[1];  //get value of model parameter (distr-param)
   xi = th[2]; //get value of model parameter (stutter-param)
   konstant = log(tau);
   konstant2 = 1/tau;

   Ai = new Col<uword>( A, *nAi, false); //insert allele evidence (allele indices)
   Hi = new Col<double>( H, *nAi, false); //insert allele evidence
   Gset = new Mat<uword>( Gvec, *nGi, 2,false); //insert genotype-combinations (allele indices)
   pGvec = new Col<double>( pG, *nGi, false); //insert genotype probabilities
   pAvec = new Col<double>( pA, *nAalli, false); //insert allele probabilities
   mkvec = new Col<int>( mkvector, *nAalli, false); //insert allele-sampled
   Abpind = new Col<uword>( allAbpind, *nAalli, false); //insert allele evidence (allele indices)

   pGenotmp = new Col<double>(*nC); //init datavector
   Yi =  new Col<double>(*nAalli); //init datavector
   Zi =  new Col<int>(*nAalli); //init datavector
   Xij = new Mat<int>(*nAalli, *nC); //init. Xij-matrix with zeros
   Yi->zeros(); //insert all as zeros
   Xij->zeros();
   Zi->zeros();
   Yi->elem(*Ai) = *Hi; //insert observed peak heights

//   Yi->submat(E) = sub( allY[CnA[i]], nA[i]); //insert observed peak heights
   bigsum = 0.0; //init big sum over all genotypes
   pGprod = 1.0; //init genotype-prob product
   //start recursion
   recUpdateXM1(0); //recursion over all genotypes within locus
   //delete when finished
   delete Ai;
   delete Hi;
   delete Gset;
   delete pGvec;
   delete pAvec;
   delete mkvec;

   delete Abpind;
   delete pGenotmp;
   delete Yi;
   delete Zi;
   delete Xij;
  } //end constructor
}; //end recursiveClassStutter


extern "C" {

//function for calculating likelihood of data with genotypes marginalised out
void loglikgammaC(double *logPE, double *theta, int *nC, int *nL, int *nA, double *allY, uword *allA , int *CnA,uword *allAbpind, int *nAall, int *CnAall, uword *Gvec, int *nG, int *CnG,int *CnG2, double *pG,  double *pA,  double *pC, int* condRef, double *t0, double *fst, int *mkvec, int *nkval, double *lambda) {

 //logPE  - logged density value of evidence
 //theta  - model parameters
 //nC  - number of contributors
 //nL  - number of loci
 //nA - number of alleles on each loci (in evidence)
 //allY - observed peak height of alleles
 //allA - Observed alleles
 //CnA - Cumulative number of alleles on each loci (in evidence)
 //allAbpind - index of corresponding 1-bp ahead twin of allele (index start from '1', '0' means no index)
 //nAall - number of alleles on each loci (all in population)
 //CnAall - Cumulative number of alleles on each loci (all in population)
 //Gvec - matrix of genotypes
 //nG  - number of genotypes on each loci
 //CnG  - cumulative position of genotype-vector
 //CnG2  - cumulative position of genotype-matrix
 //pG  -  genotype probability vector
 //pA  -  allele-probability vector
 //pC  -  drop-in probability
 //pG  -  genotype-probability 
 //condRef - matrix of values of conditioned genotypes (-1: no restriction.)
 //t0 - threshold used in drop-out
 //fst - theta-correction 
 //mkvec - full vector of #sampled alleles (each alleles in population) 
 //nkval - total number of sampled alleles (for each loci)
 int i;
 double Li; //likelihood for a given locus
 Col<double> *omega; //mixture proportion 
 Mat<int> *condMatrix; //conditional matrix for each contributor (values equal Gset-indices)
 omega = new Col<double>( theta, *nC-1, false); //insert mixture proportion: (C-1)-vector
 condMatrix = new Mat<int>(condRef, *nL, *nC,false); //insert condRef-matrix directly
 Col<double> mvec; //vector of mixture proportions. Update in X easier
 mvec = Col<double>(*nC); //initialize vector 
 if(*nC==1) {
  mvec.at(0) = 1; //insert full mix-prop if one contributor 
 } else {
  mvec.subvec(0,*nC-2) = *omega;
  mvec.at(*nC-1) = 1-sum(*omega);  //restrict last mix-prop as sum of the others
 }

 //for each loci we need to run the recursion:
 for(i=0; i<*nL; i++) { 
   recurseClassStutter *rec = new recurseClassStutter(pC, &pG[CnG[i]],&pA[CnAall[i]] ,condMatrix->row(i),&allA[CnA[i]],&allY[CnA[i]], &Gvec[CnG2[i]], nC, &nA[i], &nAall[i],&nG[i],  &mvec, &allAbpind[CnAall[i]], &theta[*nC-1]  , t0, fst, &mkvec[CnAall[i]], &nkval[i], lambda); //create object of recursion
   Li = rec->bigsum; //extract likelihood
   delete rec; //delete recurse object
   *logPE = (*logPE) + log(Li);
   if(Li==0) {
     break; //finished if the likelihood hits 0
   }//end for each loci i:
 }  
 delete omega;
 delete condMatrix;
} //end function

} 
