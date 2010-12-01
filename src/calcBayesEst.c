/*
#######################################################################
## This program is Open Source Software: you can redistribute it
## and/or modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see http://www.gnu.org/licenses/.
*/

#include<R.h>
#include<Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Applic.h>

/* Function for calculation of Bayes dose-response estimates */

/*
Notation a matrix is stored as a double
in order a11, a12, a21, a22, all matrices
assumed to be 2x2
*/

/* calculate determinant */
double dtm(double *X){
  return X[0]*X[3]-X[1]*X[2];
}

/* invert matrix */
void invert(double *X, double *Xinv){
  double det=0;
  det = dtm(X);
  Xinv[0] = X[3]/det;
  Xinv[3] = X[0]/det;
  Xinv[1] = -X[1]/det;
  Xinv[2] = -X[2]/det;
}

/* Calculate X'DX. where D is diagonal matrix with sample sizes
   on the diagonal and X = (1,x) */
void XDX(double *x, int *n, int length, double *res){
  int i;
  double sumn=0, sumnx=0, sumnx2=0;
  for(i=0; i<length; i++){
    sumn += n[i];
    sumnx += n[i]*x[i];
    sumnx2 += n[i]*pow(x[i], 2);
  }
  res[0] = sumn;
  res[1] = res[2] = sumnx;
  res[3] = sumnx2;
}

/* Calculate XDybar, where D diagonal */
void XDy(double *x, int *n, int length, double *ybar, double *res){
  int i;
  double sumn=0, sumnxy=0;
  for(i=0; i<length; i++){
    sumn += n[i]*ybar[i];
    sumnxy += x[i]*n[i]*ybar[i];
  }
  res[0] = sumn;
  res[1] = sumnxy;
}

/* Calculate quadratic form x'Ax */
double xAx(double *x, double *A){
  return A[0]*pow(x[0], 2)+(A[1]+A[2])*x[0]*x[1]+A[3]*pow(x[1],2);
}

/* Calculate matrix vector product */
void Ax(double *A, double *x, double *res){
  res[0] = A[0]*x[0] + A[1]*x[1];
  res[1] = A[2]*x[0] + A[3]*x[1];
}

/* Calculate matrix factorization of sym. 2x2 matrix */
void mf(double *A, double *V){
  // returns upper triangular matrix V
  // such that VV'= A
  V[3] = pow(A[3], 0.5);
  V[1] = A[1]/V[3];
  V[2] = 0;
  V[0] = pow(A[0]-pow(A[1], 2)/A[3], 0.5);
}

/* standardized model functions */
double emax(double x, double ed50){
  return x/(x+ed50);
}

double exponential(double x, double delta){
  return exp(x/delta)-1;
}

double logistic(double x, double ed50, double delta){
  return 1/(1+exp( (ed50 - x)/ delta ));
}

double betaMod(double x, double delta1, double delta2,
               double scal){
  double maxDens=0;
  double temp=0;
  maxDens = pow(delta1, delta1)*pow(delta2,delta2);
  maxDens = maxDens/(pow(delta1 + delta2,delta1 + delta2));
  temp = x/scal;
  return 1/maxDens*(pow(temp,delta1)*pow(1-temp,delta2));
}

/* evaluate target densities */
double evalDens(double *x, int *length, int *n, int ntotal, double *ybar, double yDinvy, double s2,
          double pra, double prd, double *Binv, double *m, double numScal, double *mean){
  int i;  
  double xdx[4]={0, 0, 0, 0},Vs[4]={0, 0, 0, 0}, as=0,Vsinv[4]={0, 0, 0, 0};
  double Binvm[2]={0, 0},xdy[2]={0, 0},quadF=0;
          
  /* Calculate Vs */
  XDX(x, n, *length, xdx);
  for(i=0;i<4;i++){
    Vsinv[i] = xdx[i]+Binv[i];
  }
  invert(Vsinv, Vs);
  Ax(Binv, m, Binvm);
  XDy(x, n, *length, ybar, xdy);
  Binvm[0] = xdy[0] + Binvm[0];
  Binvm[1] = xdy[1] + Binvm[1];

  /* posterior mean */
  Ax(Vs, Binvm, mean);
  quadF = mean[0]*Binvm[0]+mean[1]*Binvm[1];
  as = pra + xAx(m, Binv) + yDinvy - quadF + s2;
  return pow(dtm(Vs), 0.5)*pow(as/numScal,-(prd + ntotal)/2);
}

/* emax model */
/* prior = (a, d, m,    V,    alpha, beta)*/
/*          0  1 2,3  4,5,6,7  8,     9 */
void NumIntemax(double *nodes, int *nodlgth, double *aBound, double *bBound, 
   int *length, double *grpmean, double *yDinvy, double *s2,
   double *prior, double *doses, int *n, int *ntotal, double *numScal, double *x,
   double *normConst, double *postM, int *meanInd){
  int i,j;
  double c=0,priordens=0;
  double B[4]={0, 0, 0, 0};
  double Binv[4]={0, 0, 0, 0};
  double m[2]={0, 0}; 
  double mean[2]={0, 0}, intLik=0;
  double mult=0,bma=0,maxPost=0;
  double nod=0;
  bma = *bBound - *aBound;
  mult = (prior[1]-2)/ prior[0];
  

  for(i=0; i<*nodlgth; i++){
    /* Translate knowledge to NIG prior */
    nod = nodes[i]*bma+*aBound;
    c = emax(doses[*length-1], nod); /* assumes ordered doses */
    B[0] = mult * prior[4];
    B[1] = B[2] = mult * prior[5]/c;
    B[3] = mult * prior[7]/pow(c, 2);
    m[0] = prior[2];
    m[1] = prior[3]/c;
    invert(B, Binv);
    for(j=0; j<*length; j++){
      x[j] =  emax(doses[j], nod);
    } 
    /* eval target density */
    intLik = evalDens(x, length, n, *ntotal, grpmean, *yDinvy, *s2, prior[0], prior[1], Binv, m,
              *numScal, mean)/pow(dtm(B), 0.5);
    priordens = dbeta(nodes[i], prior[8], prior[9], 0);
    *normConst += intLik*priordens;
    nodes[i] = intLik*priordens;
    if(*meanInd){
      postM[0] += mean[0]*intLik*priordens;
      postM[1] += mean[1]*intLik*priordens;   
      postM[2] += nod*intLik*priordens;
    } else {
      if(intLik*priordens > maxPost){
        maxPost = intLik*priordens;
        postM[0] = mean[0];
        postM[1] = mean[1];   
        postM[2] = nod;        
      }
    }
  }
  if(*meanInd){
    postM[0] = postM[0]/ *normConst;
    postM[1] = postM[1]/ *normConst;
    postM[2] = postM[2]/ *normConst;
  }
  *normConst = *normConst/ *nodlgth;
}

/* linear model */
/* prior = a, d,  m,     V */
/*         0  1  2,3  4,5,6,7*/ 
void NumIntlinear(int *length, double *grpmean, double *yDinvy, double *s2,
   double *prior, double *doses, int *n, int *ntotal, double *numScal, double *x, 
   double *normConst, double *postM){
  int j;
  double c=0;
  double B[4]={0, 0, 0, 0};
  double Binv[4]={0, 0, 0, 0};
  double m[2]={0, 0}; 
  double mean[2]={0, 0}, intLik=0, mult=0;
  mult = (prior[1]-2)/ prior[0];

  c = doses[*length-1]; /* assumes ordered doses */
  B[0] = mult* prior[4];
  B[1] = B[2] = mult*prior[5]/c;
  B[3] = mult*prior[7]/pow(c, 2);
  m[0] = prior[2];
  m[1] = prior[3]/c;
  invert(B, Binv);
  for(j=0; j<*length; j++){
    x[j] =  doses[j];
  } 
  /* eval target density */
  intLik = evalDens(x, length, n, *ntotal, grpmean, *yDinvy, *s2, prior[0], prior[1], Binv, m,
            *numScal, mean)/pow(dtm(B), 0.5);

  *normConst = intLik;
  postM[0] = mean[0];
  postM[1] = mean[1];   
}

/* linlog model */
/* prior = a, d,  m,     V */
/*         0  1  2,3  4,5,6,7*/ 
void NumIntlinlog(int *length, double *grpmean, double *yDinvy, double *s2,
   double *prior, double *doses, int *n, int *ntotal, double *numScal, 
   double *x, double *normConst, double *postM){
  /* numScal contains numScal and off parameter*/
  int j;
  double c=0;
  double B[4]={0, 0, 0, 0};
  double Binv[4]={0, 0, 0, 0};
  double m[2]={0, 0}; 
  double mean[2]={0, 0}, intLik=0, mult=0;
  mult = (prior[1]-2)/ prior[0];

  c = log(doses[*length-1] + numScal[1])-log(numScal[1]); /* assumes ordered doses */
  B[0] = mult* (prior[4] -2*log(numScal[1])/c*prior[5]+pow(log(numScal[1])/c,2)*prior[7]);
  B[1] = B[2] = mult*(prior[5]/c - log(numScal[1])/pow(c,2)*prior[7]);
  B[3] = mult*prior[7]/pow(c, 2);
  m[0] = prior[2]-prior[3]*log(numScal[1])/c;
  m[1] = prior[3]/c;
  invert(B, Binv);
  for(j=0; j<*length; j++){
    x[j] =  log(doses[j] + numScal[1]);
  } 
  /* eval target density */
  intLik = evalDens(x, length, n, *ntotal, grpmean, *yDinvy, *s2, prior[0], prior[1], Binv, m,
            numScal[0], mean)/pow(dtm(B), 0.5);

  *normConst = intLik;
  postM[0] = mean[0];
  postM[1] = mean[1];   
}

/* Beta model */
/* prior a, d,  m,     V,   alpha1,beta1  alpha2,beta2*/
/*       0  1  2,3  4,5,6,7   8,9   10,11  */
void NumIntbetaMod(double *nodes1, double *nodes2, int *nodlgth, double *aBound, double *bBound, int *length, double *grpmean,
   double *yDinvy, double *s2, double *prior, double *doses, int *n, int *ntotal, 
   double *numScal, double *x, double *normConst, double *postM, int *meanInd){
  /* numScal contains numScal and usual scal parameter */
  int i,k;
  double priordens=0;
  double B[4]={0, 0, 0, 0};
  double Binv[4]={0, 0, 0, 0};
  double m[2]={0, 0}; 
  double mean[2]={0, 0}, intLik=0;
  double bma1=0, bma2=0, mult=0, nodi=0, nodj=0, c=0, maxPost=0;

  bma1 = (bBound[0]-aBound[0]);
  bma2 = (bBound[1]-aBound[1]);
  mult = (prior[1]-2)/ prior[0];

  for(i=0; i<*nodlgth; i++){
    nodi = nodes1[i]*bma1+aBound[0];
    nodj = nodes2[i]*bma2+aBound[1];
    /* Translate knowledge to NIG prior */
    if(nodi/(nodi+nodj) > doses[*length-1]/ numScal[1]){
      c = betaMod(doses[*length-1], nodi, nodj, numScal[1]);
      B[0] = mult * prior[4];       
      B[1] = B[2] = mult*prior[5]/c;
      B[3] = mult*prior[7]/pow(c, 2);
      m[0] = prior[2];
      m[1] = prior[3]/c;
    } else {
      B[0] = mult * prior[4];
      B[1] = B[2] = mult*prior[5];
      B[3] = mult*prior[7];
      m[0] = prior[2];
      m[1] = prior[3];
    }  
    invert(B, Binv);
    
    for(k=0; k<*length; k++){
      x[k] = betaMod(doses[k], nodi, nodj, numScal[1]);
    } 
   /* eval target density */
    intLik = evalDens(x, length, n, *ntotal, grpmean, *yDinvy, *s2, prior[0], prior[1], Binv, m,
              numScal[0], mean)/pow(dtm(B),0.5);
    priordens = dbeta(nodes1[i], prior[8], prior[9], 0);
    priordens = priordens*dbeta(nodes2[i], prior[10], prior[11], 0);
    *normConst += intLik*priordens;
    if(*meanInd){
      postM[0] += mean[0]*intLik*priordens;
      postM[1] += mean[1]*intLik*priordens;
      postM[2] += nodi*intLik*priordens;
      postM[3] += nodj*intLik*priordens;
    } else {
      if(intLik*priordens > maxPost){
        maxPost = intLik*priordens;
        postM[0] = mean[0];
        postM[1] = mean[1];
        postM[2] = nodi;
        postM[3] = nodj;
      }
    }
  }
  if(*meanInd){
    postM[0] = postM[0]/ *normConst;
    postM[1] = postM[1]/ *normConst;
    postM[2] = postM[2]/ *normConst;
    postM[3] = postM[3]/ *normConst;
  }
  *normConst = *normConst/ *nodlgth;
}

/* Logistic model */
/* prior a, d,  m,     V, alpha1,beta1  alpha2,beta2 */
/*       0  1  2,3  4,5,6,7   8,9   10,11  */
void NumIntlogistic(double *nodes1, double *nodes2, int *nodlgth, double *aBound, double *bBound,
   int *length, double *grpmean, double *yDinvy, double *s2, double *prior, double *doses, int *n, 
   int *ntotal, double *numScal, double *x, double *normConst, double *postM, int *meanInd){
  int i,k;
  double priordens=0;
  double B[4]={0, 0, 0, 0};
  double Binv[4]={0, 0, 0, 0};
  double m[2]={0, 0}; 
  double mean[2]={0, 0}, intLik=0;
  double bma1 = 0, bma2 = 0, mult=0, nodi=0, nodj=0;
  double c1=0,c2=0,maxPost=0;

  bma1 = bBound[0]-aBound[0];
  bma2 = bBound[1]-aBound[1];
  mult = (prior[1]-2)/ prior[0];

  for(i=0; i< *nodlgth; i++){
    nodi = nodes1[i]*bma1+aBound[0];
    nodj = nodes2[i]*bma2+aBound[1];
    /* Translate knowledge to NIG prior */
    c1 = logistic(doses[*length-1], nodi, nodj)-logistic(0, nodi, nodj);
    c2 = logistic(0, nodi, nodj)/c1;
    B[0] = mult * (prior[4] - 2*c2*prior[5] + c2*c2*prior[7]);
    B[1] = B[2] = mult*(prior[5]/c1 - c2/c1*prior[7]);
    B[3] = mult*prior[7]/(c1*c1);
    m[0] = prior[2]-prior[3]*c2;
    m[1] = prior[3]/c1;
    invert(B, Binv);      
    for(k=0; k<*length; k++){
      x[k] = logistic(doses[k], nodi, nodj);
    } 
    /* eval target density */
    intLik = evalDens(x, length, n, *ntotal, grpmean, *yDinvy, *s2, prior[0], prior[1], Binv, m,
          *numScal, mean)/pow(dtm(B), 0.5);
    priordens = dbeta(nodes1[i], prior[8], prior[9], 0);
    priordens = priordens*dbeta(nodes2[i], prior[10], prior[11], 0);
    *normConst += intLik*priordens;
    if(*meanInd){
      postM[0] += mean[0]*intLik*priordens;
      postM[1] += mean[1]*intLik*priordens;
      postM[2] += nodi*intLik*priordens;
      postM[3] += nodj*intLik*priordens;
    } else {
      if(intLik*priordens > maxPost){
        maxPost = intLik*priordens;
        postM[0] = mean[0];
        postM[1] = mean[1];
        postM[2] = nodi;
        postM[3] = nodj;
      }
    }
  }
  if(*meanInd){
    postM[0] = postM[0]/ *normConst;
    postM[1] = postM[1]/ *normConst;
    postM[2] = postM[2]/ *normConst;
    postM[3] = postM[3]/ *normConst;
  }
  *normConst = *normConst/ *nodlgth;
}

/* exponential model */
/* prior = (a, d, m,    V,    alpha, beta)*/
/*          0  1 2,3  4,5,6,7  8,     9 */
void NumIntexponential(double *nodes, int *nodlgth, double *aBound, double *bBound, 
   int *length, double *grpmean, double *yDinvy, double *s2,
   double *prior, double *doses, int *n, int *ntotal, double *numScal, double *x,
   double *normConst, double *postM, int *meanInd){
  int i,j;
  double c=0,priordens=0;
  double B[4]={0, 0, 0, 0};
  double Binv[4]={0, 0, 0, 0};
  double m[2]={0, 0}; 
  double mean[2]={0, 0}, intLik=0;
  double bma = 0, mult=0;
  double nod=0,maxPost=0;
  bma = *bBound-*aBound;
  mult = (prior[1]-2)/ prior[0];
  

  for(i=0; i<*nodlgth; i++){
    /* Translate knowledge to NIG prior */
    nod = nodes[i]*bma+*aBound;
    c = exponential(doses[*length-1], nod); /* assumes ordered doses */
    B[0] = mult * prior[4];
    B[1] = B[2] = mult * prior[5]/c;
    B[3] = mult * prior[7]/pow(c, 2);
    m[0] = prior[2];
    m[1] = prior[3]/c;
    invert(B, Binv);

    for(j=0; j<*length; j++){
      x[j] =  exponential(doses[j], nod);
    } 
    /* eval target density */
    intLik = evalDens(x, length, n, *ntotal, grpmean, *yDinvy, *s2, prior[0], prior[1], Binv, m,
              *numScal, mean)/pow(dtm(B), 0.5);
    priordens = dbeta(nodes[i], prior[8], prior[9], 0);
    *normConst += intLik*priordens;
    if(*meanInd){
      postM[0] += mean[0]*intLik*priordens;
      postM[1] += mean[1]*intLik*priordens;   
      postM[2] += nod*intLik*priordens;
    } else {
      if(intLik*priordens > maxPost){
        maxPost = intLik*priordens;
        postM[0] = mean[0];
        postM[1] = mean[1];   
        postM[2] = nod;        
      }
    }
  }
  if(*meanInd){
    postM[0] = postM[0]/ *normConst;
    postM[1] = postM[1]/ *normConst;
    postM[2] = postM[2]/ *normConst;
  }
  *normConst = *normConst/ *nodlgth;
}
