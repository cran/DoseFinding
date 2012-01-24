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
#include<Rdefines.h>
#include<Rinternals.h>
#include<R_ext/Lapack.h>
#include<R_ext/BLAS.h>
#include<R_ext/Applic.h>


/* Calculate */

void rank1vec(double *x, int *dim, double *alpha, double *A){
  // calculates alpha*x*x'+A
  char *uplo="U";
  int incx=1;
  F77_CALL(dsyr)(uplo, dim, alpha, x, &incx, A, dim);
}

// calculate design matrix
void calcSum(double *x, int *p, double *xsub, double *w, 
             int *k, double *A, int *incrx){
  // k - number of doses (length of w)
  // p - row (and column) nedimsion of A
  int i,j=0;
  for(i=0;i<*k;i++){
    for(j=0;j<*p;j++){
      xsub[j] = x[*incrx+4*i+j]; // general version *incrx+*p*i+j
    }
    rank1vec(xsub, p, &w[i], A);
  }
  // recover complete symmetric matrix from upper triang. part
  for(i=0;i<*p;i++){
    for(j=0;j<i;j++){
      A[*p*j+i] = A[*p*i+j];
    }
  }
}

void calcDetGinv(double *X, int *dim, double *s, double *U,
         double *VT, double *work, double *tol, 
         int *type, double *resD){
  // dimension work 10*dim
  int i,j,k,nonzero=*dim;
  char jobu = 'A';
  char jobvt = 'A';
  int info,lwork;

  lwork = 4*10; // in general: *dim*10

  // Calculate singular value decomposition
  F77_CALL(dgesvd)(&jobu, &jobvt, dim, dim, X, dim,
          s, U, dim, VT, dim, work, &lwork, &info);

  if((*type == 1) || (*type == 3)){ // calculate g-inverse
    for (i = 1; i < *dim; i++){
      if (s[i] < *tol*s[0]){
    nonzero = i;
    break;
      }
    }
    for (i = 0; i < *dim; i++){
      for (j = 0; j < nonzero; j++){
    U[j**dim + i] = U[j**dim+i] * 1.0/s[j];
      }
    }
    for (i = 0; i < *dim; i++){
      for (j = i; j < *dim; j++){
    X[j**dim+i] = 0.0;
    for (k=0; k < nonzero; k++){
      X[j**dim+i] += VT[i**dim+k] * U[k**dim+j];
    }
      }
    }
  }
  if((*type == 2) || (*type == 3)){ // calculate determinant
    *resD = 1.0;
    for (i = 0; i < *dim; i++){
      *resD *= s[i];
    }
  }
} 

void calcQuadform(double *beta, double *Q, int *dim, double *out, int *incrbeta){
  // calculates quadratic form beta'Qbeta
  // Q = (Q11,Q12,...,Q22,Q23,... (only upper triangular part of sym. matrix)
  int i,j;
  for(i=0;i<*dim;i++){
    for(j=i;j<*dim;j++){
      if(i==j){
        *out += Q[*dim*j+i]*beta[*incrbeta+i]*beta[*incrbeta+i];
      } else {
        *out += 2*Q[*dim*j+i]*beta[*incrbeta+i]*beta[*incrbeta+j];
      }
    }
  }
}

void getweights(double *w2, double *n2, double *nold, int *k){
  int i;
  double n1=0;
  for(i=0;i<*k;i++){
    n1 += nold[i];
  }
  for(i=0;i<*k;i++){
    w2[i] = (*n2*w2[i] + nold[i])/(*n2 + n1);
  }
}

void setzero(double *x, int dim){
  int i;
  for(i=0;i<dim;i++){
    x[i] = 0.0;
  }
}

void critfunc(double *x, int *p, double *xsub, int *k, double *probs, int *M,
              double *w, double *n2, double *nold,
              double *A, double *s, double *U, double *VT,
              double *work, double *tol, double *bvec, int *type,
	      int *stand, double *res){
  // x - contains gradient vectors (4 cells reserved for each model)
  // p - number of parameters (dim A)
  // k - number of dose-levels
  // w - design
  // xsub - double of length k
  // s,U,VT, work, tol - needed for SVD decomp.
  //    - dim(U)=dim(S): p*p, dim(work): (10p)*(10p) 
  // type - 1: MED, 2: Dopt, 3: MED&Dopt
  int m,incx,incb;
  double resM=0,resD=0,fracp=0;
  *res = 0.0;
  // calculate weight vector
  getweights(w, n2, nold, k);
  for(m=0;m<*M;m++){
    incx = *k*m*4;incb = m*4;
    setzero(A, 16);resM = 0.0;
    // calulate matrix 
    calcSum(x, &p[m], xsub, w, k, A, &incx);
    // calculate det and/or  MP-Inverse 
    calcDetGinv(A, &p[m], s, U, VT, work, tol, type, &resD);
    if(*type == 1){
      // calculate quadratic form (for MED designs)
      calcQuadform(bvec, A, &p[m], &resM, &incb);
      *res += probs[m]*log(resM);
    }
    if(*type == 2){
      if(*stand == 1){
	fracp = (double) p[m];
	*res += probs[m]*(-log(resD)/fracp);	
      } else {
	*res += probs[m]*(-log(resD));
      }
    }
    if(*type == 3){
      // calculate quadratic form (for MED designs)
      calcQuadform(bvec, A, &p[m], &resM, &incb);
      if(*stand == 1){
	fracp = (double) p[m];
	*res += probs[m]*(-0.5*log(resD)/fracp+0.5*log(resM));
      } else {
	*res += probs[m]*(-0.5*log(resD)+0.5*log(resM));
      }
    }
  }
}
