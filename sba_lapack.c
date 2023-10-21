/////////////////////////////////////////////////////////////////////////////////
//// 
////  Linear algebra operations for the sba package
////  Copyright (C) 2004  Manolis Lourakis (lourakis@ics.forth.gr)
////  Institute of Computer Science, Foundation for Research & Technology - Hellas
////  Heraklion, Crete, Greece.
////
////  This program is free software; you can redistribute it and/or modify
////  it under the terms of the GNU General Public License as published by
////  the Free Software Foundation; either version 2 of the License, or
////  (at your option) any later version.
////
////  This program is distributed in the hope that it will be useful,
////  but WITHOUT ANY WARRANTY; without even the implied warranty of
////  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
////  GNU General Public License for more details.
////
///////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sba.h"

/* declarations of LAPACK routines employed */

/* QR decomposition */
extern int dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);
extern int dorgqr_(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

/* solution of triangular system */
extern int dtrtrs_(char *uplo, char *trans, char *diag, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info);

/* cholesky decomposition */
extern int dpotf2_(char *uplo, int *n, double *a, int *lda, int *info);
extern int dpotrf_(char *uplo, int *n, double *a, int *lda, int *info); /* block version of dpotf2 */

/* LU decomposition, linear system solution and matrix inversion */
extern int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
extern int dgetrs_(char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
extern int dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);

/* SVD */
extern int dgesvd_(char *jobu, char *jobvt, int *m, int *n,
           double *a, int *lda, double *s, double *u, int *ldu,
           double *vt, int *ldvt, double *work, int *lwork,
           int *info);

/* lapack 3.0 routine, faster than dgesvd_() */
extern int dgesdd_(char *jobz, int *m, int *n, double *a, int *lda,
           double *s, double *u, int *ldu, double *vt, int *ldvt,
           double *work, int *lwork, int *iwork, int *info);


/* Bunch-Kaufman factorization of a real symmetric matrix A and solution of linear systems */
extern int dsytrf_(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
extern int dsytrs_(char *uplo, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);


/*
 * This function returns the solution of Ax = b
 *
 * The function is based on QR decomposition with explicit computation of Q:
 * If A=Q R with Q orthogonal and R upper triangular, the linear system becomes
 * Q R x = b or R x = Q^T b.
 *
 * A is mxn, b is mx1. Argument iscolmaj specifies whether A is
 * stored in column or row major order. Note that if iscolmaj==1
 * this function modifies A!
 *
 * The function returns 0 in case of error, 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 */
int sba_Axb_QR(double *A, double *B, double *x, int m, int iscolmaj)
{
static double *buf=NULL;
static int buf_sz=0;

double *a, *qtb, *r, *tau, *work;
int a_sz, qtb_sz, r_sz, tau_sz, tot_sz;
register int i, j;
int info, worksz, nrhs=1;
register double sum;
   
    /* calculate required memory size */
    a_sz=(iscolmaj)? 0 : m*m;
    qtb_sz=m;
    r_sz=m*m; /* only the upper triangular part really needed */
    tau_sz=m;
    worksz=3*m; /* this is probably too much */
    tot_sz=a_sz + qtb_sz + r_sz + tau_sz + worksz;

    if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
      if(buf) free(buf); /* free previously allocated memory */

      buf_sz=tot_sz;
      buf=(double *)malloc(buf_sz*sizeof(double));
      if(!buf){
        fprintf(stderr, "memory allocation in sba_Axb_QR() failed!\n");
        exit(1);
      }
    }

    if(!iscolmaj){
    	a=buf;
	    /* store A (column major!) into a */
	    for(i=0; i<m; i++)
		    for(j=0; j<m; j++)
			    a[i+j*m]=A[i*m+j];
    }
    else a=A; /* no copying required */

    qtb=buf+a_sz;
    r=qtb+qtb_sz;
    tau=r+r_sz;
    work=tau+tau_sz;

  /* QR decomposition of A */
  dgeqrf_((int *)&m, (int *)&m, a, (int *)&m, tau, work, (int *)&worksz, (int *)&info);
  /* error treatment */
  if(info!=0){
    if(info<0){
      fprintf(stderr, "LAPACK error: illegal value for argument %d of dgeqrf in sba_Axb_QR()\n", -info);
      exit(1);
    }
    else{
      fprintf(stderr, "Unknown LAPACK error %d for dgeqrf in sba_Axb_QR()\n", info);
      return 0;
    }
  }

  /* R is now stored in the upper triangular part of a; copy it in r so that dorgqr() below won't destroy it */
  for(i=0; i<r_sz; i++)
    r[i]=a[i];

  /* compute Q using the elementary reflectors computed by the above decomposition */
  dorgqr_((int *)&m, (int *)&m, (int *)&m, a, (int *)&m, tau, work, (int *)&worksz, (int *)&info);
  if(info!=0){
    if(info<0){
      fprintf(stderr, "LAPACK error: illegal value for argument %d of dorgqr in sba_Axb_QR()\n", -info);
      exit(1);
    }
    else{
      fprintf(stderr, "Unknown LAPACK error (%d) in sba_Axb_QR()\n", info);
      return 0;
    }
  }

  /* Q is now in a; compute Q^T b in qtb */
  for(i=0; i<m; i++){
    for(j=0, sum=0.0; j<m; j++)
      sum+=a[i*m+j]*B[j];
    qtb[i]=sum;
  }

  /* solve the linear system R x = Q^t b */
  dtrtrs_("U", "N", "N", (int *)&m, (int *)&nrhs, r, (int *)&m, qtb, (int *)&m, &info);
  /* error treatment */
  if(info!=0){
    if(info<0){
      fprintf(stderr, "LAPACK error: illegal value for argument %d of dtrtrs in sba_Axb_QR()\n", -info);
      exit(1);
    }
    else{
      fprintf(stderr, "LAPACK error: the %d-th diagonal element of A is zero (singular matrix) in sba_Axb_QR()\n", info);
      return 0;
    }
  }

	/* copy the result in x */
	for(i=0; i<m; i++)
    x[i]=qtb[i];

	return 1;
}

/*
 * This function returns the solution of Ax = b
 *
 * The function is based on QR decomposition without computation of Q:
 * If A=Q R with Q orthogonal and R upper triangular, the linear system becomes
 * (A^T A) x = A^T b or (R^T Q^T Q R) x = A^T b or (R^T R) x = A^T b.
 * This amounts to solving R^T y = A^T b for y and then R x = y for x
 * Note that Q does not need to be explicitly computed
 *
 * A is mxn, b is mx1. Argument iscolmaj specifies whether A is
 * stored in column or row major order. Note that if iscolmaj==1
 * this function modifies A!
 *
 * The function returns 0 in case of error, 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 */
int sba_Axb_QRnoQ(double *A, double *B, double *x, int m, int iscolmaj)
{
static double *buf=NULL;
static int buf_sz=0;

double *a, *atb, *tau, *work;
int a_sz, atb_sz, tau_sz, tot_sz;
register int i, j;
int info, worksz, nrhs=1;
register double sum;
   
    /* calculate required memory size */
    a_sz=(iscolmaj)? 0 : m*m;
    atb_sz=m;
    tau_sz=m;
    worksz=3*m; /* this is probably too much */
    tot_sz=a_sz + atb_sz + tau_sz + worksz;

    if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
      if(buf) free(buf); /* free previously allocated memory */

      buf_sz=tot_sz;
      buf=(double *)malloc(buf_sz*sizeof(double));
      if(!buf){
        fprintf(stderr, "memory allocation in sba_Axb_QRnoQ() failed!\n");
        exit(1);
      }
    }

    if(!iscolmaj){
    	a=buf;
	/* store A (column major!) into a */
	for(i=0; i<m; i++)
		for(j=0; j<m; j++)
			a[i+j*m]=A[i*m+j];
    }
    else a=A; /* no copying required */

    atb=buf+a_sz;
    tau=atb+atb_sz;
    work=tau+tau_sz;

  /* compute A^T b in atb */
  for(i=0; i<m; i++){
    for(j=0, sum=0.0; j<m; j++)
      sum+=a[i*m+j]*B[j];
    atb[i]=sum;
  }

  /* QR decomposition of A */
  dgeqrf_((int *)&m, (int *)&m, a, (int *)&m, tau, work, (int *)&worksz, (int *)&info);
  /* error treatment */
  if(info!=0){
    if(info<0){
      fprintf(stderr, "LAPACK error: illegal value for argument %d of dgeqrf in sba_Axb_QRnoQ()\n", -info);
      exit(1);
    }
    else{
      fprintf(stderr, "Unknown LAPACK error %d for dgeqrf in sba_Axb_QRnoQ()\n", info);
      return 0;
    }
  }

  /* R is stored in the upper triangular part of a */

  /* solve the linear system R^T y = A^t b */
  dtrtrs_("U", "T", "N", (int *)&m, (int *)&nrhs, a, (int *)&m, atb, (int *)&m, &info);
  /* error treatment */
  if(info!=0){
    if(info<0){
      fprintf(stderr, "LAPACK error: illegal value for argument %d of dtrtrs in sba_Axb_QRnoQ()\n", -info);
      exit(1);
    }
    else{
      fprintf(stderr, "LAPACK error: the %d-th diagonal element of A is zero (singular matrix) in sba_Axb_QRnoQ()\n", info);
      return 0;
    }
  }

  /* solve the linear system R x = y */
  dtrtrs_("U", "N", "N", (int *)&m, (int *)&nrhs, a, (int *)&m, atb, (int *)&m, &info);
  /* error treatment */
  if(info!=0){
    if(info<0){
      fprintf(stderr, "LAPACK error: illegal value for argument %d of dtrtrs in sba_Axb_QRnoQ()\n", -info);
      exit(1);
    }
    else{
      fprintf(stderr, "LAPACK error: the %d-th diagonal element of A is zero (singular matrix) in sba_Axb_QRnoQ()\n", info);
      return 0;
    }
  }

	/* copy the result in x */
	for(i=0; i<m; i++)
    x[i]=atb[i];

	return 1;
}

/*
 * This function returns the solution of Ax=b
 *
 * The function assumes that A is symmetric & postive definite and employs
 * the Cholesky decomposition:
 * If A=U^T U with U upper triangular, the system to be solved becomes
 * (U^T U) x = b
 * This amount to solving U^T y = b for y and then U x = y for x
 *
 * A is mxn, b is mx1. Argument iscolmaj specifies whether A is
 * stored in column or row major order. Note that if iscolmaj==1
 * this function modifies A and B!
 *
 * The function returns 0 in case of error, 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 */
int sba_Axb_Chol(double *A, double *B, double *x, int m, int iscolmaj)
{
static double *buf=NULL;
static int buf_sz=0;

double *a, *b;
int a_sz, b_sz, tot_sz;
register int i, j;
int info, nrhs=1;
   
    /* calculate required memory size */
    a_sz=(iscolmaj)? 0 : m*m;
    b_sz=(iscolmaj)? 0 : m;
    tot_sz=a_sz + b_sz;

    if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
      if(buf) free(buf); /* free previously allocated memory */

      buf_sz=tot_sz;
      buf=(double *)malloc(buf_sz*sizeof(double));
      if(!buf){
        fprintf(stderr, "memory allocation in sba_Axb_Chol() failed!\n");
        exit(1);
      }
    }

    if(!iscolmaj){
    	a=buf;
    	b=a+a_sz;

  /* store A (column major!) into a anb B into b */
	for(i=0; i<m; i++){
		for(j=0; j<m; j++)
			a[i+j*m]=A[i*m+j];

    	b[i]=B[i];
	}
    }
    else{ /* no copying is necessary */
      a=A;
      b=B;
    }

  /* Cholesky decomposition of A */
  dpotf2_("U", (int *)&m, a, (int *)&m, (int *)&info);
  /* error treatment */
  if(info!=0){
    if(info<0){
      fprintf(stderr, "LAPACK error: illegal value for argument %d of dpotf2 in sba_Axb_Chol()\n", -info);
      exit(1);
    }
    else{
      fprintf(stderr, "LAPACK error: the leading minor of order %d is not positive definite,\nthe factorization could not be completed for dpotf2 in sba_Axb_Chol()\n", info);
      return 0;
    }
  }

  /* solve the linear system U^T y = b */
  dtrtrs_("U", "T", "N", (int *)&m, (int *)&nrhs, a, (int *)&m, b, (int *)&m, &info);
  /* error treatment */
  if(info!=0){
    if(info<0){
      fprintf(stderr, "LAPACK error: illegal value for argument %d of dtrtrs in sba_Axb_Chol()\n", -info);
      exit(1);
    }
    else{
      fprintf(stderr, "LAPACK error: the %d-th diagonal element of A is zero (singular matrix) in sba_Axb_Chol()\n", info);
      return 0;
    }
  }

  /* solve the linear system U x = y */
  dtrtrs_("U", "N", "N", (int *)&m, (int *)&nrhs, a, (int *)&m, b, (int *)&m, &info);
  /* error treatment */
  if(info!=0){
    if(info<0){
      fprintf(stderr, "LAPACK error: illegal value for argument %d of dtrtrs in sba_Axb_Chol()\n", -info);
      exit(1);
    }
    else{
      fprintf(stderr, "LAPACK error: the %d-th diagonal element of A is zero (singular matrix) in sba_Axb_Chol()\n", info);
      return 0;
    }
  }

	/* copy the result in x */
	for(i=0; i<m; i++)
    x[i]=b[i];

	return 1;
}

/*
 * This function returns the solution of Ax = b
 *
 * The function employs LU decomposition:
 * If A=L U with L lower and U upper triangular, then the original system
 * amounts to solving
 * L y = b, U x = y
 *
 * A is mxn, b is mx1. Argument iscolmaj specifies whether A is
 * stored in column or row major order. Note that if iscolmaj==1
 * this function modifies A and B!
 *
 * The function returns 0 in case of error,
 * 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 */
int sba_Axb_LU(double *A, double *B, double *x, int m, int iscolmaj)
{
static double *buf=NULL;
static int buf_sz=0;

int a_sz, ipiv_sz, b_sz, tot_sz;
register int i, j;
int info, *ipiv, nrhs=1;
double *a, *b;
   
    /* calculate required memory size */
    ipiv_sz=m;
    a_sz=(iscolmaj)? 0 : m*m;
    b_sz=(iscolmaj)? 0 : m;
    tot_sz=ipiv_sz*sizeof(int) + (a_sz + b_sz)*sizeof(double);

    if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
      if(buf) free(buf); /* free previously allocated memory */

      buf_sz=tot_sz;
      buf=(double *)malloc(buf_sz);
      if(!buf){
        fprintf(stderr, "memory allocation in sba_Axb_LU() failed!\n");
        exit(1);
      }
    }

    ipiv=(int *)buf;
    if(!iscolmaj){
    	a=(double *)(ipiv + ipiv_sz);
    	b=a+a_sz;

    /* store A (column major!) into a and B into b */
	  for(i=0; i<m; i++){
		  for(j=0; j<m; j++)
        	a[i+j*m]=A[i*m+j];

      	b[i]=B[i];
    	}
    }
    else{ /* no copying is necessary */
      a=A;
      b=B;
    }

  /* LU decomposition for A */
	dgetrf_((int *)&m, (int *)&m, a, (int *)&m, ipiv, (int *)&info);  
	if(info!=0){
		if(info<0){
			fprintf(stderr, "argument %d of dgetrf illegal in sba_Axb_LU()\n", -info);
			exit(1);
		}
		else{
			fprintf(stderr, "singular matrix A for dgetrf in sba_Axb_LU()\n");
			return 0;
		}
	}

  /* solve the system with the computed LU */
  dgetrs_("N", (int *)&m, (int *)&nrhs, a, (int *)&m, ipiv, b, (int *)&m, (int *)&info);
	if(info!=0){
		if(info<0){
			fprintf(stderr, "argument %d of dgetrs illegal in sba_Axb_LU()\n", -info);
			exit(1);
		}
		else{
			fprintf(stderr, "unknown error for dgetrs in sba_Axb_LU()\n");
			return 0;
		}
	}

	/* copy the result in x */
	for(i=0; i<m; i++){
		x[i]=b[i];
	}

	return 1;
}

/*
 * This function returns the solution of Ax = b
 *
 * The function is based on SVD decomposition:
 * If A=U D V^T with U, V orthogonal and D diagonal, the linear system becomes
 * (U D V^T) x = b or x=V D^{-1} U^T b
 * Note that V D^{-1} U^T is the pseudoinverse A^+
 *
 * A is mxn, b is mx1. Argument iscolmaj specifies whether A is
 * stored in column or row major order. Note that if iscolmaj==1
 * this function modifies A!
 *
 * The function returns 0 in case of error, 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 */
int sba_Axb_SVD(double *A, double *B, double *x, int m, int iscolmaj)
{
static double *buf=NULL;
static int buf_sz=0;
static double eps=-1.0;

register int i, j;
double *a, *u, *s, *vt, *work;
int a_sz, u_sz, s_sz, vt_sz, tot_sz;
double thresh, one_over_denom;
register double sum;
int info, rank, worksz, *iwork, iworksz;
   
  /* calculate required memory size */
  worksz=16*m; /* more than needed */
  iworksz=8*m;
  a_sz=(!iscolmaj)? m*m : 0;
  u_sz=m*m; s_sz=m; vt_sz=m*m;

  tot_sz=iworksz*sizeof(int) + (a_sz + u_sz + s_sz + vt_sz + worksz)*sizeof(double);

  if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
    if(buf) free(buf); /* free previously allocated memory */

    buf_sz=tot_sz;
    buf=(double *)malloc(buf_sz);
    if(!buf){
      fprintf(stderr, "memory allocation in sba_Axb_SVD() failed!\n");
      exit(1);
    }
  }

  iwork=(int *)buf;
  if(!iscolmaj){
    a=(double *)(iwork+iworksz);
    /* store A (column major!) into a */
    for(i=0; i<m; i++)
      for(j=0; j<m; j++)
        a[i+j*m]=A[i*m+j];
  }
  else{
    a=A; /* no copying required */
  }

  u=((double *)(iwork+iworksz)) + a_sz;
  s=u+u_sz;
  vt=s+s_sz;
  work=vt+vt_sz;

  /* SVD decomposition of A */
  dgesvd_("A", "A", (int *)&m, (int *)&m, a, (int *)&m, s, u, (int *)&m, vt, (int *)&m, work, (int *)&worksz, &info);
  //dgesdd_("A", (int *)&m, (int *)&m, a, (int *)&m, s, u, (int *)&m, vt, (int *)&m, work, (int *)&worksz, iwork, &info);

  /* error treatment */
  if(info!=0){
    if(info<0){
      fprintf(stderr, "LAPACK error: illegal value for argument %d of dgesdd/dgesvd in sba_Axb_SVD()\n", -info);
      exit(1);
    }
    else{
      fprintf(stderr, "LAPACK error: dgesdd (dbdsdc)/dgesvd (dbdsqr) failed to converge in sba_Axb_SVD() [info=%d]\n", info);

      return 0;
    }
  }

  if(eps<0.0){
    /* compute machine epsilon */
    for(eps=1.0; eps+1.0!=1.0; eps*=0.5)
                    ;
    eps*=2.0;
  }

  /* compute the pseudoinverse in a */
  memset(a, 0, m*m*sizeof(double)); /* initialize to zero */
  for(rank=0, thresh=eps*s[0]; rank<m && s[rank]>thresh; rank++){
    one_over_denom=1.0/s[rank];

    for(j=0; j<m; j++)
      for(i=0; i<m; i++)
        a[i*m+j]+=vt[rank+i*m]*u[j+rank*m]*one_over_denom;
  }

	/* compute A^+ b in x */
	for(i=0; i<m; i++){
	  for(j=0, sum=0.0; j<m; j++)
      sum+=a[i*m+j]*B[j];
    x[i]=sum;
  }

	return 1;
}

/*
 * This function returns the solution of Ax = b for a real symmetric matrix A
 *
 * The function is based on Bunch-Kaufman factorization:
 * A is factored as U*D*U^T where U is upper triangular and
 * D symmetric and block diagonal
 *
 * A is mxn, b is mx1. Argument iscolmaj specifies whether A is
 * stored in column or row major order. Note that if iscolmaj==1
 * this function modifies A and B!
 *
 * The function returns 0 in case of error,
 * 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 */
int sba_Axb_BK(double *A, double *B, double *x, int m, int iscolmaj)
{
static double *buf=NULL;
static int buf_sz=0;

int a_sz, ipiv_sz, b_sz, work_sz, tot_sz;
register int i, j;
int info, *ipiv, nrhs=1;
double *a, *b, *work;
   
    /* calculate required memory size */
    ipiv_sz=m;
    a_sz=(iscolmaj)? 0 : m*m;
    b_sz=(iscolmaj)? 0 : m;
    work_sz=16*m; /* this is probably too much */
    tot_sz=ipiv_sz*sizeof(int) + (a_sz + b_sz + work_sz)*sizeof(double);

    if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
      if(buf) free(buf); /* free previously allocated memory */

      buf_sz=tot_sz;
      buf=(double *)malloc(buf_sz);
      if(!buf){
        fprintf(stderr, "memory allocation in sba_Axb_BK() failed!\n");
        exit(1);
      }
    }

    ipiv=(int *)buf;
    if(!iscolmaj){
    	a=(double *)(ipiv + ipiv_sz);
    	b=a+a_sz;
    	work=b+b_sz;

    /* store A (column major!) into a and B into b */
	  for(i=0; i<m; i++){
		  for(j=0; j<m; j++)
        	a[i+j*m]=A[i*m+j];

      	b[i]=B[i];
    	}
    }
    else{ /* no copying is necessary */
      a=A;
      b=B;
    	work=(double *)(ipiv + ipiv_sz);
    }

  /* factorization for A */
	dsytrf_("U", (int *)&m, a, (int *)&m, ipiv, work, (int *)&work_sz, (int *)&info);
	if(info!=0){
		if(info<0){
			fprintf(stderr, "argument %d of dsytrf illegal in sba_Axb_BK()\n", -info);
			exit(1);
		}
		else{
			fprintf(stderr, "singular block diagonal matrix D for dsytrf in sba_Axb_BK() [D(%d, %d) is zero]\n", info, info);
			return 0;
		}
	}

  /* solve the system with the computed factorization */
  dsytrs_("U", (int *)&m, (int *)&nrhs, a, (int *)&m, ipiv, b, (int *)&m, (int *)&info);
	if(info!=0){
		if(info<0){
			fprintf(stderr, "argument %d of dsytrs illegal in sba_Axb_BK()\n", -info);
			exit(1);
		}
		else{
			fprintf(stderr, "unknown error for dsytrs in sba_Axb_BK()\n");
			return 0;
		}
	}

	/* copy the result in x */
	for(i=0; i<m; i++){
		x[i]=b[i];
	}

	return 1;
}

/*
 * This function computes the inverse of a square matrix A into B
 * using LU decomposition
 *
 * The function returns 0 in case of error (e.g. A is singular),
 * 1 if successfull
 *
 * This function is often called repetitively to solve problems of identical
 * dimensions. To avoid repetitive malloc's and free's, allocated memory is
 * retained between calls and free'd-malloc'ed when not of the appropriate size.
 */
int sba_mat_invert(double *A, double *B, int m)
{
static double *buf=NULL;
static int buf_sz=0;

int a_sz, ipiv_sz, work_sz, tot_sz;
register int i, j;
int info, *ipiv;
double *a, *work;
   
    /* calculate required memory size */
	  ipiv_sz=m;
    a_sz=m*m;
    work_sz=16*m; /* this is probably too much */
    tot_sz=ipiv_sz*sizeof(int) + (a_sz + work_sz)*sizeof(double); 

    if(tot_sz>buf_sz){ /* insufficient memory, allocate a "big" memory chunk at once */
      if(buf) free(buf); /* free previously allocated memory */

      buf_sz=tot_sz;
      buf=(double *)malloc(buf_sz);
      if(!buf){
        fprintf(stderr, "memory allocation in sba_mat_invert() failed!\n");
        exit(1);
      }
    }

	  ipiv=(int *)buf;
    a=(double *)(ipiv + ipiv_sz);
    work=a+a_sz;

  /* store A (column major!) into a */
	for(i=0; i<m; i++)
		for(j=0; j<m; j++)
			a[i+j*m]=A[i*m+j];

  /* LU decomposition for A */
	dgetrf_((int *)&m, (int *)&m, a, (int *)&m, ipiv, (int *)&info);  
	if(info!=0){
		if(info<0){
			fprintf(stderr, "argument %d of dgetrf illegal in sba_mat_invert()\n", -info);
			exit(1);
		}
		else{
			fprintf(stderr, "singular matrix A for dgetrf in sba_mat_invert()\n");
			return 0;
		}
	}

  /* (A)^{-1} from LU */
	dgetri_((int *)&m, a, (int *)&m, ipiv, work, (int *)&work_sz, (int *)&info);
	if(info!=0){
		if(info<0){
			fprintf(stderr, "argument %d of dgetri illegal in sba_mat_invert()\n", -info);
			exit(1);
		}
		else{
			fprintf(stderr, "singular matrix A for dgetri in sba_mat_invert()\n");
			return 0;
		}
	}

	/* store (A)^{-1} in B */
	for(i=0; i<m; i++)
		for(j=0; j<m; j++)
      B[i*m+j]=a[i+j*m];

	return 1;
}
