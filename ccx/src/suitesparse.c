/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                     */
/*     This subroutine                                                   */
/*              Copyright (C) 2013 Peter A. Gustafson                    */
/*                                                                       */
/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                                                                       */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#ifdef SUITESPARSE

#include "cholmod.h"
#include "suitesparse.h"
#include "SuiteSparseQR_C.h"

// #include "SuiteSparseQR.hpp" // Could be used for c++ version
// #include <iostream>
// extern "C"

int suitesparsecholmod (double *ad, double *au, double *adb, double *aub, double *sigma, 
			double *b, int *icol, int *irow, int *neq, int *nzs)
{
  
  cholmod_common Common, *ccc ;
  cholmod_sparse *A ;
  cholmod_dense *X, *B;
  double one[2] = {1,0}, mone [2] = {-1,0};

  int symmetric = 1; typedef int itype;
  // int symmetric = 0; typedef long itype;

  // start CHOLMOD
  ccc = &Common;
  cholmod_start (ccc);
  
  // load A 
  // Square symmetric upper triangular. 
  cholmod_triplet *AUT;
  cholmod_triplet *ADT;
  AUT = cholmod_allocate_triplet(*neq, *neq, *nzs, 1, CHOLMOD_REAL, ccc);
  ADT = cholmod_allocate_triplet(*neq, *neq, *neq, 1, CHOLMOD_REAL, ccc);
  
  itype i,j,l,m;
  l=0; // row index
  m=0; // column tracker index
  for (i = 0; i < *neq; i++){
    for (j = 0; j < icol[i]; j++){
      ((itype*)AUT->i)[AUT->nnz] = l; 
      ((itype*)AUT->j)[AUT->nnz] = (irow[m]-1); 
      ((double*)AUT->x)[AUT->nnz] = au[m++];
      (AUT->nnz)++;
    } l++;
  }
  
  // Now add the Diagonal matrix
  for (i = 0; i < *neq; i++){
    ((itype*)ADT->i)[ADT->nnz] = i; 
    ((itype*)ADT->j)[ADT->nnz] = i;
    ((double*)ADT->x)[ADT->nnz] = ad[i];
    (ADT->nnz)++;
  }
  
  cholmod_sparse *AU;
  A = cholmod_triplet_to_sparse (ADT, 0, ccc) ;
  AU = cholmod_triplet_to_sparse (AUT, 0, ccc) ;
  cholmod_free_triplet (&AUT, ccc) ;
  cholmod_free_triplet (&ADT, ccc) ;
  
  A = cholmod_add (AU, A, one, one, 1, 1, ccc);
  cholmod_free_sparse (&AU, ccc) ;

  // B = ones (size (A,1),1)
  B = cholmod_zeros (A->nrow, 1, A->xtype, ccc) ;

  // Copy the rhs to the BB array
  for (i = 0; i < *neq; i++){
    ((double*)B->x)[i] = b[i];
  }
  
  int num_cpus=1;
#if USE_MT
  char *env;
  env=getenv("OMP_NUM_THREADS");
  if(env) {
    num_cpus=atoi(env);}
  else{
    num_cpus=1;
  }
#endif

  // Blas itself should obey OMP_NUM_THREADS.  Set SPQR threads
  // Set the number of threads
  // ccc->SPQR_nthreads = num_cpus;

  // Force the use of supernodal if you want to force CUDA
  // ccc->supernodal = CHOLMOD_SUPERNODAL;
  // Force the use of SIMPLICIAL if the matrix is indefinate with well conditioned leading minors
  // ccc->supernodal = CHOLMOD_SIMPLICIAL;

  // The system of equations has already been ordered by ccx.  Thus,
  // do not reorder the equations yet again.
  // ccc->nmethods = 1;
  // ccc->method[0].ordering = CHOLMOD_NATURAL;
  // ccc->method[0].ordering = CHOLMOD_AMD;
  // ccc->postorder = 0;

  // Set to 2, 4, or 8 depending on available memory.  (8 requires more)
  // ccc->maxrank=8;

  cholmod_factor *L ;
  L = cholmod_analyze (A, ccc) ;
  cholmod_factorize (A, L, ccc) ;

  if (ccc->status==CHOLMOD_NOT_POSDEF) {
    printf("  Stiffnessmatrix is not positive definite in column %i.\n\n", L->minor);
  }
  printf("  Factoring with cholmod_solve using %i threads\n\n", num_cpus);
  X = cholmod_solve (CHOLMOD_A, L, B, ccc) ;
  cholmod_print_common("Common", ccc);
  cholmod_free_factor (&L, ccc) ;

  printf ("\n");
  if (ccc->supernodal == CHOLMOD_SUPERNODAL){
    printf("    The use of supernodal (thus CUDA if available) is hard coded.  See suitesparse.c\n");
  }else if (ccc->supernodal == CHOLMOD_SIMPLICIAL){
    printf("    The use of simplicial is hard coded (thus CUDA is not used).  See suitesparse.c\n");
  }else{
    float sw=(ccc->fl)/(ccc->lnz);
    printf("  The supernodal switch (thus CUDA if available) is:\n    (flops/lnz) = (%0.5g/%0.5g) = %0.5g",
	   ccc->fl, ccc->lnz, sw);
    if (sw>(ccc->supernodal_switch)){
      printf(" >= %0.5g.\n    Thus supernodal is used.\n", ccc->supernodal_switch);
    }else{
      printf(" < %0.5g.\n    Thus supernodal is not used.\n", ccc->supernodal_switch);
    }
  }

  printf ("\n\n");
  // Copy the rhs to the BB array
  for (i = 0; i < *neq; i++){
    b[i] = ((double*)X->x)[i];
  }
  // free and finish CHOLMOD
  cholmod_free_sparse (&A, ccc) ;
  cholmod_free_dense (&X, ccc) ;
  cholmod_free_dense (&B, ccc) ;
  cholmod_finish (ccc) ;

  /* Note, the memory used during solve can be modified with Common->maxrank */  
  return (0) ;
}


int suitesparseqr (double *ad, double *au, double *adb, double *aub, double *sigma, 
		   double *b, int *icol, int *irow, int *neq, int *nzs)
{
  
  cholmod_common Common, *ccc ;
  cholmod_sparse *A ;
  cholmod_dense *X, *B;
  double one[2] = {1,0}, mone [2] = {-1,0};

  // int symmetric = 1; typedef int itype;
  int symmetric = 0; typedef long itype;

  // start CHOLMOD
  ccc = &Common;
  cholmod_l_start (ccc);

  // load A 
  // Square symmetric upper triangular. 
  cholmod_triplet *AUT;
  cholmod_triplet *ADT;
  AUT = cholmod_l_allocate_triplet(*neq, *neq, *nzs, 0, CHOLMOD_REAL, ccc);
  ADT = cholmod_l_allocate_triplet(*neq, *neq, *neq, 0, CHOLMOD_REAL, ccc);
  
  itype i,j,l,m;
  l=0; // row index
  m=0; // column tracker index
  for (i = 0; i < *neq; i++){
    for (j = 0; j < icol[i]; j++){
      ((itype*)AUT->i)[AUT->nnz] = l; 
      ((itype*)AUT->j)[AUT->nnz] = (irow[m]-1); 
      ((double*)AUT->x)[AUT->nnz] = au[m++];
      (AUT->nnz)++;
    } l++;
  }
  
  // Now add the Diagonal matrix
  for (i = 0; i < *neq; i++){
    ((itype*)ADT->i)[ADT->nnz] = i; 
    ((itype*)ADT->j)[ADT->nnz] = i;
    ((double*)ADT->x)[ADT->nnz] = ad[i];
    (ADT->nnz)++;
  }
  
  cholmod_sparse *AU;
  A = cholmod_l_triplet_to_sparse (ADT, 0, ccc) ;
  AU = cholmod_l_triplet_to_sparse (AUT, 0, ccc) ;
  cholmod_l_free_triplet (&AUT, ccc) ;
  cholmod_l_free_triplet (&ADT, ccc) ;
  
  A = cholmod_l_add (AU, A, one, one, 1, 1, ccc);
  A = cholmod_l_add (A, cholmod_l_transpose (AU, 2, ccc), one, one, 1, 1, ccc);
  cholmod_l_free_sparse (&AU, ccc) ;
  
  // B = ones (size (A,1),1)
  B = cholmod_l_zeros (A->nrow, 1, A->xtype, ccc) ;

  // Copy the rhs to the BB array
  for (i = 0; i < *neq; i++){
    ((double*)B->x)[i] = b[i];
  }
  
  int num_cpus=1;
#if USE_MT
  char *env;
  env=getenv("OMP_NUM_THREADS");
  if(env) {
    num_cpus=atoi(env);}
  else{
    num_cpus=1;
  }
#endif

  // Blas itself should obey OMP_NUM_THREADS.  Set SPQR threads
  // Set the number of threads
  ccc->SPQR_nthreads = num_cpus;

  // Force the use of supernodal if you want to force CUDA
  // ccc->supernodal = CHOLMOD_SUPERNODAL;

  // Set to 2, 4, or 8 depending on available memory.  (8 requires more)
  // ccc->maxrank=8;

  printf("  Factoring with SuiteSparseQR using %i threads\n", num_cpus);
  // X = SuiteSparseQR <double> (A, B, ccc); // Could be used for c++ version
  X = SuiteSparseQR_C_backslash_default (A, B, ccc);
  printf ("    Number of equations %i: rank of the matrix %ld\n\n", *neq, ccc->SPQR_istat[4]) ;
  cholmod_l_print_common("Common", ccc);

  printf ("\n");
  if (ccc->supernodal == CHOLMOD_SUPERNODAL){
    printf("    The use of supernodal (thus CUDA if available) is hard coded.  See suitesparse.c\n");
  }else{
    float sw=(ccc->fl)/(ccc->lnz);
    printf("  The supernodal switch (thus CUDA if available) is:\n    (flops/lnz) = (%0.5g/%0.5g) = %0.5g",
	   ccc->fl, ccc->lnz, sw);
    if (sw>(ccc->supernodal_switch)){
      printf(" >= %0.5g.\n    Thus supernodal was used.\n", ccc->supernodal_switch);
    }else{
      printf(" < %0.5g.\n    Thus supernodal was not used.\n", ccc->supernodal_switch);
    }
  }

  printf ("\n\n");
  // Copy the rhs to the BB array
  for (i = 0; i < *neq; i++){
    b[i] = ((double*)X->x)[i];
  }
  
  // free and finish CHOLMOD
  cholmod_l_free_sparse (&A, ccc) ;
  cholmod_l_free_dense (&X, ccc) ;
  cholmod_l_free_dense (&B, ccc) ;
  cholmod_finish (ccc) ;
  
  /* Note, the memory used during solve can be modified with Common->maxrank */
  
  return (0) ;
}

#endif