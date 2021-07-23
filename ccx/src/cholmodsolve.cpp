#include <cholmod.h>
#include "SuiteSparseQR.hpp"
#include <iostream>


// NOTE TO SELF: This was a c++ version of the cholmod solver
// interface but didn't appear to be better or faster than the c
// version.  Thus, for now, it is left as unmaintained code.

extern "C"
int cholmodsolve (double *ad, double *au, double *adb, double *aub, double *sigma, 
		  double *b, ITG *icol, ITG *irow, ITG *neq, ITG *nzs)
{
  cholmod_common com ;
  cholmod_sparse *A ;
  cholmod_dense *X, *B, *XX, *Residual ;
  double bb[*neq];
  
  // start CHOLMOD
  cholmod_start (&com) ;
  
  // load A
  // Square symmetric upper triangular.  Could be done better with lower triangular from ccx but leaving for now
  cholmod_triplet *AU;
  AU = cholmod_allocate_triplet(*neq, *neq, *nzs+*neq, 1, CHOLMOD_REAL, &com);
  // cholmod_print_triplet (AU, "AU", &com);
  
  ITG i,j,l,m;
  ITG *AUi=(ITG*)AU->i;
  ITG *AUj=(ITG*)AU->j;
  double *AUx=(double*)AU->x;

  l=0; // row index
  m=0; // column tracker index
  for (i = 0; i < *neq; i++){
    for (j = 0; j < icol[i]; j++){
      AUi[AU->nnz] = l; 
      AUj[AU->nnz] = irow[m]-1; 
      AUx[AU->nnz] = au[m++];
      (AU->nnz)++;
    }
    l++;
  }

  // Now add the Diagonal matrix
  for (i = 0; i < *neq; i++){
    AUi[AU->nnz] = i; 
    AUj[AU->nnz] = i;
    AUx[AU->nnz] = ad[i];
    (AU->nnz)++;
  }

  A = cholmod_triplet_to_sparse (AU, 0, &com) ;
  // cholmod_write_sparse(stdout, A, 0, 0,&com);

  // B = ones (size (A,1),1)
  B = cholmod_ones (A->nrow, 1, A->xtype, &com) ;
  // Copy the rhs to the BB array
  for (i = 0; i < *neq; i++){
    bb[i] = b[i];
    ((double*)B->x)[i] = bb[i];
  }
  // cholmod_write_dense(stdout, B, 0, &com);


  cholmod_factor *L ;
  L = cholmod_analyze (A, &com) ;
  cholmod_factorize (A, L, &com) ;
  X = cholmod_solve (CHOLMOD_A, L, B, &com) ;
  
  // X = A\B
  // X = SuiteSparseQR <double> (A, B, &com) ;
  // // rnorm = norm (B-A*XX)
  // Residual = cholmod_copy_dense (B, &com) ;
  // double rnorm, one [2] = {1,0}, minusone [2] = {-1,0} ;
  // cholmod_sdmult (A, 0, minusone, one, X, Residual, &com) ;
  // rnorm = cholmod_norm_dense (Residual, 2, &com) ;
  // printf ("2-norm of residual: %8.1e\n", rnorm) ;
  // printf ("rank %ld\n", &com->SPQR_istat[4]) ;
  
  // Copy the rhs to the BB array
  for (i = 0; i < *neq; i++){
    bb[i] = ((double*)X->x)[i];
  }

  // free everything and finish CHOLMOD
  // cholmod_free_dense (&Residual, &com) ;
  cholmod_free_factor (&L, &com) ;
  cholmod_free_triplet (&AU, &com) ;
  cholmod_free_sparse (&A, &com) ;
  cholmod_free_dense (&X, &com) ;
  cholmod_free_dense (&B, &com) ;
  cholmod_finish (&com) ;

  // Point the value of b at bb
  for (i = 0; i < *neq; i++){
    b[i] = bb[i];
  }
  
  return (0) ;
}

