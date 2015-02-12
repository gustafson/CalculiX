/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                     */
/*     This subroutine                                                   */
/*              Copyright (C) 2013-2015 Peter A. Gustafson               */
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

#ifdef CUDACUSP

#include <cusp/hyb_matrix.h>
#include <cusp/dia_matrix.h>
// #include <cusp/gallery/poisson.h>
#include <cusp/krylov/cg.h>
// #include <cusp/krylov/cg_m.h>
// #include <cusp/krylov/bicg.h>
// #include <cusp/krylov/bicgstab.h>
#include <cusp/version.h>
#include <cusp/print.h>
#include <cusp/array1d.h>
#include <cusp/multiply.h>
#include <cusp/precond/ainv.h> 
#include <iostream>
#include <cusp/precond/aggregation/smoothed_aggregation.h>
// #include <cusp/krylov/gmres.h>
// #include <cusp/detail/format_utils.h>
#include <thrust/copy.h>
#include <thrust/transform.h>
// #include <cusp/ell_matrix.h>
#include <time.h>

#ifdef LONGLONG
#define ITG long long
#define ITGFORMAT "lld"
#else
#define ITG int
#define ITGFORMAT "d"
#endif

// which floating point type to use
typedef ITG IndexType;
typedef double ValueType;
// typedef cusp::host_memory MemorySpace;
typedef cusp::device_memory MemorySpace;

template <typename T>
struct invsqr : public thrust::unary_function<T,T>
{
  __host__ __device__
  T operator()(const T& v) 
  {
    return T (1.0)/sqrt(v);
  }
};

template <typename T>
struct absolute : public thrust::unary_function<T,T>
{
    __host__ __device__
    T operator()(T x)
  {
    return x < 0 ? -x : x;
  }
};

extern "C"
int cudacusp(double *ad, double *au, double *adb, double *aub, double *sigma, 
	     double *b, ITG *icol, ITG *irow, ITG *neq, ITG *nzs, 
	     int *symmetryflag, int *inputformat, ITG *jq, ITG *nzs3)
{
  int cuda_major =  CUDA_VERSION / 1000;
  int cuda_minor = (CUDA_VERSION % 1000) / 10;

  int thrust_major = THRUST_MAJOR_VERSION;
  int thrust_minor = THRUST_MINOR_VERSION;

  int cusp_major = CUSP_MAJOR_VERSION;
  int cusp_minor = CUSP_MINOR_VERSION;

  clock_t timeb;
  clock_t timee;

  std::cout << " Using CUDA based on CUSP CG SOLVER\n";
  std::cout << "   CUDA   v" << cuda_major   << "." << cuda_minor   << "\n";
  std::cout << "   Thrust v" << thrust_major << "." << thrust_minor << "\n";
  std::cout << "   Cusp   v" << cusp_major   << "." << cusp_minor   << "\n";


  timeb = clock();
  // Test for non zero values
  int nvals=0;
  for (int i=0; i<*neq; i++){if (ad[i]<0) nvals++;}
  if (nvals) {thrust::transform(ad, ad+*neq, ad, absolute<ValueType>());}

  /* Fill the matrix.  
     The off diagonal triangle is columnar from ccx
     irow() identifies the row within the column
     icol() identifies the number of non zeros within the column
     Move the the next column after achieving icol() within a column. */
   
  cusp::coo_matrix<int, ValueType, cusp::host_memory> A(*neq,*neq,2*(*nzs)+*neq);
  // ASSEMBLE FULL MATRIX.  No symmetric matrix defined in CUSP //
  // Scope for off-diagonal matrix assembly
  int k=*neq; 
  int l=0; 
  // This is somewhat expensive... can it be parallelized.  Attempted below.
  for (int i = 0; i < *neq; i++){
    // i acts as a column index
    A.row_indices[i] = i; 
    A.column_indices[i] = i; 
    A.values[i] = ad[i];
    for (int j = 0; j < icol[i]; j++){
      // Looping cols
      int nrow = irow[l]-1;
      A.row_indices[k] = nrow; 
      A.column_indices[k] = i; 
      A.values[k++] = au[l];
      // Symmetry
      A.row_indices[k] = i; 
      A.column_indices[k] = nrow; 
      A.values[k++] = au[l++];
    }
  }
  

// WORKING OMP BUT NOT FASTER //   // Perform a cumsum on the column index to make a conventional csr index
// WORKING OMP BUT NOT FASTER //   thrust::exclusive_scan(icol, icol+*neq+1, icol);
// WORKING OMP BUT NOT FASTER //   {// Scope
// WORKING OMP BUT NOT FASTER //     int i,j,k,nrow;
// WORKING OMP BUT NOT FASTER // #pragma omp parallel for private(i,j,k,nrow)
// WORKING OMP BUT NOT FASTER //     for (i = 0; i < *neq; i++){
// WORKING OMP BUT NOT FASTER //       // Diagonal elements
// WORKING OMP BUT NOT FASTER //       A.row_indices[i] = i; 
// WORKING OMP BUT NOT FASTER //       A.column_indices[i] = i; 
// WORKING OMP BUT NOT FASTER //       A.values[i] = ad[i];
// WORKING OMP BUT NOT FASTER //       k=*neq+icol[i]*2;
// WORKING OMP BUT NOT FASTER //       for (j = icol[i]; j < icol[i+1]; j++){
// WORKING OMP BUT NOT FASTER // 	nrow = irow[j]-1;
// WORKING OMP BUT NOT FASTER // 	A.row_indices[k] = nrow; 
// WORKING OMP BUT NOT FASTER // 	A.column_indices[k] = i; 
// WORKING OMP BUT NOT FASTER // 	A.values[k++] = au[j];
// WORKING OMP BUT NOT FASTER // 	// Symmetry
// WORKING OMP BUT NOT FASTER // 	A.row_indices[k] = i; 
// WORKING OMP BUT NOT FASTER // 	A.column_indices[k] = nrow; 
// WORKING OMP BUT NOT FASTER // 	A.values[k++] = au[j];
// WORKING OMP BUT NOT FASTER //       }
// WORKING OMP BUT NOT FASTER //     }
// WORKING OMP BUT NOT FASTER //   }

  A.sort_by_row_and_column();
  // cusp::print(A);
  cusp::hyb_matrix<int, ValueType, MemorySpace> AA;
  try {AA = A;}
  catch(std::bad_alloc &e)
    {
      std::cerr << "bad_alloc during transfer of A to GPU" << std::endl;
      exit(-1);
    }

  
  timee = clock();
  std::cout << "  Assembled stiffness matrix on CUDA device in = " << 
    (double(timee)-double(timeb))/double(CLOCKS_PER_SEC) << " seconds\n\n";

  timeb = clock();
  // printf ("Smoothed aggregation algebraic multigrid preconditioner\n");
  // cusp::precond::aggregation::smoothed_aggregation<IndexType, ValueType, MemorySpace> MM(AA);
  printf ("Diagnonal preconditioner\n");
  cusp::precond::diagonal<ValueType, MemorySpace> MM(AA);
  // int nunsc=15;
  // printf ("Scaled bridson with %i non-zeros per row\n", nunsc);
  // cusp::precond::scaled_bridson_ainv<ValueType, MemorySpace> MM(AA, 0, nunsc);
  // printf ("Unscaled bridson with %i non-zeros per row\n", nunsc);
  // cusp::precond::bridson_ainv<ValueType, MemorySpace> MM(AA, 0, nunsc);

  timee = clock();
  std::cout << "  Preconditioning time = " << 
    (double(timee)-double(timeb))/double(CLOCKS_PER_SEC) << " seconds\n\n";
  
  // allocate storage for and copy right hand side (BB). 
  cusp::array1d<ValueType, MemorySpace> BB(*neq, 0.0);
  thrust::copy (b, b+*neq, BB.begin());

  timeb = clock();
  
  int i=50000;
  // if ((*b)<0.0){
  if (nvals){
    // Non-positive definite.  Give up quickly after spawning an answer
    // thrust::copy (ad, ad+*neq, DD.begin());
    // thrust::transform(DD.begin(), DD.end(), DD.begin(), absolute<ValueType>());
    i=0;
    printf ("There are %i negative values on the diagonal.  The attempt is abandoned.\n", nvals);
  }

  // set stopping criteria 
  // http://docs.cusp-library.googlecode.com/hg/classcusp_1_1default__monitor.html
  // ||b - A x|| <= absolute_tolerance + relative_tolerance * ||b||
  // cusp::default_monitor<ValueType> monitor(BB, i, 5e-3);
  // Abaqus uses a relative tolerance of 1e-3
  cusp::default_monitor<ValueType> monitor(BB, i, 1e-6);

  try 
    {
      // solve the linear system AA * XX = BB 
      cusp::krylov::cg(AA, BB, BB, monitor, MM); //Conjugate Gradient method
      timee = clock();
    }
  catch(std::bad_alloc &e)
    {
      std::cerr << "Couldn't solve system due to memory limits" << std::endl;
      exit(-1);
    }

  std::cout << "  CUDA iterative solver time = " << 
    (double(timee)-double(timeb))/double(CLOCKS_PER_SEC) << " seconds\n\n";

  // Copy the result to the b array
  thrust::copy (BB.begin(), BB.end(), b);

  if (monitor.converged()){
    std::cout << "Solver converged to " << monitor.relative_tolerance() << " relative tolerance";
    std::cout << " after " << monitor.iteration_count() << " iterations" << std::endl;
  }else{
    std::cout << "Solver reached iteration limit " << monitor.iteration_limit() << " before converging";
    std::cout << " to " << monitor.relative_tolerance() << " relative tolerance " << std::endl;
  }
  return 0;
}
#endif
