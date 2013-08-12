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

#ifdef CUDACUSP

#include <cusp/hyb_matrix.h>
#include <cusp/dia_matrix.h>
// #include <cusp/gallery/poisson.h>
#include <cusp/krylov/cg.h>
// #include <cusp/krylov/cg_m.h>
// #include <cusp/krylov/bicg.h>
// #include <cusp/krylov/bicgstab.h>
#include <cusp/version.h>
// #include <cusp/print.h>
#include <cusp/array1d.h>
#include <cusp/multiply.h>
#include <cusp/precond/ainv.h> 
#include <iostream>
#include <cusp/precond/smoothed_aggregation.h>
// #include <cusp/krylov/gmres.h>
// #include <cusp/detail/format_utils.h>
#include <thrust/copy.h>
#include <thrust/transform.h>
// #include <cusp/ell_matrix.h>


// which floating point type to use
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
	     double *b, int *icol, int *irow, int *neq, int *nzs, 
	     int *symmetryflag, int *inputformat, int *jq, int *nzs3)
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
  { // Scope for matrix assembly
    int k=0; // data index
    int l=0; // row index
    int m=0; // column tracker index
  
    // This is somewhat expensive... can it be parallelized.
    for (int i = 0; i < *neq; i++){
      // This is for the diagonal
      A.row_indices[k] = i; 
      A.column_indices[k] = i; 
      A.values[k++] = ad[i];

      for (int j = 0; j < icol[i]; j++){
      // This is for the off-diagonals
	int n = irow[m]-1;
	A.row_indices[k] = l; 
	A.column_indices[k] = n; 
	A.values[k++] = au[m];
	A.row_indices[k] = n; 
	A.column_indices[k] = l; 
	A.values[k++] = au[m++];
      }
      l++;
    }
  }

  A.sort_by_row_and_column();
  // cusp::print(A);
  cusp::hyb_matrix<int, ValueType, MemorySpace> AA = A;
  timee = clock();
  std::cout << "  Assembled stiffness matrix on CUDA device in = " << 
    (double(timee)-double(timeb))/double(CLOCKS_PER_SEC) << "\n\n";

  timee = clock();
  
  timeb = clock();
  printf ("Diagnonal preconditioner\n");
  cusp::precond::diagonal<ValueType, MemorySpace> MM(AA);
  timee = clock();
  std::cout << "  Preconditioning time = " << 
    (double(timee)-double(timeb))/double(CLOCKS_PER_SEC) << "\n\n";
  
  // allocate storage for and copy right hand side (BB). 
  cusp::array1d<ValueType, MemorySpace> BB(*neq, 0.0);
  thrust::copy (b, b+*neq, BB.begin());
  
  timeb = clock();
  // set stopping criteria 
  // http://docs.cusp-library.googlecode.com/hg/classcusp_1_1default__monitor.html
  // ||b - A x|| <= absolute_tolerance + relative_tolerance * ||b||
  
  int i=50000;
  // if ((*b)<0.0){
  if (nvals){
    // Non-positive definite.  Give up quickly after spawning an answer
    // thrust::copy (ad, ad+*neq, DD.begin());
    // thrust::transform(DD.begin(), DD.end(), DD.begin(), absolute<ValueType>());
    i=0;
    printf ("There are %i negative values on the diagonal.  The attempt is abandoned.\n", nvals);
  }
  cusp::verbose_monitor<ValueType> monitor(BB, i, 1e-6);


  // solve the linear system AA * XX = BB 
  cusp::krylov::cg(AA, BB, BB, monitor, MM); //Conjugate Gradient method
  timee = clock();

  std::cout << "  CUDA iterative solver time = " << 
    (double(timee)-double(timeb))/double(CLOCKS_PER_SEC) << "\n\n";

  // Copy the result to the b array
  thrust::copy (BB.begin(), BB.end(), b);

  if (!monitor.converged()){
    printf (" WARNING: Cuda Cusp did not find a solution.\n");
  }
  return 0;
}
#endif

