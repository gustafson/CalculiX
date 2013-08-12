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


template <typename Monitor>
void report_status(Monitor& monitor)
{
  if (monitor.converged())
    {
      std::cout << "  Solver converged to " << monitor.tolerance() << " tolerance";
      std::cout << " after " << monitor.iteration_count() << " iterations";
      std::cout << " (" << monitor.residual_norm() << " final residual)" << "\n";
    }
  else
    {
      std::cout << "  Solver reached iteration limit " << monitor.iteration_limit() << " before converging";
      std::cout << " to " << monitor.tolerance() << " tolerance ";
      std::cout << " (" << monitor.residual_norm() << " final residual)" << "\n";
    }
  std::cout <<  "\n\n";
}


// which floating point type to use
typedef double ValueType;
// typedef cusp::host_memory MemorySpace;
typedef cusp::device_memory MemorySpace;
// int global_recalc_cuda_M = 1;
// Can create pointers to precond matrices... can't transfer pointers to device and back as of 7/17/2013
// cusp::precond::bridson_ainv<ValueType, MemorySpace> *MM;
// cusp::precond::bridson_ainv<ValueType, cusp::host_memory> *M;




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
  /* Fill the matrix.  
     The off diagonal triangle is columnar from ccx
     irow() identifies the row within the column
     icol() identifies the number of non zeros within the column
     Move the the next column after achieving icol() within a column. /*

     Since cusp need be row sorted, we enter the transpose */
  // Create the row and column indices

  int i,j,k,l,m,n;
  int nvals=0;

  // Test for non zero values
  for (i=0; i<*neq; i++){if (ad[i]<0) nvals++;}
  if (nvals) {thrust::transform(ad, ad+*neq, ad, absolute<ValueType>());}
  
  k=0; // data index
  l=0; // row index
  m=0; // column tracker index

  // ASSEMBLE UPPER ONLY // cusp::coo_matrix<int, ValueType, cusp::host_memory> AU(*neq,*neq,*nzs);
  // ASSEMBLE UPPER ONLY // for (i = 0; i < *neq; i++){
  // ASSEMBLE UPPER ONLY //   for (j = 0; j < icol[i]; j++){
  // ASSEMBLE UPPER ONLY //     n = irow[m]-1;
  // ASSEMBLE UPPER ONLY //     AU.row_indices[k] = l; 
  // ASSEMBLE UPPER ONLY //     AU.column_indices[k] = n; 
  // ASSEMBLE UPPER ONLY //     AU.values[k++] = au[m++];
  // ASSEMBLE UPPER ONLY //   }
  // ASSEMBLE UPPER ONLY //   l++;
  // ASSEMBLE UPPER ONLY // }
  // ASSEMBLE UPPER ONLY // timee = clock();
  // ASSEMBLE UPPER ONLY // std::cout << "  Assemble upper triangular time = " << 
  // ASSEMBLE UPPER ONLY //   (double(timee)-double(timeb))/double(CLOCKS_PER_SEC) << "\n\n";

  cusp::coo_matrix<int, ValueType, cusp::host_memory> A(*neq,*neq,2*(*nzs)+*neq);
  // ASSEMBLE FULL MATRIX //
  for (i = 0; i < *neq; i++){
    A.row_indices[k] = i; 
    A.column_indices[k] = i; 
    A.values[k++] = ad[i];
    for (j = 0; j < icol[i]; j++){
      n = irow[m]-1;
      A.row_indices[k] = l; 
      A.column_indices[k] = n; 
      A.values[k++] = au[m];
      A.row_indices[k] = n; 
      A.column_indices[k] = l; 
      A.values[k++] = au[m++];
    }
    l++;
  }

  // cusp::print(A);
  A.sort_by_row_and_column();
  // cusp::print(A);
  cusp::hyb_matrix<int, ValueType, MemorySpace> AA = A;
  timee = clock();
  std::cout << "  Assembled stiffness matrix on CUDA device in = " << 
    (double(timee)-double(timeb))/double(CLOCKS_PER_SEC) << "\n\n";

  timee = clock();
  // CONVERT UPPER TO FULL MATRIX
  //
  // Start on device version
  // ON DEVICE // cusp::hyb_matrix<int, ValueType, MemorySpace> AA = AU;
  // ON DEVICE // // Bring the matrices together limiting scope as much as possible
  // ON DEVICE // {
  // ON DEVICE //   cusp::hyb_matrix<int, ValueType, MemorySpace> AAT;
  // ON DEVICE //   cusp::transpose(AA,AAT);
  // ON DEVICE //   cusp::add(AA,AAT,AA);
  // ON DEVICE // } // free AAT
  // ON DEVICE // {
  // ON DEVICE //   cusp::coo_matrix<int, ValueType, MemorySpace> DD(*neq,*neq,*neq);
  // ON DEVICE //   // Potentially not the most efficient possible
  // ON DEVICE //   thrust::sequence (DD.row_indices.begin(),DD.row_indices.end());
  // ON DEVICE //   thrust::sequence (DD.column_indices.begin(),DD.column_indices.end());
  // ON DEVICE //   thrust::copy (ad, ad+*neq, DD.values.begin());
  // ON DEVICE //   cusp::add(AA,DD,AA);
  // ON DEVICE // }
  // End on device version
  // Start on host version
  // ON HOST // {
  // ON HOST //   cusp::hyb_matrix<int, ValueType, cusp::host_memory> AAT;
  // ON HOST //   cusp::transpose(AU,AAT);
  // ON HOST //   cusp::add(AU,AAT,AU);
  // ON HOST // } // free AAT
  // ON HOST // {
  // ON HOST //   cusp::coo_matrix<int, ValueType, cusp::host_memory> DD(*neq,*neq,*neq);
  // ON HOST //   // Potentially not the most efficient possible
  // ON HOST //   thrust::sequence (DD.row_indices.begin(),DD.row_indices.end());
  // ON HOST //   thrust::sequence (DD.column_indices.begin(),DD.column_indices.end());
  // ON HOST //   thrust::copy (ad, ad+*neq, DD.values.begin());
  // ON HOST //   cusp::add(AU,DD,AU);
  // ON HOST // }
  // ON HOST // cusp::hyb_matrix<int, ValueType, MemorySpace> AA = AU;
  // End on host version

  // timee = clock();
  // std::cout << "  Time to assemble AA = " << 
  //   (double(timee)-double(timeb))/double(CLOCKS_PER_SEC) << "\n\n";
  
  // cusp::hyb_matrix<int, ValueType, cusp::host_memory> M;
  // cusp::hyb_matrix<int, ValueType, MemorySpace> MMM;
  // if (inccholpre) {
  //   int ier;
  //   AT = PreConditionCudaCusp (A,neq,nzs,&ier);  // Lower triangular, reuse AT
  //   // printf ("A");
  //   // cusp::print(A);
  //   cusp::transpose(AT,AU);
  //   cusp::multiply(AT,AU,M);
  //   
  //   // printf ("M");
  //   // cusp::print(M);
  //   MMM=M;
  // } 
  
  timeb = clock();
  // set preconditioners
  // cusp::identity_operator<ValueType, MemorySpace> MM(A.num_rows, A.num_rows);
  // AINV preconditioner, using standard drop tolerance strategy 
  // cusp::precond::scaled_bridson_ainv<ValueType, MemorySpace> MM(AA, .1);
  // printf ("Scaled bridson with .1 drop tolerarance\n");
  // cusp::precond::scaled_bridson_ainv<ValueType, MemorySpace> MM(AA, .1);
  // int nunsc = 15;
  // printf ("Scaled bridson with %i non-zeros per row\n", nunsc);
  // cusp::precond::scaled_bridson_ainv<ValueType, MemorySpace> MM(AA, 0, nunsc);
  // printf ("Unscaled bridson with %i non-zeros per row\n", nunsc);
  // cusp::precond::bridson_ainv<ValueType, MemorySpace> MM(AA, 0, nunsc);

  // The compiler warns about race conditions with smoothed.
  // printf ("Smoothed aggregation on device\n");
  // cusp::precond::smoothed_aggregation<int, ValueType, MemorySpace> MM(AA);
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
  
  i=50000;
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
  // cusp::krylov::cg(AA, BB, BB, monitor); //Conjugate Gradient method
  // cusp::krylov::bicgstab(AA, BB, BB, monitor, MM); //BiConjugate Gradient Stabilized method
  // cusp::krylov::bicg(AA, AA, BB, BB, monitor, MM, MM); //BiConjugate Gradient
  
  // solve the linear system AA * XX = BB with the GMRES.
  // Cost grows as O(n^2) with n being number of iterations.  Thus,
  // can restart with a fraction of the answer as an initial guess
  // thrust::fill( BB.begin(), BB.end(), ValueType(0) );
  // int restart = 20;
  // timeb = clock();
  // cusp::krylov::gmres(AA, BB, BB, restart, monitor, MM);
  // cusp::krylov::gmres(AA, BB, BB, restart, monitor);
  // timee = clock();
  
  timee = clock();

  std::cout << "  CUDA iterative solver time = " << 
    (double(timee)-double(timeb))/double(CLOCKS_PER_SEC) << "\n\n";

  // report status
  // report_status(monitor);
  
  // Works only with smoothed_aggregation
  // std::cout << "\nPreconditioner statistics" << "\n";
  // M.print();
  
  // Copy the result to the b array
  thrust::copy (BB.begin(), BB.end(), b);

  if (!monitor.converged()){
    printf (" WARNING: Cuda Cusp did not find a solution.\n");
  }
  return 0;
}
#endif

