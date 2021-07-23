/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                     */
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
// #include <cusp/ell_matrix.h>
#include <cusp/krylov/cg.h>
// #include <cusp/krylov/cg_m.h>
// #include <cusp/krylov/bicg.h>
// #include <cusp/krylov/bicgstab.h>
// #include <cusp/krylov/gmres.h>
#include <cusp/version.h>
#include <cusp/array1d.h>
#include <cusp/precond/diagonal.h> 
// #include <cusp/precond/ainv.h> 
#include <cusp/precond/aggregation/smoothed_aggregation.h>
// #include <cusp/detail/format_utils.h>
#include <thrust/copy.h>
#include <thrust/transform.h>
// #include <cusp/print.h>
#include <iostream>

#ifdef LONGLONG
#define ITG long long
#define ITGFORMAT "lld"
#else
#define ITG int
#define ITGFORMAT "d"
#endif

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
int cudacusp_thrustassembly(double *ad, double *au, double *adb, double *aub, double *sigma, 
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
  /* Fill the matrix.  ccx stores in modified compressesed sparse row
     format.  Instead of storing the pivot locations in icol, it
     stores the distance between pivots.  To make a conventional csr
     format, you must cumsum the icol vector. */

  ITG nvals=0;

  // Test for non zero values
  for (ITG i=0; i<*neq; i++){if (ad[i]<0) nvals++;}
  if (nvals) {thrust::transform(ad, ad+*neq, ad, absolute<ValueType>());}

  // Change to a zero based vector by subtracting 1
  thrust::transform(irow, irow+*nzs, thrust::make_constant_iterator(-1), irow, thrust::plus<ITG>());
  // Perform a cumsum on the column index to make a conventional csr index
  thrust::exclusive_scan(icol, icol+*neq+1, icol);

  // Create a set of "views" which act like pointers to existing memory.
  typedef typename cusp::array1d_view<ITG *> HostIndexArrayView;
  typedef typename cusp::array1d_view<ValueType *> HostValueArrayView;

  HostIndexArrayView row_offsets(icol, icol+*neq+1);
  HostIndexArrayView column_indices(irow, irow+*nzs);
  HostValueArrayView values(au, au+*nzs);
  // combine the three array1d_views into a csr_matrix_view
  typedef cusp::csr_matrix_view<HostIndexArrayView,HostIndexArrayView,HostValueArrayView> HostView;
  HostView A(*neq, *neq, *nzs, row_offsets, column_indices, values);
  
  // TRANSPOSE AND ADD ON HOST //
  cusp::coo_matrix<ITG, ValueType, cusp::host_memory> AT;
  {
    cusp::transpose(A,AT);
    cusp::add(A,AT,AT);
    
    // Create a diagonal matrix and add it to the A matrix
    // Store result in AT because A is just a matrix view to the original memory
    cusp::dia_matrix<ITG, ValueType, cusp::host_memory> D(*neq,*neq,*neq,1);
    D.diagonal_offsets[0]=0;
    for (ITG i=0; i<*neq; i++){D.values(i,0)=ad[i];}
    cusp::add(AT,D,AT);
    // Free DD
  }
  // Move to the device
  AT.sort_by_row_and_column();
  cusp::hyb_matrix<ITG, ValueType, MemorySpace> AA = AT;

  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS //
  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS // // Move to the device
  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS // cusp::hyb_matrix<ITG, ValueType, MemorySpace> AA = A;
  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS // // Bring the matrices together limiting scope as much as possible
  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS // {
  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS //   cusp::hyb_matrix<ITG, ValueType, MemorySpace> AAT;
  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS //   cusp::transpose(AA,AAT);
  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS //   cusp::add(AA,AAT,AA);
  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS // } // free AAT
  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS // {
  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS //   // Create a diagonal matrix and add it to the A matrix
  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS //   cusp::dia_matrix<ITG, ValueType, MemorySpace> DD(*neq,*neq,*neq,1);
  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS //   DD.diagonal_offsets[0]=0;
  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS //   for (ITG i=0; i<*neq; i++){DD.values(i,0)=ad[i];}
  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS //   cusp::add(AA,DD,AA);
  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS //   // Free DD
  // TRANSPOSE AND ADD ON DEVICE: EXHAUSTS DEVICE MEMORY FOR LARGE MODELS // }

  // cusp::print(AA);

  timee = clock();
  std::cout << "  Assembled stiffness matrix on CUDA device in = " << 
    (double(timee)-double(timeb))/double(CLOCKS_PER_SEC) << " seconds\n\n";

  timeb = clock();
  // set preconditioner
  printf ("Diagnonal preconditioner\n");
  cusp::precond::diagonal<ValueType, MemorySpace> MM(AA);
  timee = clock();
  std::cout << "  Preconditioning time = " << 
    (double(timee)-double(timeb))/double(CLOCKS_PER_SEC) << " seconds\n\n";
  
  // allocate storage for and copy right hand side (BB). 
  cusp::array1d<ValueType, MemorySpace> BB(*neq, 0);
  thrust::copy (b, b+*neq, BB.begin());
  
  // set stopping criteria 
  ITG i=50000;
  if (nvals){
    // Non-positive definite.  Give up quickly after spawning an answer
    i=0;
    printf ("There are %i negative values on the diagonal.  The attempt is abandoned.\n", nvals);
  }
  cusp::verbose_monitor<ValueType> monitor(BB, i, 1e-6);
    
  // solve the linear system AA * XX = BB 
  timeb = clock();
  cusp::krylov::cg(AA, BB, BB, monitor, MM); //Conjugate Gradient method
  timee = clock();

  std::cout << "  CUDA iterative solver time = " << 
    (double(timee)-double(timeb))/double(CLOCKS_PER_SEC) << " seconds\n\n";

  thrust::copy (BB.begin(), BB.end(), b);

  if (!monitor.converged()){
    printf (" WARNING: Cuda Cusp did not find a solution.\n");
  }
  return 0;
}
#endif

