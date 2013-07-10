/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2005 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#ifdef PARDISO

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"
#include "pardiso.h"

int *irowpardisoas=NULL,*pointersas=NULL,iparmas[64];
long long ptas[64];
double *aupardisoas=NULL;

/* factorization, solving and cleaning with PARDISO for
   real unsymmetric matrices */

void pardiso_factor_as(double *ad, double *au, double *adb, double *aub, 
                double *sigma,int *icol, int *irow, 
                int *neq, int *nzs, int *jq){

  char *env;
  int i,j,k,maxfct=1,mnum=1,mtype=11,phase=12,nrhs=1,*perm=NULL,
      msglvl=0,error=0,ifortran,lfortran,index,id;
  long long ndim;
  double *b=NULL,*x=NULL;

  printf(" Factoring the system of equations using the asymmetric pardiso solver\n");

  iparmas[0]=0;
  env=getenv("OMP_NUM_THREADS");
  if(env) {
    iparmas[2]=atoi(env);}
  else{
    iparmas[2]=1;
  }

  printf(" number of threads =% d\n\n",iparmas[2]);

  for(i=0;i<64;i++){ptas[i]=0;}

  ndim=*neq+2**nzs;

  pointersas=NNEW(int,*neq+1);
  irowpardisoas=NNEW(int,ndim);
  aupardisoas=NNEW(double,ndim);

  k=0;

  /* PARDISO requires the matrix to be stored row by row
     aupardisoas contains the entries
     irowpardisoas the corresponding column numbers
         (for each row in ascending order)
     pointersas(i) points to the first entry for row i in 
           field aupardisoas */

  if(*sigma==0.){
    for(i=0;i<*neq;i++){
      pointersas[i]=k+1;
      
      /* lower left triangular matrix */

      for(j=0;j<i;j++){
	ifortran=i+1;
	lfortran=jq[j+1]-jq[j];
	FORTRAN(nident,(&irow[jq[j]-1],&ifortran,&lfortran,&id));
	if(id>0){
	  index=jq[j]+id-2;
	  if(irow[index]==ifortran){
	    irowpardisoas[k]=j+1;
	    aupardisoas[k]=au[index];
	    k++;
	  }
	}
      }

      /* diagonal entry */

      irowpardisoas[k]=i+1;
      aupardisoas[k]=ad[i];
      k++;

      /* upper right triangular matrix */

      for(j=jq[i];j<jq[i+1];j++){
	irowpardisoas[k]=irow[j-1];
	aupardisoas[k]=au[j+*nzs-1];
	k++;
      }
    }
    pointersas[*neq]=k+1;
  }else{
    for(i=0;i<*neq;i++){
      pointersas[i]=k+1;
      
      /* lower left triangular matrix */

      for(j=0;j<i;j++){
	ifortran=i+1;
	lfortran=jq[j+1]-jq[j];
	FORTRAN(nident,(&irow[jq[j]-1],&ifortran,&lfortran,&id));
	if(id>0){
	  index=jq[j]+id-2;
	  if(irow[index]==ifortran){
	    irowpardisoas[k]=j+1;
	    aupardisoas[k]=au[index]-*sigma*aub[index];
	    k++;
	  }
	}
      }

      /* diagonal entry */

      irowpardisoas[k]=i+1;
      aupardisoas[k]=ad[i]-*sigma*adb[i];
      k++;

      /* upper right triangular matrix */

      for(j=jq[i];j<jq[i+1];j++){
	irowpardisoas[k]=irow[j-1];
	aupardisoas[k]=au[j+*nzs-1]-*sigma*aub[j+*nzs-1];
	k++;
      }
    }
    pointersas[*neq]=k+1;
  }

  FORTRAN(pardiso,(ptas,&maxfct,&mnum,&mtype,&phase,neq,aupardisoas,
		   pointersas,irowpardisoas,perm,&nrhs,iparmas,&msglvl,
                   b,x,&error));

  return;
}

void pardiso_solve_as(double *b, int *neq){

  char *env;
  int maxfct=1,mnum=1,mtype=11,phase=33,*perm=NULL,nrhs=1,
      msglvl=0,i,error=0;
  double *x=NULL;

  printf(" Solving the system of equations using the asymmetric pardiso solver\n");

  iparmas[0]=0;
  env=getenv("OMP_NUM_THREADS");
  if(env) {
    iparmas[2]=atoi(env);}
  else{
    iparmas[2]=1;
  }

  printf(" number of threads =% d\n\n",iparmas[2]);

  x=NNEW(double,*neq);

  FORTRAN(pardiso,(ptas,&maxfct,&mnum,&mtype,&phase,neq,aupardisoas,
		   pointersas,irowpardisoas,perm,&nrhs,iparmas,&msglvl,
                   b,x,&error));

  for(i=0;i<*neq;i++){b[i]=x[i];}
  free(x);

  return;
}

void pardiso_cleanup_as(int *neq){

  int maxfct=1,mnum=1,mtype=11,phase=-1,*perm=NULL,nrhs=1,
      msglvl=0,error=0;
  double *b=NULL,*x=NULL;

  FORTRAN(pardiso,(ptas,&maxfct,&mnum,&mtype,&phase,neq,aupardisoas,
		   pointersas,irowpardisoas,perm,&nrhs,iparmas,&msglvl,
                   b,x,&error));

  free(irowpardisoas);
  free(aupardisoas);
  free(pointersas);

  return;
}

void pardiso_main_as(double *ad, double *au, double *adb, double *aub, double *sigma,
         double *b, int *icol, int *irow, 
         int *neq, int *nzs, int *jq){

  if(*neq==0) return;

  pardiso_factor_as(ad,au,adb,aub,sigma,icol,irow, 
             neq,nzs,jq);

  pardiso_solve_as(b,neq);

  pardiso_cleanup_as(neq);

  return;
}

#endif

