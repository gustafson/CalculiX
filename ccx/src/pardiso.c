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

int *irowpardiso=NULL,*pointers=NULL,iparm[64];
long int pt[64];
double *aupardiso=NULL;

void pardiso_factor(double *ad, double *au, double *adb, double *aub, 
                double *sigma,int *icol, int *irow, 
                int *neq, int *nzs){

  char *env;
  int i,j,k,l,maxfct=1,mnum=1,mtype=-2,phase=12,nrhs=1,*perm=NULL,
      msglvl=0,error=0;
  long long ndim;
  double *b=NULL,*x=NULL;

  printf(" Factoring the system of equations using the symmetric pardiso solver\n");

  iparm[0]=0;
  env=getenv("OMP_NUM_THREADS");
  if(env) {
    iparm[2]=atoi(env);}
  else{
    iparm[2]=1;
  }

  printf(" number of threads =% d\n\n",iparm[2]);

  for(i=0;i<64;i++){pt[i]=0;}

  ndim=*neq+*nzs;

  pointers=NNEW(int,*neq+1);
  irowpardiso=NNEW(int,ndim);
  aupardiso=NNEW(double,ndim);

  k=ndim;
  l=*nzs;

  if(*sigma==0.){
    pointers[*neq]=ndim+1;
    for(i=*neq-1;i>=0;--i){
      for(j=0;j<icol[i];++j){
	irowpardiso[--k]=irow[--l];
	aupardiso[k]=au[l];
      }
      pointers[i]=k--;
      irowpardiso[k]=i+1;
      aupardiso[k]=ad[i];
    }
  }
  else{
    pointers[*neq]=ndim+1;
    for(i=*neq-1;i>=0;--i){
      for(j=0;j<icol[i];++j){
	irowpardiso[--k]=irow[--l];
	aupardiso[k]=au[l]-*sigma*aub[l];
      }
      pointers[i]=k--;
      irowpardiso[k]=i+1;
      aupardiso[k]=ad[i]-*sigma*adb[i];
    }
  }

  /*  for(i=0;i<ndim;i++){printf("aupardiso %d %lf\n",i+1,aupardiso[i]);}
  for(i=0;i<ndim;i++){printf("irowpardiso %d %d\n",i+1,irowpardiso[i]);}
  for(i=0;i<*neq+1;i++){printf("pointers %d %d\n",i+1,pointers[i]);}
  printf("maxfct,mnum,mtype,phase,neq,nrhs,msglvl,error=%d,%d,%d,%d,%d,%d,%d,%d\n",maxfct,mnum,mtype,phase,*neq,nrhs,msglvl,error);

  printf(" pardiso_factor: before pardiso\n\n");*/

  FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
		   pointers,irowpardiso,perm,&nrhs,iparm,&msglvl,
                   b,x,&error));

  /*printf(" pardiso_factor: after pardiso\n\n");*/

  return;
}

void pardiso_solve(double *b, int *neq){

  int maxfct=1,mnum=1,mtype=-2,phase=33,*perm=NULL,nrhs=1,
      msglvl=0,i,error=0;
  double *x=NULL;

  printf(" Solving the system of equations using the symmetric pardiso solver\n");

  printf(" number of threads =% d\n\n",iparm[2]);

  x=NNEW(double,*neq);

  /*printf(" pardiso_solve: before pardiso\n\n");
    for(i=0;i<*neq;i++){printf("b %d %lf\n",i+1,b[i]);}*/

  FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
		   pointers,irowpardiso,perm,&nrhs,iparm,&msglvl,
                   b,x,&error));

  /* printf(" pardiso_solve: after pardiso\n\n");
     for(i=0;i<*neq;i++){printf("x %d %lf\n",i+1,x[i]);}*/


  for(i=0;i<*neq;i++){b[i]=x[i];}
  free(x);

  return;
}

void pardiso_cleanup(int *neq){

  int maxfct=1,mnum=1,mtype=-2,phase=-1,*perm=NULL,nrhs=1,
      msglvl=0,error=0;
  double *b=NULL,*x=NULL;

  /* printf(" pardiso_clean: before pardiso\n\n");*/


  FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
		   pointers,irowpardiso,perm,&nrhs,iparm,&msglvl,
                   b,x,&error));

  /*printf(" pardiso_clean: after pardiso\n\n");*/


  free(irowpardiso);
  free(aupardiso);
  free(pointers);

  return;
}

void pardiso_main(double *ad, double *au, double *adb, double *aub, double *sigma,
         double *b, int *icol, int *irow, 
         int *neq, int *nzs){

  if(*neq==0) return;

  pardiso_factor(ad,au,adb,aub,sigma,icol,irow, 
             neq,nzs);

  pardiso_solve(b,neq);

  pardiso_cleanup(neq);


  /* printf(" pardiso_main: geschaft\n\n");*/

  

  return;
}

#endif

