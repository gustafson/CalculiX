/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2013 Guido Dhondt                          */

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

#ifdef TAUCS

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"
#include "tau.h"
#include <taucs.h>

taucs_ccs_matrix aa[1];
void* F=NULL;
char* taufactor[]={ "taucs.factor.LLT=true","taucs.factor.mf=true",
  "taucs.factor.ordering=amd",NULL };
char* taufactorooc[]={ "taucs.factor.LLT=true","taucs.ooc=true",
                  "taucs.ooc.basename=/home/guido/scratch/scratch",
                  "taucs.ooc.memory=500000000.0",NULL };
char* tausolve[]={ "taucs.factor=false",NULL };
char* tausolveooc[]={"taucs.factor=false",NULL };
int *irowtau=NULL,*pointtau=NULL;
double *autau=NULL;
int* perm;


void tau_factor(double *ad, double **aup, double *adb, double *aub, 
                double *sigma,int *icol, int **irowp, 
                int *neq, int *nzs){

  int i,j,k,l,*irow=NULL;
  long long ndim;
  double *au=NULL;
  double memory_mb = -1.0;
  int    mb = -1;
  int ret;

  printf(" Factoring the system of equations using TAUCS\n\n");

  taucs_logfile("stdout");

  au=*aup;
  irow=*irowp;

  ndim=*neq+*nzs;

  autau= NNEW(double,ndim);
  irowtau=NNEW(int,ndim);
  pointtau=NNEW(int,*neq+1);

  k=ndim;
  l=*nzs;

  if(*sigma==0.){
    pointtau[*neq]=ndim;
    for(i=*neq-1;i>=0;--i){
      for(j=0;j<icol[i];++j){
	irowtau[--k]=irow[--l]-1;
	autau[k]=au[l];
      }
      pointtau[i]=--k;
      irowtau[k]=i;
      autau[k]=ad[i];
    }
  }
  else{
    pointtau[*neq]=ndim;
    for(i=*neq-1;i>=0;--i){
      for(j=0;j<icol[i];++j){
	irowtau[--k]=irow[--l]-1;
	autau[k]=au[l]-*sigma*aub[l];
      }
      pointtau[i]=--k;
      irowtau[k]=i;
      autau[k]=ad[i]-*sigma*adb[i];
    }
  }

  /* defining the matrix */

  aa->n = *neq;
  aa->m = *neq;
  aa->flags  = TAUCS_SYMMETRIC | TAUCS_LOWER | TAUCS_DOUBLE;
  aa->colptr = pointtau;
  aa->rowind = irowtau;
  aa->values.d = autau;

  if(*neq<50000){
      taucs_linsolve(aa,&F,0,NULL,NULL,taufactor,NULL);
  }
  else{
      /*ret = taucs_linsolve(aa,&F,0,NULL,NULL,taufactorooc,NULL);*/

     if (mb > 0)
        memory_mb = (double) mb;
      else
        memory_mb = ((double) (-mb)) * taucs_available_memory_size()/1048576.0;  
	            
      F = taucs_io_create_multifile("~/scratch/scratch");
        
      ret = taucs_ooc_factor_llt(aa,F,memory_mb*1048576.0);
            
      printf(" Return Code from Factoring %d\n\n",ret);
  }
  
  *aup=au;
  *irowp=irow;

  return;
}

void tau_solve(double *b,int *neq){

  int i;
  /*static int j;*/
  double *x=NULL;
  int ret;

  x=NNEW(double,*neq);
  
  if(*neq<150){
      taucs_linsolve(aa,&F,1,x,b,tausolve,NULL);
  }
  else{
      /*ret = taucs_linsolve(aa,&F,1,x,b,tausolveooc,NULL);*/
      
      ret = taucs_ooc_solve_llt(F, x, b);
       
      printf(" Return Code from Solving %d\n\n",ret);
      
      taucs_io_delete(F);
  }

  for(i=0;i<=*neq-1;++i){
    b[i]=x[i];
  }
  free(x);/*
  if (mb > 0)
    memory_mb = (double) mb;
  else
    memory_mb = ((double) (-mb)) * taucs_available_memory_size()/1048576.0;  
  */
  /*j++;printf("%d\n",j);*/

  return;
}

void tau_cleanup(){

  /*taucs_linsolve(NULL,&F,0,NULL,NULL,NULL,NULL);*/
  free(pointtau);
  free(irowtau);
  free(autau);

  return;
}

void tau(double *ad, double **aup, double *adb, double *aub, double *sigma,
         double *b, int *icol, int **irowp, 
         int *neq, int *nzs){

  if(*neq==0) return;

     
  tau_factor(ad,aup,adb,aub,sigma,icol,irowp, 
             neq,nzs);

  tau_solve(b,neq);
  
  tau_cleanup();
  

  return;
}

#endif
