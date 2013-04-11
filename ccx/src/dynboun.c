/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2007 Guido Dhondt                          */

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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

#ifdef SPOOLES
   #include "spooles.h"
#endif
#ifdef SGI
   #include "sgi.h"
#endif
#ifdef TAUCS
   #include "tau.h"
#endif

void dynboun(double *amta,int *namta,int *nam,double *ampli, double *time,
             double *ttime,double *dtime,double *xbounold,double *xboun,
             double *xbounact,int *iamboun,int *nboun,int *nodeboun,
             int *ndirboun, double *ad, double *au, double *adb,
             double *aub, int *icol, int *irow, int *neq, int *nzs,
             double *sigma, double *b, int *isolver,
             double *alpham, double *betam, int *nzl,
             int *init,double *bact, double *bmin, int *jq){

    int idiff[3],i,j,ic,ir;

    double *xbounmin=NULL,*xbounplus=NULL,*bplus=NULL,
	*bv=NULL,*ba=NULL,deltatime,deltatime2,deltatimesq,timemin,ttimemin,
        timeplus,ttimeplus,*aux=NULL,*b1=NULL,*b2=NULL;

#ifdef SGI
  int token=1;
#endif
    
  xbounmin=NNEW(double,*nboun);
  xbounplus=NNEW(double,*nboun);

      /* time increment for the calculation of the change of the
         particular solution (needed to account for nonzero
         SPC's) */

  deltatime=*dtime;
  deltatime2=2.*deltatime;
  deltatimesq=deltatime*deltatime;

      /* the SPC value at timemin is stored in xbounmin */

  if(*init==1){

      /* at the start of a new step it is assumed that the previous step
	 has reached steady state (at least for the SPC conditions) */

      for(i=0;i<*nboun;i++){
	  xbounmin[i]=xbounold[i];
	  xbounact[i]=xbounold[i];
      }
  }
  else{
      timemin=*time-deltatime;
      ttimemin=*ttime-deltatime;
      FORTRAN(temploadmodal,(amta,namta,nam,ampli,&timemin,&ttimemin,dtime,
	   xbounold,xboun,xbounmin,iamboun,nboun,nodeboun,ndirboun));
  }

      /* the SPC value at timeplus is stored in xbounplus */

  timeplus=*time+deltatime;
  ttimeplus=*ttime+deltatime;
  FORTRAN(temploadmodal,(amta,namta,nam,ampli,&timeplus,&ttimeplus,dtime,
	  xbounold,xboun,xbounplus,iamboun,nboun,nodeboun,ndirboun));

  bplus=NNEW(double,neq[1]);
  bv=NNEW(double,neq[1]);
  ba=NNEW(double,neq[1]);
  b1=NNEW(double,neq[1]);
  b2=NNEW(double,neq[1]);

      /* check whether boundary conditions changed 
         comparision of min with prev */

  if(*init==1){
      for(i=0;i<*nboun;i++){
	  ic=neq[1]+i;
	  for(j=jq[ic]-1;j<jq[ic+1]-1;j++){
	      ir=irow[j]-1;
	      bmin[ir]=bmin[ir]-au[j]*xbounmin[i];
	  }
      }
      if(*isolver==0){
#ifdef SPOOLES
	  spooles_solve(bmin,&neq[1]);
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	  sgi_solve(bmin,token);
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	  tau_solve(bmin,&neq[1]);
#endif
      }
  }

  /* check whether boundary conditions changed 
     comparision of act with min */
  
  idiff[1]=0;
  for(i=0;i<*nboun;i++){
      if(fabs(xbounact[i]-xbounmin[i])>1.e-10){
	  idiff[1]=1;
	  break;
      }
  }
  if(*init==1){
      for(i=0;i<*nboun;i++){
	  ic=neq[1]+i;
	  for(j=jq[ic]-1;j<jq[ic+1]-1;j++){
	      ir=irow[j]-1;
	      bact[ir]=bact[ir]-au[j]*xbounact[i];
	  }
      }
      if(*isolver==0){
#ifdef SPOOLES
	  spooles_solve(bact,&neq[1]);
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	  sgi_solve(bact,token);
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	  tau_solve(bact,&neq[1]);
#endif
      }
  }

      /* check whether boundary conditions changed 
         comparision of plus with act */
  
  idiff[2]=0;
  for(i=0;i<*nboun;i++){
      if(fabs(xbounplus[i]-xbounact[i])>1.e-10){
	  idiff[2]=1;
	  break;
      }
  }
  if(idiff[2]==1){
      for(i=0;i<*nboun;i++){
	  ic=neq[1]+i;
	  for(j=jq[ic]-1;j<jq[ic+1]-1;j++){
	      ir=irow[j]-1;
	      bplus[ir]=bplus[ir]-au[j]*xbounplus[i];
	  }
      }
      if(*isolver==0){
#ifdef SPOOLES
	  spooles_solve(bplus,&neq[1]);
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	  sgi_solve(bplus,token);
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	  tau_solve(bplus,&neq[1]);
#endif
      }
  }
  
  if((idiff[1]!=0)||(idiff[2]!=0)){
      if(idiff[2]==0){
	  for(i=0;i<neq[1];i++){bplus[i]=bact[i];}
      }
      for(i=0;i<neq[1];i++){
	  
	  /* bv is the velocity */
	  
	  bv[i]=(bplus[i]-bmin[i])/deltatime2;
	  
	  /* ba is the acceleration */
	  
	  ba[i]=(bmin[i]-2.*bact[i]+bplus[i])/deltatimesq;
	  
	  b1[i]=ba[i]+*alpham*bv[i];
	  b2[i]=*betam*bv[i];

	  bmin[i]=bact[i];
	  bact[i]=bplus[i];
      }
      FORTRAN(op,(&neq[1],aux,b1,bplus,adb,aub,
		  icol,irow,nzl));
      for(i=0;i<neq[1];i++){b[i]=b[i]-bplus[i];}
      FORTRAN(op,(&neq[1],aux,b2,bplus,ad,au,
		  icol,irow,nzl));
      for(i=0;i<neq[1];i++){b[i]=b[i]-bplus[i];}
  }
  
  free(xbounmin);free(xbounplus);
  free(bplus);free(bv);free(ba);free(b1);free(b2);
  
  return;
}
