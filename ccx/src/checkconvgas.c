/*     CalculiX - A 3-dimensional finite element program                 */
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

void checkconvgas(int *icutb, int *iin,
		  double *qamt, double *qamf, double *qamp, 
		  double *ram1t, double *ram1f, double *ram1p,
		  double *ram2t, double *ram2f, double *ram2p,
		  double *ramt, double *ramf, double *ramp,
		  int *icntrl, double *dtheta, double *ctrl){
  
  int i0,ir,ip,ic,il,ig,ia,idivergence;
  
  double c1t,c1f,c1p;
  double df,dc,db,dd,ran,can,rap,ea,cae,ral;

  i0=ctrl[0];ir=ctrl[1];ip=ctrl[2];ic=ctrl[3];il=ctrl[4];ig=ctrl[5];ia=ctrl[7];
  df=ctrl[10];dc=ctrl[11];db=ctrl[12];dd=ctrl[16];
  ran=ctrl[18];can=ctrl[19];rap=ctrl[22];
  ea=ctrl[23];cae=ctrl[24];ral=ctrl[25];
  
  /* temperature */
  
  if(*iin<=ip){c1t=0.001*ran;}
  else{c1t=0.001*rap;}
  
  /* mass flow */
  
  if(*iin<=ip){c1f=0.001*ran;}
  else{c1f=0.001*rap;}
  
  /* pressure */
  
  if(*iin<=ip){c1p=0.001*ran;}
  else{c1p=0.001*rap;}
  
  if(*ram1t<*ram2t) {*ram2t=*ram1t;}
  if(*ram1f<*ram2f) {*ram2f=*ram1f;}
  if(*ram1p<*ram2p) {*ram2p=*ram1p;}
  
  /* check for convergence or divergence */

  if(((*ramt<=c1t**qamt)||(*ramt<1.e-8))&&
     ((*ramf<=c1f**qamf)||(*ramf<1.e-15))&&
     ((*ramp<=c1p**qamp)||(*ramp<1.e-8))&&
     (*iin!=0)){
      
      /* increment convergence reached */
      
      printf("      convergence\n\n");
      *icntrl=1;
      *icutb=0;
  }
  
  else {

      idivergence=0;

      /* divergence based on temperatures */
      
      if((*iin>=5*i0)||(fabs(*ramt)>1.e20)){
	  if((*ram1t>*ram2t)&&(*ramt>*ram2t)&&(*ramt>c1t**qamt)){
	      idivergence=1;
	  }
      }

      /* divergence based on the mass flux */
      
      if((*iin>=5*i0)||(fabs(*ramt)>1.e20)){
	  if((*ram1t>*ram2t)&&(*ramt>*ram2t)&&(*ramt>c1t**qamt)){
	      idivergence=1;
	  }
      }

      /* divergence based on pressures */
      
      if((*iin>=5*i0)||(fabs(*ramt)>1.e20)){
	  if((*ram1t>*ram2t)&&(*ramt>*ram2t)&&(*ramt>c1t**qamt)){
	      idivergence=1;
	  }
      }

      /* divergence based on singular matrix or negative pressures */

      if(*iin==0) idivergence=1;
      
      if(idivergence==1){
	  *dtheta=*dtheta*df;
	  printf("\n divergence; the increment size is decreased to %e\n",*dtheta);
	  printf(" the increment is reattempted\n\n");
	  *iin=0;
	  (*icutb)++;
	  if(*icutb>ia){
	      printf("\n *ERROR: too many cutbacks\n");
	      FORTRAN(stop,());
	  }
      }else{
	  printf("      no convergence\n\n");
      }
  }
  return;
}
