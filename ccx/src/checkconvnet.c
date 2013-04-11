/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                          */

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

void checkconvnet(int *icutb, int *iin,
		  double *uamt, double *uamf, double *uamp, 
		  double *cam1t, double *cam1f, double *cam1p,
		  double *cam2t, double *cam2f, double *cam2p,
		  double *camt, double *camf, double *camp,
		  int *icntrl, double *dtheta, double *ctrl,
                  double *uama,double *cam1a,double *cam2a,double *cama,
                  double *vamt, double *vamf, double *vamp, double *vama){
  
  int i0,ir,ip,ic,il,ig,ia,idivergence;
  
  double c1t,c1f,c1p,c1a;
  double df,dc,db,dd,ran,can,rap,ea,cae,ral;

  i0=ctrl[0];ir=ctrl[1];ip=ctrl[2];ic=ctrl[3];il=ctrl[4];ig=ctrl[5];ia=ctrl[7];
  df=ctrl[10];dc=ctrl[11];db=ctrl[12];dd=ctrl[16];
  ran=ctrl[18];can=ctrl[19];rap=ctrl[22];
  ea=ctrl[23];cae=ctrl[24];ral=ctrl[25];
  
  /* temperature */
  
  if(*iin<=ip){c1t=0.0001*ran;}
  else{c1t=0.0001*rap;}
  
  /* mass flow */
  
  if(*iin<=ip){c1f=0.0001*ran;}
  else{c1f=0.0001*rap;}
  
  /* pressure */
  
  if(*iin<=ip){c1p=0.0001*ran;}
  else{c1p=0.0001*rap;}
  
  /* geometry */
  
  if(*iin<=ip){c1a=0.0001*ran;}
  else{c1a=0.0001*rap;}
  
  if(*cam1t<*cam2t) {*cam2t=*cam1t;}
  if(*cam1f<*cam2f) {*cam2f=*cam1f;}
  if(*cam1p<*cam2p) {*cam2p=*cam1p;}
  if(*cam1a<*cam2a) {*cam2a=*cam1a;}
  
  /* check for convergence or divergence 
     comparison of the latest change with
        - the largest change in the present calculation
        - the largest value in the present calculation */

  if(((*camt<=c1t**uamt)||(*camt<1.e-8**vamt))&&
     ((*camf<=c1f**uamf)||(*camf<1.e-8**vamf))&&
     ((*camp<=c1p**uamp)||(*camp<1.e-8**vamp))&&
     ((*cama<=c1p**uama)||(*cama<1.e-8**vama))&&
     (*iin>3)){
      
      /* increment convergence reached */
      
      printf("      flow network: convergence in gas iteration %d \n\n",*iin);
      *icntrl=1;
      *icutb=0;
  }
  
  else {

      idivergence=0;

      /* divergence based on temperatures */
      
      if((*iin>=20*i0)||(fabs(*camt)>1.e20)){
	  if((*cam1t>=*cam2t)&&(*camt>=*cam2t)&&(*camt>c1t**uamt)){
	      idivergence=1;
	  }
      }

      /* divergence based on the mass flux */
      
      if((*iin>=20*i0)||(fabs(*camf)>1.e20)){
	  if((*cam1f>=*cam2f)&&(*camf>=*cam2f)&&(*camf>c1f**uamf)){
	      idivergence=1;
	  }
      }

      /* divergence based on pressures */
      
      if((*iin>=20*i0)||(fabs(*camp)>1.e20)){
	  if((*cam1p>=*cam2p)&&(*camp>=*cam2p)&&(*camp>c1p**uamp)){
	      idivergence=1;
	  }
      }

      /* divergence based on geometry */
      
      if((*iin>=20*i0)||(fabs(*cama)>1.e20)){
	  if((*cam1a>=*cam2a)&&(*cama>=*cam2a)&&(*cama>c1p**uama)){
	      idivergence=1;
	  }
      }

      /* divergence based on the number of iterations */

      if(*iin>20*ic) idivergence=1;

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
