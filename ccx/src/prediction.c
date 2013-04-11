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


void prediction(double *uam, int *nmethod, double *bet, double *gam, double *dtime,
               int *ithermal, int *nk, double *veold, double *accold, double *v,
               int *iinc, int *idiscon, double *vold, int *nactdof){

    int j,k;
    double dextrapol,scal1,scal2;

    uam[0]=0.;
    uam[1]=0.;
    if(*nmethod==4){
      
	scal1=0.5*(1.-2.**bet)**dtime**dtime;
	scal2=(1.-*gam)**dtime;
	
	if(*ithermal<2){
	    for(k=0;k<*nk;++k){
		for(j=0;j<4;j++){
		    dextrapol=*dtime*veold[4*k+j]+scal1*accold[4*k+j];
		    if((fabs(dextrapol)>uam[0])&&(nactdof[4*k+j]>0)) {uam[0]=fabs(dextrapol);}
		    v[5*k+j]=vold[5*k+j]+dextrapol;
		    veold[4*k+j]=veold[4*k+j]+scal2*accold[4*k+j];
		    accold[4*k+j]=0.;
		}
	    }
	}else if(*ithermal==2){
	    for(k=0;k<*nk;++k){
		for(j=0;j<4;j++){
		    v[5*k+j]=vold[5*k+j];
		}
	    }
	    for(k=0;k<*nk;++k){
		dextrapol=*dtime*veold[4*k];
		if(fabs(dextrapol)>100.) dextrapol=100.*dextrapol/fabs(dextrapol);
		if((fabs(dextrapol)>uam[1])&&(nactdof[4*k]>0)) {uam[1]=fabs(dextrapol);}
		v[5*k]+=dextrapol;
	    }
	}else{
	    for(k=0;k<*nk;++k){
		for(j=0;j<4;++j){
		    dextrapol=*dtime*veold[4*k+j]+scal1*accold[4*k+j];
		    if((j==0)&&fabs(dextrapol)>100.) dextrapol=100.*dextrapol/fabs(dextrapol);
		    if(j==0){
			if((fabs(dextrapol)>uam[1])&&(nactdof[4*k]>0)) {uam[1]=fabs(dextrapol);}
		    }else{
			if((fabs(dextrapol)>uam[0])&&(nactdof[4*k+j]>0)) {uam[0]=fabs(dextrapol);}
		    }
		    v[5*k+j]=vold[5*k+j]+dextrapol;
		    veold[4*k+j]=veold[4*k+j]+scal2*accold[4*k+j];
		    accold[4*k+j]=0.;
		}
	    }
	}
    }
    
    /* for the static case: extrapolation of the previous increment
       (if any within the same step) */
    
    else{
	if(*iinc>1){
	    if(*ithermal<2){
		for(k=0;k<*nk;++k){
		    for(j=0;j<4;++j){
			if(*idiscon==0){
			    dextrapol=*dtime*veold[4*k+j];
			    if((fabs(dextrapol)>uam[0])&&(nactdof[4*k+j]>0)) {uam[0]=fabs(dextrapol);}	
			    v[5*k+j]=vold[5*k+j]+dextrapol;
			}else{
			    v[5*k+j]=vold[5*k+j];
			}
		    }
		}
	    }else if(*ithermal==2){
		for(k=0;k<*nk;++k){
		    for(j=0;j<4;++j){
			v[5*k+j]=vold[5*k+j];
		    }
		}
		for(k=0;k<*nk;++k){
		    if(*idiscon==0){
			dextrapol=*dtime*veold[4*k];
			if(fabs(dextrapol)>100.) dextrapol=100.*dextrapol/fabs(dextrapol);
			if((fabs(dextrapol)>uam[1])&&(nactdof[4*k]>0)) {uam[1]=fabs(dextrapol);}	
			v[5*k]+=dextrapol;
		    }
		}
	    }else{
		for(k=0;k<*nk;++k){
		    for(j=0;j<4;++j){
			if(*idiscon==0){
			    dextrapol=*dtime*veold[4*k+j];
			    if((j==0)&&fabs(dextrapol)>100.) dextrapol=100.*dextrapol/fabs(dextrapol);
			    if(j==0){
				if((fabs(dextrapol)>uam[1])&&(nactdof[4*k+j]>0)) {uam[1]=fabs(dextrapol);}
			    }else{
				if((fabs(dextrapol)>uam[0])&&(nactdof[4*k+j]>0)) {uam[0]=fabs(dextrapol);}
			    }	
			    v[5*k+j]=vold[5*k+j]+dextrapol;
			}else{
			    v[5*k+j]=vold[5*k+j];
			}
		    }
		}
	    }
	}
	else{
	    for(k=0;k<*nk;++k){
		for(j=0;j<4;++j){
		    v[5*k+j]=vold[5*k+j];
		}
	    }
	}
    }
    *idiscon=0;

  return;
}
