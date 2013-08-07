/*     CalculiX - A 3-dimensional finite element program                 */
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


void calcresidual(int *nmethod, int *neq, double *b, double *fext, double *f,
        int *iexpl, int *nactdof, double *aux1, double *aux2, double *vold,
        double *vini, double *dtime, double *accold, int *nk, double *adb,
        double *aub, int *icol, int *irow, int *nzl, double *alpha,
        double *fextini, double *fini, int *islavnode, int *nslavnode,
        int *mortar, int *ntie,double *f_cm,
	double* f_cs, int *mi,int *nzs,int *nasym){

    int j,k,mt=mi[1]+1;
    double scal1;
      
    /* residual for a static analysis */
      
    if(*nmethod!=4){
	for(k=0;k<neq[1];++k){
	    b[k]=fext[k]-f[k];
//	   printf("calcresidual dof=%d,fext=%e,f=%e,resi=%e\n",k,fext[k],f[k],b[k]);
	}
    }
      
    /* residual for implicit dynamics */
      
    else if(*iexpl<=1){
	for(k=0;k<*nk;++k){
	    if(nactdof[mt*k]!=0){
		aux2[nactdof[mt*k]-1]=(vold[mt*k]-vini[mt*k])/(*dtime);}
	    for(j=1;j<mt;++j){
		if(nactdof[mt*k+j]!=0){aux2[nactdof[mt*k+j]-1]=accold[mt*k+j];}
	    }
	}
	if(*nasym==0){
	    FORTRAN(op,(&neq[1],aux1,aux2,b,adb,aub,icol,irow,nzl)); 
	}else{
	    FORTRAN(opas,(&neq[1],aux1,aux2,b,adb,aub,icol,irow,nzl,nzs)); 
	}
	scal1=1.+*alpha;
	for(k=0;k<neq[0];++k){
	    b[k]=scal1*(fext[k]-f[k])-*alpha*(fextini[k]-fini[k])-b[k];
	} 
	for(k=neq[0];k<neq[1];++k){
	    b[k]=fext[k]-f[k]-b[k];
	} 
    }
    
    /* residual for explicit dynamics */
    
    else{
	for(k=0;k<*nk;++k){
	    if(nactdof[mt*k]!=0){
		aux2[nactdof[mt*k]-1]=(vold[mt*k]-vini[mt*k])/(*dtime);}
	    for(j=1;j<mt;++j){
		if(nactdof[mt*k+j]!=0){aux2[nactdof[mt*k+j]-1]=accold[mt*k+j];}
	    }
	}
	scal1=1.+*alpha;
	for(k=0;k<neq[0];++k){
	    b[k]=scal1*(fext[k]-f[k])-*alpha*(fextini[k]-fini[k])
		-adb[k]*aux2[k];
	} 
	for(k=neq[0];k<neq[1];++k){
	    b[k]=fext[k]-f[k]-adb[k]*aux2[k];
	} 
    }
    
    if(*mortar==1){
      /* removing the forces in the slave nodes */
      for(k=0;k<neq[1];++k){
	b[k]-=f_cs[k];
	b[k]-=f_cm[k];
      }	
    }

    return;
}
