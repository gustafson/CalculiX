/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2014 Guido Dhondt                          */

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


void calcresidual_em(int *nmethod, int *neq, double *b, double *fext, double *f,
        int *iexpl, int *nactdof, double *aux1, double *aux2, double *vold,
        double *vini, double *dtime, double *accold, int *nk, double *adb,
        double *aub, int *jq, int *irow, int *nzl, double *alpha,
        double *fextini, double *fini, int *islavnode, int *nslavnode,
        int *mortar, int *ntie,double *f_cm,
	double* f_cs, int *mi,int *nzs,int *nasym,int *ithermal){

    int j,k,mt=mi[1]+1,jstart;
      
    /* residual for a static analysis */
      
    if(*nmethod!=4){
	for(k=0;k<neq[1];++k){
	    b[k]=fext[k]-f[k];
	}
    }
      
    /* residual for implicit dynamics */
      
    else{

	if(*ithermal<2){
	    jstart=1;
	}else{
	    jstart=0;
	}

        /* calculating a pseudo-velocity */

	for(k=0;k<*nk;++k){
	    for(j=jstart;j<mt;++j){
		if(nactdof[mt*k+j]!=0){aux2[nactdof[mt*k+j]-1]=(vold[mt*k+j]-vini[mt*k+j])/(*dtime);}
	    }
	}

        /* calculating "capacity"-matrix times pseudo-velocity */

	if(*nasym==0){
	    FORTRAN(op,(&neq[1],aux2,b,adb,aub,jq,irow)); 
	}else{
	    FORTRAN(opas,(&neq[1],aux2,b,adb,aub,jq,irow,nzs)); 
	}

	for(k=0;k<neq[1];++k){
	    b[k]=fext[k]-f[k]-b[k];
	} 
    }

    return;
}
