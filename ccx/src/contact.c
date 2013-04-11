/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2007 Guido Dhondt                     */

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

void contact(int *ncont, int *ntie, char *tieset, int *nset, char *set,
	     int *istartset, int *iendset, int *ialset, int *itietri,
	     char *lakon, int *ipkon, int *kon, int *koncont, int *ne,
	     double *cg, double *straight, int *ifree, double *co,
	     double *vold, int *ielmat, double *cs, double *elcon,
             int *istep,int *iinc,int *iit,int *ncmat_,int *ntmat_,
             int *ifcont1, int *ifcont2, int *ne0, double *vini,
             int *nmethod){
    
    int i,ntrimax,*nx=NULL,*ny=NULL,*nz=NULL;
    
    double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL;
    
    FORTRAN(updatecont,(koncont,ncont,co,vold,
			cg,straight));
    
    /* determining the size of the auxiliary fields */
    
    ntrimax=0;
    for(i=0;i<*ntie;i++){
	if(itietri[2*i+1]-itietri[2*i]+1>ntrimax)
	    ntrimax=itietri[2*i+1]-itietri[2*i]+1;
    }
    xo=NNEW(double,ntrimax);
    yo=NNEW(double,ntrimax);
    zo=NNEW(double,ntrimax);
    x=NNEW(double,ntrimax);
    y=NNEW(double,ntrimax);
    z=NNEW(double,ntrimax);
    nx=NNEW(int,ntrimax);
    ny=NNEW(int,ntrimax);
    nz=NNEW(int,ntrimax);
    
    FORTRAN(gencontelem,(tieset,ntie,itietri,ne,ipkon,kon,lakon,set,
       istartset,iendset,ialset,cg,straight,ifree,koncont,
       co,vold,xo,yo,zo,x,y,z,nx,ny,nz,nset,ielmat,cs,elcon,istep,
       iinc,iit,ncmat_,ntmat_,ifcont1,ifcont2,ne0,vini,nmethod));

    free(xo);free(yo);free(zo);free(x);free(y);free(z);free(nx);
    free(ny);free(nz);
  
    return;
}
