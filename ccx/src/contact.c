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

void contact(int *ncont, int *ntie, char *tieset,int *nset,char *set,
	     int *istartset, int *iendset, int *ialset,int *itietri,
	     char *lakon, int *ipkon, int *kon, int *koncont, int *ne,
	     double *cg, double *straight, int *ifree, double *co,
	     double *vold, int *ielmat, double *cs, double *elcon,
             int *istep,int *iinc,int *iit,int *ncmat_,int *ntmat_,
             int *ne0, double *vini,
             int *nmethod, int *nmpc, int *mpcfree, int *memmpc_,
             int **ipompcp, char **labmpcp, int **ikmpcp, int **ilmpcp,
             double **fmpcp, int **nodempcp, double **coefmpcp,
             int *iperturb, int *ikboun, int *nboun, int *mi,
             int *imastop,int *nslavnode,int *islavnode,int *islavsurf,
             int *itiefac,double *areaslav,int *iponoels,int *inoels,
             double *springarea, double *tietol, double *reltime,
	     int *imastnode, int *nmastnode, double *xmastnor,
	     double *xnormastface, char *filab, int *mcs, int *ics,
             int *nasym,double *xnoels){
    
    char *labmpc=NULL;

    int i,ntrimax,*nx=NULL,*ny=NULL,*nz=NULL,*ipompc=NULL,*ikmpc=NULL,
	*ilmpc=NULL,*nodempc=NULL,nmpc_,im;
    
    double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL,
        *fmpc=NULL, *coefmpc=NULL;

    ipompc=*ipompcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
    fmpc=*fmpcp;nodempc=*nodempcp;coefmpc=*coefmpcp;
    nmpc_=*nmpc;

    DMEMSET(xmastnor,0,3*nmastnode[*ntie],0.);
    
    FORTRAN(updatecontpen,(koncont,ncont,co,vold,
			cg,straight,mi,imastnode,nmastnode,xmastnor,
			ntie,tieset,nset,set,istartset,
			iendset,ialset,ipkon,lakon,kon,cs,mcs,ics));
    
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
    
    FORTRAN(gencontelem,(tieset,ntie,itietri,ne,ipkon,kon,lakon,
       cg,straight,ifree,koncont,
       co,vold,xo,yo,zo,x,y,z,nx,ny,nz,ielmat,cs,elcon,istep,
       iinc,iit,ncmat_,ntmat_,ne0,vini,nmethod,mi,
       imastop,nslavnode,islavnode,islavsurf,itiefac,areaslav,iponoels,
       inoels,springarea,ikmpc,ilmpc,nmpc,ipompc,nodempc,coefmpc,
       set,nset,istartset,iendset,ialset,tietol,reltime,xmastnor,
       xnormastface,imastnode,nmastnode,filab,nasym,xnoels));

    free(xo);free(yo);free(zo);free(x);free(y);free(z);free(nx);
    free(ny);free(nz);

    *ipompcp=ipompc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
    *fmpcp=fmpc;*nodempcp=nodempc;*coefmpcp=coefmpc;
  
    return;
}
