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

void remeshcontactq(int *ntie, char *tieset, int *nset, char *set,
               int *istartset, int *iendset, int **ialsetp,
               char **lakonp, int **ipkonp, int **konp,
	       int *nalset,int *nmpc, int *mpcfree, int *memmpc_,
               int **ipompcp, char **labmpcp, int **ikmpcp, int **ilmpcp,
               double **fmpcp, int **nodempcp, double **coefmpcp,
	       double **cop,int *nmpc_, int *mi, int *nk, int *nkon, 
	       int *ne,int *nk_, int *ithermal, int **ielmatp, 
	       int **ielorienp,double **t0p,double **voldp,double **veoldp,
	       int *ncont,double **xstatep,int *nstate_,double **prestrp,
	       int *iprestr,int *nxstate,int **iamt1p){

  char *labmpc=NULL,*lakon=NULL;

  int k,index,*ipompc=NULL,*ikmpc=NULL,*nodface=NULL,im,
      *ilmpc=NULL,*nodempc=NULL,neqterms,*ipoface=NULL,nface,
      nquadface,ninterface,*ipkon=NULL,*ielmat=NULL,
      *kon=NULL,*ialset=NULL,nk0,*ielorien=NULL,mt=mi[1]+1,
      ntets2remesh,*iamt1=NULL;

  double *fmpc=NULL,*coefmpc=NULL,*co=NULL,*t0=NULL,
      *vold=NULL,*veold=NULL,*xstate=NULL,*prestr=NULL;

  ipompc=*ipompcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  fmpc=*fmpcp;nodempc=*nodempcp;coefmpc=*coefmpcp;ipkon=*ipkonp;
  kon=*konp;co=*cop;lakon=*lakonp;ialset=*ialsetp;ielmat=*ielmatp;
  ielorien=*ielorienp;t0=*t0p;vold=*voldp;veold=*veoldp;xstate=*xstatep;
  prestr=*prestrp;iamt1=*iamt1p;
      
  /* allocating the field to catalogue the faces external to the 
     set A of all C3D20(R) elements adjacent to contact surfaces
     for field nodface: 9 entries per face, 6 faces per element */
  
  nodface=NNEW(int,9*6**ne);
  ipoface=NNEW(int,*nk);

  /* allocating fields for listing all C3D10 elements belonging to
     one and the same end node */
  
  FORTRAN(remeshsurfq,(tieset,ntie,set,nset,istartset,
		      iendset,ialset,ipkon,kon,lakon,nodface,ipoface,&nface,
		      nk,ne));

  /* if no quadratic faces: return */

  if(nface==0){
      free(nodface);free(ipoface);
      *ncont=0;
      return;
  }
  
  RENEW(nodface,int,9*nface);
      
  /* one extra node per quadratic contact face
     (in the middle of the face) */
      
  RENEW(co,double,3*(*nk+nface));
  
  RENEW(t0,double,*nk+nface);
  RENEW(iamt1,int,*nk+nface);
  for(k=*nk_;k<*nk+nface;k++){iamt1[k]=0;}
  
  RENEW(vold,double,mt*(*nk+nface));
  DMEMSET(vold,mt**nk,mt*(*nk+nface),0.);
  
  RENEW(veold,double,mt*(*nk+nface));

  /* each element adjacent to a quadratic contact
     face gets more nodes:
     C3D20(R) -> C3D26(R)
     C3D10    -> C3D14
     C3D15    -> C3D19
     they are moved to a new location in field kon
     without freeing the previous space
 */

  RENEW(kon,int,*nkon+26*nface);
  
  nk0=*nk;

  FORTRAN(remeshcontactelq,(tieset,ntie,set,nset,istartset,
	  iendset,ialset,ipkon,kon,nkon,lakon,nodface,ipoface,
	  nk,ipompc,nodempc,ikmpc,ilmpc,nmpc,nmpc_,labmpc,coefmpc,
	  mpcfree,nalset,co,ithermal,&nk0,ne,ielmat,ielorien,mi,
	  t0,vold,veold,xstate,nstate_,prestr,iprestr));

  free(nodface);free(ipoface);
  
  if(nface>0){
      *nk_+=(*nk-nk0);
  }
  
  *ipompcp=ipompc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*nodempcp=nodempc;*coefmpcp=coefmpc;*ipkonp=ipkon;
  *konp=kon;*lakonp=lakon;*cop=co;*ialsetp=ialset;*ielmatp=ielmat;
  *ielorienp=ielorien;*t0p=t0;*voldp=vold;*veoldp=veold;*xstatep=xstate;
  *prestrp=prestr;*iamt1p=iamt1;
  
  return;
}
