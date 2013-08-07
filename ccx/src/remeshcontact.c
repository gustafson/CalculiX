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

void remeshcontact(int *ntie, char *tieset, int *nset, char *set,
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
      nquadface,ninterface,ntotface,*ipkon=NULL,*ielmat=NULL,
      *kon=NULL,*ialset=NULL,nk0,*ielorien=NULL,mt=mi[1]+1,
      *iponoel=NULL,*inoel=NULL,ntets2remesh,*iamt1=NULL;

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

  iponoel=NNEW(int,*nk);
  inoel=NNEW(int,2*4**ne);
  
  FORTRAN(remeshsurf,(tieset,ntie,set,nset,istartset,
		      iendset,ialset,ipkon,kon,lakon,nodface,ipoface,&nface,
		      &nquadface,&ninterface,&ntotface,nk,ne,iponoel,inoel,
                      &ntets2remesh));

  /* if no quadratic faces: return */

  if((ntets2remesh==0)&&(ntotface==0)){
      free(iponoel);free(inoel);free(nodface);free(ipoface);
      *ncont=0;
      return;
  }

  /* if no quadratic tets: free fields iponoel and inoel */

  if(ntets2remesh==0){free(iponoel);free(inoel);}
  
  RENEW(nodface,int,9*ntotface);

  if(ithermal[1]==0){
      ninterface*=3;
  }else if((ithermal[1]==1)||(ithermal[1]>2)){
      ninterface*=4;
  }

  if(nquadface>0){
      
      /* nquadface: upper limit for the C3D20(R) elements
	 in set A (25 extra nodes per element: 1 in the center
         of the element and 3 dummy nodes for each of the 8 subelements
         for the incompatible modes)
	 ntotface: all faces in set A */
      
      RENEW(co,double,3*(*nk+25*nquadface+ntotface));

      RENEW(t0,double,*nk+25*nquadface+ntotface);
      RENEW(iamt1,int,*nk+25*nquadface+ntotface);
      for(k=*nk_;k<*nk+25*nquadface+ntotface;k++){iamt1[k]=0;}

      RENEW(vold,double,mt*(*nk+25*nquadface+ntotface));
      DMEMSET(vold,mt**nk,mt*(*nk+25*nquadface+ntotface),0.);

      RENEW(veold,double,mt*(*nk+25*nquadface+ntotface));
      
      /* reallocating the MPC fields for the new MPC's
	 ninterface: all external faces in set A except
	 those belonging to the contact surface */  
      
      RENEW(ipompc,int,*nmpc_+ninterface);
      RENEW(labmpc,char,20*(*nmpc_+ninterface)+1);
      RENEW(ikmpc,int,*nmpc_+ninterface);
      RENEW(ilmpc,int,*nmpc_+ninterface);
      RENEW(fmpc,double,*nmpc_+ninterface);
      
      /* determining the maximum number of terms;
	 expanding nodempc and coefmpc to accommodate
	 those terms */
      
      neqterms=9*ninterface;
      index=*memmpc_;
      (*memmpc_)+=neqterms;
      RENEW(nodempc,int,3**memmpc_);
      RENEW(coefmpc,double,*memmpc_);
      for(k=index;k<*memmpc_;k++){
	  nodempc[3*k-1]=k+1;
      }
      nodempc[3**memmpc_-1]=0;
      
      (*nmpc_)+=ninterface;
  }
  
  /* reallocating the fields for the topology 
     (upper limits) */
  
  RENEW(ipkon,int,*ne+8*nface+8*ntets2remesh);
  RENEW(ielmat,int,mi[2]*(*ne+8*nface+8*ntets2remesh));
  RENEW(ielorien,int,mi[2]*(*ne+8*nface+8*ntets2remesh));
  if(*iprestr>0){RENEW(prestr,double,6*mi[0]*(*ne+8*nface+8*ntets2remesh));}
  if(*nstate_>0){
      RENEW(xstate,double,*nstate_*mi[0]*(*ne+8*nface+8*ntets2remesh));
      for(k=*nxstate;k<*nstate_*mi[0]*(*ne+8*nface+8*ntets2remesh);k++){
	  xstate[k]=0.;}
      *nxstate=*ne+8*nface+8*ntets2remesh;
  }

  /* 8 new elements with max 11 nodes (C3D8I) 
     8 new elements with 4 nodes (C3D4) */

  RENEW(kon,int,*nkon+88*nface+32*ntets2remesh);
  RENEW(lakon,char,8*(*ne+8*nface+8*ntets2remesh));
  
  /* reallocating the fields for the contact surface */
  
  RENEW(ialset,int,*nalset+4*nface);
  
  /* for C3D20(R) and C3D15 elements:
     generating additional nodes
     setting up the multiple point constraints
     in general:
     generating the new topology
     updating the contact surfaces */
  
  nk0=*nk;

  FORTRAN(remeshcontactel,(tieset,ntie,set,nset,istartset,
	  iendset,ialset,ipkon,kon,nkon,lakon,nodface,ipoface,
	  nk,ipompc,nodempc,ikmpc,ilmpc,nmpc,nmpc_,labmpc,coefmpc,
	  mpcfree,nalset,co,ithermal,&nk0,ne,ielmat,ielorien,mi,
	  t0,vold,veold,iponoel,inoel,xstate,nstate_,prestr,iprestr));

  free(nodface);free(ipoface);

  if(ntets2remesh!=0){free(iponoel);free(inoel);}
  
  if(nquadface>0){
      *nk_+=(*nk-nk0);
  }
  
  *ipompcp=ipompc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*nodempcp=nodempc;*coefmpcp=coefmpc;*ipkonp=ipkon;
  *konp=kon;*lakonp=lakon;*cop=co;*ialsetp=ialset;*ielmatp=ielmat;
  *ielorienp=ielorien;*t0p=t0;*voldp=vold;*veoldp=veold;*xstatep=xstate;
  *prestrp=prestr;*iamt1p=iamt1;
  
  /*  for(i=0;i<*nmpc;i++){
      j=i+1;
      FORTRAN(writempc,(ipompc,nodempc,coefmpc,labmpc,&j));
      }*/
  
  return;
}
