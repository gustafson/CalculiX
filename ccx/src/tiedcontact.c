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

void tiedcontact(int *ntie, char *tieset, int *nset, char *set,
               int *istartset, int *iendset, int *ialset,
               char *lakon, int *ipkon, int *kon,
	       double *tietol,
               int *nmpc, int *mpcfree, int *memmpc_,
               int **ipompcp, char **labmpcp, int **ikmpcp, int **ilmpcp,
               double **fmpcp, int **nodempcp, double **coefmpcp,
	       int *ithermal, double *co, double *vold, int *cfd,
	       int *nmpc_, int *mi, int *nk,int *istep,int *ikboun,
	       int *nboun,char *kind1,char *kind2){

  char *labmpc=NULL;

  int *itietri=NULL,*koncont=NULL,nconf,i,k,*nx=NULL,im,
      *ny=NULL,*nz=NULL,*ifaceslave=NULL,*istartfield=NULL,
      *iendfield=NULL,*ifield=NULL,ntrimax,index,
      ncont,ncone,*ipompc=NULL,*ikmpc=NULL,
      *ilmpc=NULL,*nodempc=NULL,ismallsliding=0,neq,neqterms,
      nmpctied,mortar=0,*ipe=NULL,*ime=NULL,*imastop=NULL,ifreeme;

  double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL,
    *cg=NULL,*straight=NULL,*fmpc=NULL,*coefmpc=NULL;

  ipompc=*ipompcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  fmpc=*fmpcp;nodempc=*nodempcp;coefmpc=*coefmpcp;

  /* identifying the slave surfaces as nodal or facial surfaces */

  ifaceslave=NNEW(int,*ntie);

  FORTRAN(identifytiedface,(tieset,ntie,set,nset,ifaceslave,kind1));

  /* determining the number of triangles of the triangulation
     of the master surface and the number of entities on the
     slave side */

  FORTRAN(allocont,(&ncont,ntie,tieset,nset,set,istartset,iendset,
	  ialset,lakon,&ncone,tietol,&ismallsliding,kind1,
	  kind2,&mortar,istep));

  if(ncont==0) return;

  /* allocation of space for the triangulation; 
     koncont(1..3,i): nodes belonging to triangle i
     koncont(4,i): face label to which the triangle belongs =
     10*element+side number */

  itietri=NNEW(int,2**ntie);
  koncont=NNEW(int,4*ncont);

  /* triangulation of the master surface */

  FORTRAN(triangucont,(&ncont,ntie,tieset,nset,set,istartset,iendset,
	  ialset,itietri,lakon,ipkon,kon,koncont,kind1,kind2,co,nk));
  
  /* catalogueing the neighbors of the master triangles */
  
  RENEW(ipe,int,*nk);
  RENEW(ime,int,12*ncont);
  DMEMSET(ipe,0,*nk,0.);
  DMEMSET(ime,0,12*ncont,0.);
  imastop=NNEW(int,3*ncont);

  FORTRAN(trianeighbor,(ipe,ime,imastop,&ncont,koncont,
		        &ifreeme));

  free(ipe);free(ime);

  /* allocation of space for the center of gravity of the triangles
     and the 4 describing planes */

  cg=NNEW(double,3*ncont);
  straight=NNEW(double,16*ncont);
  
  FORTRAN(updatecont,(koncont,&ncont,co,vold,cg,straight,mi));
  
  /* determining the nodes belonging to the slave face surfaces */

  istartfield=NNEW(int,*ntie);
  iendfield=NNEW(int,*ntie);
  ifield=NNEW(int,8*ncone);

  FORTRAN(nodestiedface,(tieset,ntie,ipkon,kon,lakon,set,istartset,
       iendset,ialset,nset,ifaceslave,istartfield,iendfield,ifield,
       &nconf,&ncone,kind1));

  /* determining the maximum number of equations neq */

  if(*cfd==1){
    if(ithermal[1]<=1){
      neq=4;
    }else{
      neq=5;
    }
  }else{
    if(ithermal[1]<=1){
      neq=3;
    }else if(ithermal[1]==2){
      neq=1;
    }else{
      neq=4;
    }
  }
  neq*=(ncone+nconf);

  /* reallocating the MPC fields for the new MPC's
     ncone: number of MPC'S due to nodal slave surfaces
     nconf: number of MPC's due to facal slave surfaces */  

  RENEW(ipompc,int,*nmpc_+neq);
  RENEW(labmpc,char,20*(*nmpc_+neq)+1);
  RENEW(ikmpc,int,*nmpc_+neq);
  RENEW(ilmpc,int,*nmpc_+neq);
  RENEW(fmpc,double,*nmpc_+neq);

  /* determining the maximum number of terms;
     expanding nodempc and coefmpc to accommodate
     those terms */
  
  neqterms=9*neq;
  index=*memmpc_;
  (*memmpc_)+=neqterms;
  RENEW(nodempc,int,3**memmpc_);
  RENEW(coefmpc,double,*memmpc_);
  for(k=index;k<*memmpc_;k++){
      nodempc[3*k-1]=k+1;
  }
  nodempc[3**memmpc_-1]=0;

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
  
  /* generating the tie MPC's */

  FORTRAN(gentiedmpc,(tieset,ntie,itietri,ipkon,kon,
	  lakon,set,istartset,iendset,ialset,cg,straight,
	  koncont,co,xo,yo,zo,x,y,z,nx,ny,nz,nset,
	  ifaceslave,istartfield,iendfield,ifield,
	  ipompc,nodempc,coefmpc,nmpc,&nmpctied,mpcfree,ikmpc,ilmpc,
	  labmpc,ithermal,tietol,cfd,&ncont,imastop,ikboun,nboun,kind1));

  (*nmpc_)+=nmpctied;
  
  free(xo);free(yo);free(zo);free(x);free(y);free(z);free(nx);
  free(ny);free(nz);free(imastop);

  free(ifaceslave);free(istartfield);free(iendfield);free(ifield);
  free(itietri);free(koncont);free(cg);free(straight);

  /* reallocating the MPC fields */

  /*  RENEW(ipompc,int,nmpc_);
  RENEW(labmpc,char,20*nmpc_+1);
  RENEW(ikmpc,int,nmpc_);
  RENEW(ilmpc,int,nmpc_);
  RENEW(fmpc,double,nmpc_);*/

  *ipompcp=ipompc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*nodempcp=nodempc;*coefmpcp=coefmpc;

  /*  for(i=0;i<*nmpc;i++){
    j=i+1;
    FORTRAN(writempc,(ipompc,nodempc,coefmpc,labmpc,&j));
    }*/

  return;
}
