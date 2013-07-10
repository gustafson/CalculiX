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
#include <string.h>
#include "CalculiX.h"

void inicont(int * nk,int *ncont, int *ntie, char *tieset, int *nset, char *set,
               int *istartset, int *iendset, int *ialset, int **itietrip,
               char *lakon, int *ipkon, int *kon, int **koncontp,
               int *nslavs, double *tietol, int *ismallsliding, int **itiefacp,
               int **islavsurfp, int **islavnodep, int **imastnodep,
               int **nslavnodep, int **nmastnodep, int *mortar,
               int **imastopp,int *nkon,int **iponoelsp,int **inoelsp,
	       int **ipep, int **imep, int *ne, int *ifacecount,
               int *nmpc, int *mpcfree, int *memmpc_,
               int **ipompcp, char **labmpcp, int **ikmpcp, int **ilmpcp,
               double **fmpcp, int **nodempcp, double **coefmpcp,
	       int *iperturb, int *ikboun, int *nboun, double *co,
	       int *istep,double **xnoelsp){
    
  char kind1[2]="C",kind2[2]="-",*tchar1=NULL,*tchar3=NULL,*labmpc=NULL;
    
  int *itietri=NULL,*koncont=NULL,*itiefac=NULL, *islavsurf=NULL,im,
      *islavnode=NULL,*imastnode=NULL,*nslavnode=NULL,*nmastnode=NULL,
      nmasts,*ipe=NULL,*ime=NULL,*imastop=NULL,
      *iponoels=NULL,*inoels=NULL,ifreenoels,ifreeme,*ipoface=NULL,
      *nodface=NULL,iface,*ipompc=NULL,*ikmpc=NULL,
      *ilmpc=NULL,*nodempc=NULL,nmpc_,i,j,k,ncone;
    
  double *fmpc=NULL,*coefmpc=NULL,*xnoels=NULL;
    
  itietri=*itietrip;koncont=*koncontp;itiefac=*itiefacp;islavsurf=*islavsurfp;
  islavnode=*islavnodep;imastnode=*imastnodep;nslavnode=*nslavnodep;
  nmastnode=*nmastnodep;imastop=*imastopp,iponoels=*iponoelsp;
  inoels=*inoelsp;ipe=*ipep;ime=*imep;xnoels=*xnoelsp;

  ipompc=*ipompcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  fmpc=*fmpcp;nodempc=*nodempcp;coefmpc=*coefmpcp;
  nmpc_=*nmpc;

  /* determining the number of slave entities (nodes or faces, ncone),
     and the number of master triangles (ncont) */

  FORTRAN(allocont,(ncont,ntie,tieset,nset,set,istartset,iendset,
	  ialset,lakon,&ncone,tietol,ismallsliding,kind1,kind2,mortar,
          istep));
  if(ncont==0) return;

  itietri=NNEW(int,2**ntie);
  koncont=NNEW(int,4**ncont);
  
  /* triangulation of the master side */
  
  FORTRAN(triangucont,(ncont,ntie,tieset,nset,set,istartset,iendset,
	  ialset,itietri,lakon,ipkon,kon,koncont,kind1,kind2,co,nk));
    
  RENEW(ipe,int,*nk);
  RENEW(ime,int,12**ncont);
//  memset(&ipe[0],0,sizeof(int)**nk);
//  memset(&ime[0],0,sizeof(int)*12**ncont);
  DMEMSET(ipe,0,*nk,0.);
  DMEMSET(ime,0,12**ncont,0.);
  imastop=NNEW(int,3**ncont);

  FORTRAN(trianeighbor,(ipe,ime,imastop,ncont,koncont,
		        &ifreeme));

  if(*mortar==0){free(ipe);free(ime);}
  else{RENEW(ime,int,4*ifreeme);}

  /* catalogueing the external faces (only for node-to-face
     contact with a nodal slave surface */

  ipoface=NNEW(int,*nk);
  nodface=NNEW(int,5*6**ne);
  FORTRAN(findsurface,(ipoface,nodface,ne,ipkon,kon,lakon,ntie,
		 tieset));
    
  itiefac=NNEW(int,2**ntie);
  islavsurf=NNEW(int,2*6**ne);
  islavnode=NNEW(int,8*ncone);
  nslavnode=NNEW(int,*ntie+1);
  iponoels=NNEW(int,*nk);
  inoels=NNEW(int,2**nkon);
  xnoels=NNEW(double,*nkon);
  
  imastnode=NNEW(int,3**ncont);
  nmastnode=NNEW(int,*ntie+1);
  
  /* catalogueing the slave faces and slave nodes 
     catalogueing the master nodes (only for Mortar contact) */

  FORTRAN(tiefaccont,(lakon,ipkon,kon,ntie,tieset,nset,set,
       istartset,iendset,ialset,itiefac,islavsurf,islavnode,
       imastnode,nslavnode,nmastnode,nslavs,&nmasts,ifacecount,
       iponoels,inoels,&ifreenoels,mortar,ipoface,nodface,nk,
       xnoels));

  RENEW(islavsurf,int,2**ifacecount+2);
  RENEW(islavnode,int,*nslavs);
  RENEW(inoels,int,2*ifreenoels);
  RENEW(xnoels,double,ifreenoels);
  free(ipoface);free(nodface);
  
  RENEW(imastnode,int,nmasts);
  
  /* constraining the middle nodes for the slave surfaces (not
     for modal dynamics calculations) */

  if(*iperturb>1){
      for(i=0;i<*ifacecount;i++){
	  iface=islavsurf[2*i];
	  gencontmpc(ne,&iface,lakon,ipkon,kon,nmpc,&ikmpc,&ilmpc,&ipompc,
               mpcfree,&fmpc,&labmpc,&nodempc,memmpc_,&coefmpc,&nmpc_,ikboun,
	       nboun);
      }

  /* constraining the middle nodes for the master surfaces (not
     for modal dynamics calculations) */

      tchar1=NNEW(char,81);
      tchar3=NNEW(char,81);
      for(i=0; i<*ntie; i++){
	if(tieset[i*(81*3)+80]=='C'){
	  //a contact constraint was found, so increase nalset
	  memcpy(tchar3,&tieset[i*(81*3)+81+81],81);
	  tchar3[80]='\0';
	  for(j=0; j<*nset; j++){
	    memcpy(tchar1,&set[j*81],81);
	    tchar1[80]='\0';
	    if(strcmp(tchar1,tchar3)==0){
	      //independent element face surface was found
		for(k=istartset[j]-1;k<iendset[j];k++){
		    iface=ialset[k];
		    gencontmpc(ne,&iface,lakon,ipkon,kon,nmpc,&ikmpc,&ilmpc,
                               &ipompc,mpcfree,&fmpc,&labmpc,&nodempc,
                               memmpc_,&coefmpc,&nmpc_,ikboun,nboun);
		}
		    
	    }
	  }
	}
      }
      free(tchar1);
      free(tchar3);

      RENEW(ipompc,int,*nmpc);
      RENEW(labmpc,char,20**nmpc+1);
      RENEW(ikmpc,int,*nmpc);
      RENEW(ilmpc,int,*nmpc);
      RENEW(fmpc,double,*nmpc);
    
      /*   for(i=0;i<*nmpc;i++){
	j=i+1;
	FORTRAN(writempc,(ipompc,nodempc,coefmpc,labmpc,&j));
	}*/

  }

  *itietrip=itietri;*koncontp=koncont;
  *itiefacp=itiefac;*islavsurfp=islavsurf;
  *islavnodep=islavnode;*imastnodep=imastnode;
  *nslavnodep=nslavnode;*nmastnodep=nmastnode;
  *imastopp=imastop;*iponoelsp=iponoels;*inoelsp=inoels;
  *ipep=ipe;*imep=ime;*xnoelsp=xnoels;

  *ipompcp=ipompc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*nodempcp=nodempc;*coefmpcp=coefmpc;
  
  return;
}
