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

void inicont(int * nk,int *ncont, int *ntie, char *tieset, int *nset, char *set,
               int *istartset, int *iendset, int *ialset, int **itietrip,
               char *lakon, int *ipkon, int *kon, int **koncontp,
               int *ncone, double *tietol, int *ismallsliding, int **itiefacp,
               int **islavsurfp, int **islavnodep, int **imastnodep,
               int **nslavnodep, int **nmastnodep, int *mortar,
               int **imastopp,int *nkon,int **iponoelsp,int **inoelsp){

  char kind[2]="C";

  int *itietri=NULL,*koncont=NULL, *itiefac=NULL, *islavsurf=NULL,
      *islavnode=NULL,*imastnode=NULL,*nslavnode=NULL,*nmastnode=NULL,
      nslavs, nmasts, ifacecount, *ipe=NULL, *ime=NULL, *imastop=NULL,
      *iponoels=NULL,*inoels=NULL,ifreenoels;

  itietri=*itietrip;koncont=*koncontp;itiefac=*itiefacp;islavsurf=*islavsurfp;
  islavnode=*islavnodep;imastnode=*imastnodep;nslavnode=*nslavnodep;
  nmastnode=*nmastnodep;imastop=*imastopp,iponoels=*iponoelsp,
  inoels=*inoelsp;;

  /* determining the number of slave entities (nodes or faces, ncone),
     and the number of master triangles (ncont) */

  FORTRAN(allocont,(ncont,ntie,tieset,nset,set,istartset,iendset,
          ialset,lakon,ncone,tietol,ismallsliding,kind,mortar));
  if(ncont==0) return;

  itietri=NNEW(int,2**ntie);
  koncont=NNEW(int,4**ncont);
  
  /* triangulation of the master side */
  
  FORTRAN(triangucont,(ncont,ntie,tieset,nset,set,istartset,iendset,
          ialset,itietri,lakon,ipkon,kon,koncont,kind));
    
		
  if(*mortar==1){
    
    itiefac=NNEW(int,2**ntie);
    islavsurf=NNEW(int,2**ncone);
    islavnode=NNEW(int,8**ncone);
    imastnode=NNEW(int,3**ncont);
    nslavnode=NNEW(int,*ntie+1);
    nmastnode=NNEW(int,*ntie+1);
    ipe=NNEW(int,*nk);
    ime=NNEW(int,12**ncont);
    imastop=NNEW(int,3**ncont);
    iponoels=NNEW(int,*nk);
    inoels=NNEW(int,3**nkon);
    
    FORTRAN(tiefaccont,(lakon,ipkon,kon,ntie,tieset,nset,set,
       istartset,iendset,ialset,itiefac,islavsurf,islavnode,
       imastnode,nslavnode,nmastnode,&nslavs,&nmasts,&ifacecount,
       ipe,ime,imastop,ncont,koncont,iponoels,inoels,&ifreenoels));

    RENEW(islavsurf, int, 2*ifacecount);
    RENEW(islavnode, int, nslavs);
    RENEW(imastnode, int, nmasts);
    RENEW(inoels,int,3*ifreenoels);
    free(ipe);free(ime);
  }
  
  *itietrip=itietri;*koncontp=koncont;
  *itiefacp=itiefac;*islavsurfp=islavsurf;
  *islavnodep=islavnode;*imastnodep=imastnode;
  *nslavnodep=nslavnode;*nmastnodep=nmastnode;
  *imastopp=imastop;*iponoelsp=iponoels;*inoelsp=inoels;
  
  return;
}
