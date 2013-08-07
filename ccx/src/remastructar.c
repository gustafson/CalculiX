/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2013 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                                                                       */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"

void remastructar(int *ipompc, double **coefmpcp, int **nodempcp, int *nmpc,
              int *mpcfree, int *nodeboun, int *ndirboun, int *nboun,
              int *ikmpc, int *ilmpc, int *ikboun, int *ilboun,
              char *labmpc, int *nk,
              int *memmpc_, int *icascade, int *maxlenmpc,
              int *kon, int *ipkon, char *lakon, int *ne, int *nnn,
              int *nactdof, int *icol, int *jq, int **irowp, int *isolver,
              int *neq, int *nzs,int *nmethod, int *ithermal,
	      int *iperturb, int *mass, int *mi, int *ics, double *cs,
	      int *mcs){

    /* reconstructs the nonzero locations in the stiffness and mass
       matrix after a change in MPC's or the generation of contact
       spring elements: version for frequency calculations (called
       by arpack and arpackcs)  */

    int *nodempc=NULL,*npn=NULL,*adj=NULL,*xadj=NULL,*iw=NULL,*mmm=NULL,
	*xnpn=NULL,*mast1=NULL,*ipointer=NULL,mpcend,mpcmult,
        callfrommain,i,*irow=NULL;

    double *coefmpc=NULL;
    
    nodempc=*nodempcp;coefmpc=*coefmpcp;irow=*irowp;

    /* decascading the MPC's */

    printf(" Decascading the MPC's\n\n");
   
    callfrommain=0;
    cascade(ipompc,&coefmpc,&nodempc,nmpc,
	    mpcfree,nodeboun,ndirboun,nboun,ikmpc,
	    ilmpc,ikboun,ilboun,&mpcend,&mpcmult,
	    labmpc,nk,memmpc_,icascade,maxlenmpc,
            &callfrommain,iperturb,ithermal);
    
    for(i=1;i<=*nk;++i) nnn[i-1]=i;

    /* renumbering the nodes */

    /*printf(" Renumbering the nodes to decrease the profile:\n");
    
    npn=NNEW(int,20**ne+mpcend);
    adj=NNEW(int,380**ne+mpcmult);
    xadj=NNEW(int,*nk+1);
    iw=NNEW(int,4**nk+1);
    mmm=NNEW(int,*nk);
    xnpn=NNEW(int,*ne+*nmpc+1);
    
     FORTRAN(renumber,(nk,kon,ipkon,lakon,ne,ipompc,nodempc,nmpc,nnn,
	npn,adj,xadj,iw,mmm,xnpn));
    
	free(npn);free(adj);free(xadj);free(iw);free(mmm);free(xnpn);*/

    /* determining the matrix structure */
    
    printf(" Determining the structure of the matrix:\n");
 
    if(nzs[1]<10) nzs[1]=10;   
    mast1=NNEW(int,nzs[1]);
    RENEW(irow,int,nzs[1]);for(i=0;i<nzs[1];i++) irow[i]=0;
  
    if((*mcs==0)||(cs[1]<0)){

	ipointer=NNEW(int,4**nk);
    
	mastruct(nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,nboun,ipompc,
	     nodempc,nmpc,nactdof,icol,jq,&mast1,&irow,isolver,neq,nnn,
	     ikmpc,ilmpc,ipointer,nzs,nmethod,ithermal,
             ikboun,ilboun,iperturb,mi);

    }else{
      
      ipointer=NNEW(int,8**nk);
      
      mastructcs(nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,nboun,
		 ipompc,nodempc,nmpc,nactdof,icol,jq,&mast1,&irow,isolver,
		 neq,nnn,ikmpc,ilmpc,ipointer,nzs,nmethod,
		 ics,cs,labmpc,mcs,mi);
    }

    free(ipointer);free(mast1);
    RENEW(irow,int,nzs[2]);
    
    *nodempcp=nodempc;*coefmpcp=coefmpc;*irowp=irow;

    return;
}


