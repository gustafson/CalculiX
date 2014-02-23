/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2014 Guido Dhondt                     */

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
#include <time.h>
#include "CalculiX.h"

void precontact(int *ncont, int *ntie, char *tieset, int *nset, char *set,
        int *istartset, int *iendset, int *ialset, int *itietri,
        char *lakon, int *ipkon, int *kon, int *koncont, int *ne,
        double *cg, double *straight, double *co,double *vold,
        int *istep,int *iinc,int *iit,int *itiefac,
        int *islavsurf, int *islavnode, int *imastnode,
        int *nslavnode, int *nmastnode,int *imastop,int *mi,
	int *ipe, int *ime,double *tietol,int *iflagact,
	int *nintpoint,double **pslavsurfp,double *xmastnor,double *cs,
	int *mcs,int *ics,double *clearini,int *nslavs){

    /* authors: S. Rakotonanahary, S. Sitzmann and J. Hokkanen */

    int i,j,ntrimax,*nx=NULL,*ny=NULL,*nz=NULL,im,
        l,nstart,kflag,ntri,ii;
    
    double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL,
        *pslavsurf=NULL,*clearslavnode=NULL;
    
    pslavsurf=*pslavsurfp;
    
    /* update the location of the center of gravity of 
       the master triangles and the coefficients of their
       bounding planes */
    
    DMEMSET(xmastnor,0,3*nmastnode[*ntie],0.);
    
    FORTRAN(updatecontpen,(koncont,ncont,co,vold,
			   cg,straight,mi,imastnode,nmastnode,xmastnor,
			   ntie,tieset,nset,set,istartset,
			   iendset,ialset,ipkon,lakon,kon,cs,mcs,ics));
    
/*      FORTRAN(updatecont,(koncont,ncont,co,vold,
	cg,straight,mi));*/	
    
    /* determining the size of the auxiliary fields 
       (needed for the master triangle search for any
	   given location on the slave faces */	
    
    ntrimax=0;	
    for(i=0;i<*ntie;i++){	    
	if(itietri[2*i+1]-itietri[2*i]+1>ntrimax)		
	    ntrimax=itietri[2*i+1]-itietri[2*i]+1;  	
    }
    
    /* only at the start of a new step */
    
//    if ((*iinc==1)&&(*iit<=0)){	    
    if ((*istep==1)&&(*iinc==1)&&(*iit<=0)){	    
	xo=NNEW(double,ntrimax);	    
	yo=NNEW(double,ntrimax);	    
	zo=NNEW(double,ntrimax);	    
	x=NNEW(double,ntrimax);	    
	y=NNEW(double,ntrimax);	    
	z=NNEW(double,ntrimax);	   
	nx=NNEW(int,ntrimax);	   
	ny=NNEW(int,ntrimax);	    
	nz=NNEW(int,ntrimax);

	clearslavnode=NNEW(double,3**nslavs);
	    
	FORTRAN(adjustcontactnodes,(tieset,ntie,itietri,cg,straight,
		co,vold,xo,yo,zo,x,y,z,nx,ny,nz,istep,iinc,iit,
		mi,imastop,nslavnode,islavnode,set,nset,istartset,
		iendset,ialset,tietol,clearini,clearslavnode,itiefac,
                ipkon,kon,lakon,islavsurf));
	    
	free(clearslavnode);
	free(xo);free(yo);free(zo);free(x);free(y);free(z);free(nx);	    
	free(ny);free(nz);	    
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
    
    /* Calculating the location of the matched slave/master
       integration points */
    
    RENEW(pslavsurf,double,198);
    
    /* pointer of islavsurf into field pslavsurf and
       pmastsurf */
    
    islavsurf[1]=0;	
    
    /* loop over all ties */
    
    for(i=0;i<*ntie;i++){
	ii=i+1;	   
	
	/* only active contact ties are treated */
	
	if(tieset[i*(81*3)+80]=='C'){		
	    nstart=itietri[2*i]-1;		
	    ntri=itietri[2*i+1]-nstart;		
	    for(j=0;j<ntri;j++){		    
		xo[j]=cg[(nstart+j)*3];		    
		x[j]=xo[j];		   
		nx[j]=j+1;		    
		yo[j]=cg[(nstart+j)*3+1];		    
		y[j]=yo[j];		    
		ny[j]=j+1;		    
		zo[j]=cg[(nstart+j)*3+2];		    
		z[j]=zo[j];		    
		nz[j]=j+1;		
	    }
	    kflag=2;		
	    FORTRAN(dsort,(x,nx,&ntri,&kflag));		
	    FORTRAN(dsort,(y,ny,&ntri,&kflag));		
	    FORTRAN(dsort,(z,nz,&ntri,&kflag));	
	    
	    /* loop over all slave faces belonging to the tie */
	    
	    for(l=itiefac[2*i];l<=itiefac[2*i+1];l++){
		RENEW(pslavsurf,double,3*(*nintpoint+ntri*66));		    
		FORTRAN(slavintpoints,(ntie,itietri,ipkon,kon,
			lakon,straight,nintpoint,koncont,co,vold,
                        xo,yo,zo,x,y,z,nx,ny,nz,islavsurf,
	        	islavnode,nslavnode,imastop,
	        	mi,ncont,ipe,ime,pslavsurf,&ii,&l,&ntri));
	    }	    
	}	
    }
    free(xo);free(yo);free(zo);free(x);free(y);free(z);free(nx);	    
    free(ny);free(nz);
    
    *pslavsurfp=pslavsurf;
    
    return;
}
