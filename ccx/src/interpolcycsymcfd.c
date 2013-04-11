/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */

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
#include <unistd.h>
#include "CalculiX.h"

void interpolcycsymcfd(int *nkold, double *cotet, int *neold, int *ipkon,
     int *kon, int **nodempcp, int *ipompc, int *nmpc,
     int *ikmpc, int *ilmpc, double **coefmpcp, char *labmpc,
     int *mpcfree, int *memmpc_, char *lakon,int *ncs,int *icscp,
     double *xcs,double *ycs,double *zcs,int *nslav,int *islavcp,
     double *xslav,double *yslav,double *zslav,int *ithermal)
{

    int *kontet=NULL,*ifatet=NULL,*inodfa=NULL,*ipofa=NULL,n1,n2,n3,n4,
	*nnx=NULL,*nny=NULL,*nnz=NULL,*ni=NULL,*kontyp=NULL,
	*iparent=NULL,ifreefa=1,kflag=2,ne,netet,nfield,nselect,
	indexe,nfaces,netet_,nktet=0,j,nodes[4],i,*iselect=NULL,
        *nodempc=NULL,nterms,idof,id,k,index,
        *ielemnr=NULL,imastset=0,*istartset=NULL,*iendset=NULL,*ialset=NULL,
        konl[20],idir,idirmin,idirmax;
    
    int i1[24]={3,7,8,6,4,3,8,1,3,8,5,6,3,5,8,1,2,3,5,6,2,5,3,1};
    int i2[12]={1,2,3,5,1,5,3,4,4,5,3,6};
    
    double *planfa=NULL,*cgtet=NULL,*field=NULL,*value=NULL,ratio[20],
	*x=NULL,*y=NULL,*z=NULL,*xo=NULL,*yo=NULL,*zo=NULL,*coefmpc=NULL,
        xp,yp,zp;

    nodempc=*nodempcp;coefmpc=*coefmpcp;
  
    /* number of tetrahedral nodes = number of nodes in original cfd-mesh */
    
    nktet=*nkold;
    
    /* number of elements */
    
    ne=*neold;
    
    /* storage for the frd-element-type */
    
    kontyp=NNEW(int,ne);
    ielemnr=NNEW(int,ne);
    
    /* generating the tetrahedral elements */
    
    netet=0;
    netet_=22*ne;
    
    iparent=NNEW(int,netet_);
    kontet=NNEW(int,4*netet_);
    ipofa=NNEW(int,4*netet_);
    inodfa=NNEW(int,16*netet_);
    ifatet=NNEW(int,4*netet_);
    planfa=NNEW(double,16*netet_);
    
    /* initialization of fields */
    
    FORTRAN(init,(&nktet,inodfa,ipofa,&netet_));
    
    for(i=0;i<ne;i++){
	indexe=ipkon[i]-1;

        /* check whether element exists */

	if(indexe<0) continue;

	/* check whether fluid element */

	if(strcmp1(&lakon[8*i],"F")!=0) continue;
	ielemnr[i]=i+1;
	if(strcmp1(&lakon[8*i+3],"8")==0){
	    
	    /* C3D8* */
	    
	    kontyp[i]=1;
	    for(j=0;j<6;j++){
		nodes[0]=kon[indexe+i1[4*j]];
		nodes[1]=kon[indexe+i1[4*j+1]];
		nodes[2]=kon[indexe+i1[4*j+2]];
		nodes[3]=kon[indexe+i1[4*j+3]];
		iparent[netet]=i+1;
		netet++;
		FORTRAN(generatetet,(kontet,ifatet,&netet,inodfa,
				     &ifreefa,planfa,ipofa,nodes,cotet));
	    }
	}
	else if(strcmp1(&lakon[8*i+3],"6")==0){
	    
	    /* C3D6 */
	    
	    kontyp[i]=2;
	    for(j=0;j<3;j++){
		nodes[0]=kon[indexe+i2[4*j]];
		nodes[1]=kon[indexe+i2[4*j+1]];
		nodes[2]=kon[indexe+i2[4*j+2]];
		nodes[3]=kon[indexe+i2[4*j+3]];
		iparent[netet]=i+1;
		netet++;
		FORTRAN(generatetet,(kontet,ifatet,&netet,inodfa,
				     &ifreefa,planfa,ipofa,nodes,cotet));
	    }
	}
	else if(strcmp1(&lakon[8*i+3],"4")==0){
	    
	    /* C3D4 */
	    
	    kontyp[i]=3;
	    nodes[0]=kon[indexe+1];
	    nodes[1]=kon[indexe+2];
	    nodes[2]=kon[indexe+3];
	    nodes[3]=kon[indexe+4];
	    iparent[netet]=i+1;
	    netet++;
	    FORTRAN(generatetet,(kontet,ifatet,&netet,inodfa,
				 &ifreefa,planfa,ipofa,nodes,cotet));
	}
    }
    free(ipofa);
    
    nfaces=ifreefa-1;
    
    RENEW(ifatet,int,4*netet);
    RENEW(iparent,int,netet);
    RENEW(planfa,double,4*nfaces);
    
    /* writing the tet mesh in frd format */
    
    FORTRAN(writetetmesh,(kontet,&netet,cotet,&nktet,field,&nfield));
    
    /* calculating the center of gravity of the tetrahedra */
    
    cgtet=NNEW(double,3*netet);
    for(i=0;i<netet;i++){
	n1=kontet[4*i]-1;
	n2=kontet[4*i+1]-1;
	n3=kontet[4*i+2]-1;
	n4=kontet[4*i+3]-1;
	cgtet[3*i]=(cotet[3*n1]+cotet[3*n2]+cotet[3*n3]+cotet[3*n4])/4.;
	cgtet[3*i+1]=(cotet[3*n1+1]+cotet[3*n2+1]+cotet[3*n3+1]+cotet[3*n4+1])/4.;
	cgtet[3*i+2]=(cotet[3*n1+2]+cotet[3*n2+2]+cotet[3*n3+2]+cotet[3*n4+2])/4.;
    }
    
    /* initialization of additional fields */
    
    x=NNEW(double,netet);
    y=NNEW(double,netet);
    z=NNEW(double,netet);
    xo=NNEW(double,netet);
    yo=NNEW(double,netet);
    zo=NNEW(double,netet);
    nnx=NNEW(int,netet);
    nny=NNEW(int,netet);
    nnz=NNEW(int,netet);
    ni=NNEW(int,netet);
    for(i=0;i<netet;i++){
	nnx[i]=i+1;
	nny[i]=i+1;
	nnz[i]=i+1;
	x[i]=cgtet[3*i];
	y[i]=cgtet[3*i+1];
	z[i]=cgtet[3*i+2];
	xo[i]=x[i];
	yo[i]=y[i];
	zo[i]=z[i];
    }
    FORTRAN(dsort,(x,nnx,&netet,&kflag));
    FORTRAN(dsort,(y,nny,&netet,&kflag));
    FORTRAN(dsort,(z,nnz,&netet,&kflag));
    free(cgtet);
    
    nfield=0;nselect=0;
    free(kontet);free(inodfa);

    /* generating MPC's for the master side */

    for(i=0;i<*ncs;i++){
	xp=xcs[i];
	yp=ycs[i];
	zp=zcs[i];
	FORTRAN(basis,(x,y,z,xo,yo,zo,nnx,nny,nnz,
              planfa,ifatet,&nktet,&netet,field,&nfield,
              cotet,kontyp,ipkon,kon,iparent,
              &xp,&yp,&zp,value,ratio,iselect,&nselect,
              istartset,iendset,ialset,&imastset,ielemnr,
	      &nterms,konl));
 		
	if(*ithermal>1){
	    idirmin=0;idirmax=4;
	}else{
	    idirmin=1;idirmax=4;
	}
	for(idir=idirmin;idir<idirmax;idir++){
	    idof=8*(icscp[i]-1)+idir;
	    FORTRAN(nident,(ikmpc,&idof,nmpc,&id));
	    if(id>0){
		if(ikmpc[id-1]==idof)continue;
	    }
	    (*nmpc)++;
	    ipompc[*nmpc-1]=*mpcfree;
//	    strcpy1(&labmpc[20*(*nmpc-1)],"CONTACT             ",20);
	    for(k=*nmpc-1;k>id;k--){
		ikmpc[k]=ikmpc[k-1];
		ilmpc[k]=ilmpc[k-1];
	    }
	    ikmpc[id]=idof;
	    ilmpc[id]=*nmpc;
		    
	    /* first term */
	    
	    nodempc[3**mpcfree-3]=icscp[i];
	    nodempc[3**mpcfree-2]=idir;
	    coefmpc[*mpcfree-1]=1.;
	    index=*mpcfree;
	    *mpcfree=nodempc[3**mpcfree-1];
	    if(*mpcfree==0){
		*mpcfree=*memmpc_+1;
		nodempc[3*index-1]=*mpcfree;
		if(*memmpc_<11)*memmpc_=11;
		*memmpc_=(int)(1.1**memmpc_);
		printf("*INFO in gencontmpc: reallocating nodempc; new size = %d\n\n",*memmpc_);
		RENEW(nodempc,int,3**memmpc_);
		RENEW(coefmpc,double,*memmpc_);
		for(k=*mpcfree;k<*memmpc_;k++){
		    nodempc[3*k-1]=k+1;
		}
		nodempc[3**memmpc_-1]=0;
	    }
		    
	    /* subsequent terms */
	    
	    for(j=0;j<nterms;j++){
		nodempc[3**mpcfree-3]=konl[j];
		nodempc[3**mpcfree-2]=idir;
		coefmpc[*mpcfree-1]=-ratio[j];
		index=*mpcfree;
		*mpcfree=nodempc[3**mpcfree-1];
		if(*mpcfree==0){
		    *mpcfree=*memmpc_+1;
		    nodempc[3*index-1]=*mpcfree;
		    if(*memmpc_<11)*memmpc_=11;
		    *memmpc_=(int)(1.1**memmpc_);
		    printf("*INFO in gencontmpc: reallocating nodempc; new size = %d\n\n",*memmpc_);
		    RENEW(nodempc,int,3**memmpc_);
		    RENEW(coefmpc,double,*memmpc_);
		    for(k=*mpcfree;k<*memmpc_;k++){
			nodempc[3*k-1]=k+1;
		    }
		    nodempc[3**memmpc_-1]=0;
		}
	    }
	}
    }

    /* genercating MPC's for the slave side */

    for(i=0;i<*nslav;i++){
	xp=xslav[i];
	yp=yslav[i];
	zp=zslav[i];
	FORTRAN(basis,(x,y,z,xo,yo,zo,nnx,nny,nnz,
              planfa,ifatet,&nktet,&netet,field,&nfield,
              cotet,kontyp,ipkon,kon,iparent,
              &xp,&yp,&zp,value,ratio,iselect,&nselect,
              istartset,iendset,ialset,&imastset,ielemnr,
	      &nterms,konl));
 		
	if(*ithermal>1){
	    idirmin=0;idirmax=4;
	}else{
	    idirmin=1;idirmax=4;
	}
	for(idir=idirmin;idir<idirmax;idir++){
	    idof=8*(islavcp[i]-1)+idir;
	    FORTRAN(nident,(ikmpc,&idof,nmpc,&id));
	    if(id>0){
		if(ikmpc[id-1]==idof)continue;
	    }
	    (*nmpc)++;
	    ipompc[*nmpc-1]=*mpcfree;
//	    strcpy1(&labmpc[20*(*nmpc-1)],"CONTACT             ",20);
	    for(k=*nmpc-1;k>id;k--){
		ikmpc[k]=ikmpc[k-1];
		ilmpc[k]=ilmpc[k-1];
	    }
	    ikmpc[id]=idof;
	    ilmpc[id]=*nmpc;
		    
	    /* first term */
	    
	    nodempc[3**mpcfree-3]=islavcp[i];
	    nodempc[3**mpcfree-2]=idir;
	    coefmpc[*mpcfree-1]=1.;
	    index=*mpcfree;
	    *mpcfree=nodempc[3**mpcfree-1];
	    if(*mpcfree==0){
		*mpcfree=*memmpc_+1;
		nodempc[3*index-1]=*mpcfree;
		if(*memmpc_<11)*memmpc_=11;
		*memmpc_=(int)(1.1**memmpc_);
		printf("*INFO in gencontmpc: reallocating nodempc; new size = %d\n\n",*memmpc_);
		RENEW(nodempc,int,3**memmpc_);
		RENEW(coefmpc,double,*memmpc_);
		for(k=*mpcfree;k<*memmpc_;k++){
		    nodempc[3*k-1]=k+1;
		}
		nodempc[3**memmpc_-1]=0;
	    }
		    
	    /* subsequent terms */
	    
	    for(j=0;j<nterms;j++){
		nodempc[3**mpcfree-3]=konl[j];
		nodempc[3**mpcfree-2]=idir;
		coefmpc[*mpcfree-1]=-ratio[j];
		index=*mpcfree;
		*mpcfree=nodempc[3**mpcfree-1];
		if(*mpcfree==0){
		    *mpcfree=*memmpc_+1;
		    nodempc[3*index-1]=*mpcfree;
		    if(*memmpc_<11)*memmpc_=11;
		    *memmpc_=(int)(1.1**memmpc_);
		    printf("*INFO in gencontmpc: reallocating nodempc; new size = %d\n\n",*memmpc_);
		    RENEW(nodempc,int,3**memmpc_);
		    RENEW(coefmpc,double,*memmpc_);
		    for(k=*mpcfree;k<*memmpc_;k++){
			nodempc[3*k-1]=k+1;
		    }
		    nodempc[3**memmpc_-1]=0;
		}
	    }
	}
    }
	
    free(nnx);free(nny);free(nnz);free(ifatet);free(kontyp);free(iparent);

    free(x);free(y);free(z);free(xo);free(yo);free(zo);
    free(planfa);free(field);

    *nodempcp=nodempc;*coefmpcp=coefmpc;
    
    return;

}
