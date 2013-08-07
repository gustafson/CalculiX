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

void gencontmpc(int *ne, int *iface, char *lakon, int *ipkon, int *kon,
		int *nmpc, int **ikmpcp, int **ilmpcp, int **ipompcp, 
                int *mpcfree,
                double **fmpcp, char **labmpcp, int **nodempcp, int *memmpc_,
                double **coefmpcp, int *nmpc_, int *ikboun, int *nboun){
    
  /*   creates contact MPC's for the middle nodes of the
       dependent surface*/
    
    char *labmpc=NULL;

    int i,j,k,indexe,
        node,id,idir,idof,node1,node2,index,*ipompc=NULL,*ikmpc=NULL,
	*ilmpc=NULL,*nodempc=NULL,jface,jmin,jmax,ij;
		int kk;
    
    int nonei6[9]={7,13,14,8,14,15,9,15,13};
    
    int nonei8[12]={9,17,18,10,18,19,11,19,20,12,20,17};
    
    int nonei10[18]={5,1,2,6,2,3,7,3,1,8,1,4,9,2,4,10,3,4};
    
    int nonei15[27]={7,1,2,8,2,3,9,3,1,10,4,5,11,5,6,12,6,4,
		     13,1,4,14,2,5,15,3,6};
    
    int nonei20[36]={9,1,2,10,2,3,11,3,4,12,4,1,
		     13,5,6,14,6,7,15,7,8,16,8,5,
		     17,1,5,18,2,6,19,3,7,20,4,8};
    
    int ifaceq[48]={4,3,2,1,11,10,9,12,
		    5,6,7,8,13,14,15,16,
		    1,2,6,5,9,18,13,17,
		    2,3,7,6,10,19,14,18,
		    3,4,8,7,11,20,15,19,
		    4,1,5,8,12,17,16,20};

    int ifacet[24]={1,3,2,7,6,5,
		    1,2,4,5,9,8,
		    2,3,4,6,10,9,
		    1,4,3,8,10,7};

    int ifacew[40]={1,3,2,9,8,7,0,0,
		    4,5,6,10,11,12,0,0,
		    1,2,5,4,7,14,10,13,
		    2,3,6,5,8,15,11,14,
		    3,1,4,6,9,13,12,15};
    
    double *fmpc=NULL, *coefmpc=NULL;

    ipompc=*ipompcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
    fmpc=*fmpcp;nodempc=*nodempcp;coefmpc=*coefmpcp;

    /* determining which nodes are middle nodes */

    i=(int)(*iface/10.);
    jface=*iface-i*10;
    
    /* C-convention */
    
    i--;
    
    indexe=ipkon[i];
    if(indexe<0) return;
    if(strcmp1(&lakon[8*i+3],"2")==0){
	if(strcmp1(&lakon[8*i+6]," ")==0){
	    
	    /* genuine 20-node element */
	    
	    for(ij=4;ij<8;ij++){
		j=ifaceq[(jface-1)*8+ij]-1;
		node=kon[indexe+j];
		
		/* create a MPC between node and the 
		   corresponding end nodes */
		
		node1=kon[indexe+nonei20[(j-8)*3+1]-1];
		node2=kon[indexe+nonei20[(j-8)*3+2]-1];
		
		/* create a MPC between node, node1 and node2 */
		
		for(idir=1;idir<4;idir++){
		    idof=8*(node-1)+idir;
		    FORTRAN(nident,(ikboun,&idof,nboun,&id));
		    if(id>0){
			if(ikboun[id-1]==idof)continue;
		    }
		    FORTRAN(nident,(ikmpc,&idof,nmpc,&id));
		    if(id>0){
			if(ikmpc[id-1]==idof)continue;
		    }
		    (*nmpc)++;
		    if(*nmpc>*nmpc_){
			if(*nmpc_<11)*nmpc_=11;
			*nmpc_=(int)(1.1**nmpc_);
			RENEW(ipompc,int,*nmpc_);
			RENEW(labmpc,char,20**nmpc_+1);
			RENEW(ikmpc,int,*nmpc_);
			RENEW(ilmpc,int,*nmpc_);
			RENEW(fmpc,double,*nmpc_);
		    }
		    ipompc[*nmpc-1]=*mpcfree;
		    strcpy1(&labmpc[20*(*nmpc-1)],"CONTACT             ",20);
		    for(k=*nmpc-1;k>id;k--){
			ikmpc[k]=ikmpc[k-1];
			ilmpc[k]=ilmpc[k-1];
		    }
		    ikmpc[id]=idof;
		    ilmpc[id]=*nmpc;
		    
		    /* first term */
		    
		    nodempc[3**mpcfree-3]=node;
		    nodempc[3**mpcfree-2]=idir;
		    coefmpc[*mpcfree-1]=2.;
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
		    
		    /* second term */
		    
		    nodempc[3**mpcfree-3]=node1;
		    nodempc[3**mpcfree-2]=idir;
		    coefmpc[*mpcfree-1]=-1.;
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
		    
		    /* third term */
		    
		    nodempc[3**mpcfree-3]=node2;
		    nodempc[3**mpcfree-2]=idir;
		    coefmpc[*mpcfree-1]=-1.;
		    index=*mpcfree;
		    *mpcfree=nodempc[3**mpcfree-1];
		    nodempc[3*index-1]=0;
		    if(*mpcfree==0){
			*mpcfree=*memmpc_+1;
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
	    } /* j */
		
	}else if(strcmp1(&lakon[8*i+6],"B")!=0){
	    
	    /* plane strain, plane stress, axisymmetric elements
	       or shell elements */
	    
	    if(jface<3){
		jmax=8;
	    }else{
		jmax=5;
	    }
	    for(ij=4;ij<jmax;ij++){
		j=ifaceq[(jface-1)*8+ij]-1;
		if(jface>2){
		    node=(kon[indexe+j]+kon[indexe+j+4])/2;
		    node1=kon[indexe+nonei8[(j-8)*3+1]-1];
		    node2=kon[indexe+nonei8[(j-8)*3+2]-1];
		}else{
		    node=kon[indexe+j];
		    node1=kon[indexe+nonei20[(j-8)*3+1]-1];
		    node2=kon[indexe+nonei20[(j-8)*3+2]-1];
		}
		
		/* create a MPC between node and the 
		   corresponding end nodes */
		
//		printf("gencontmpc %d %d %d\n",j,jface,ij);
//		kk=nonei8[(j-8)*3+1];
//		printf("gencontmpc %d \n",kk);
//		printf("gencontmpc %d \n",indexe);
//		printf("gencontmpc %d \n",indexe+kk-1);
		
		/* create a MPC between node, node1 and node2 */
		
		for(idir=1;idir<3;idir++){
		    idof=8*(node-1)+idir;
		    FORTRAN(nident,(ikboun,&idof,nboun,&id));
		    if(id>0){
			if(ikboun[id-1]==idof)continue;
		    }
		    FORTRAN(nident,(ikmpc,&idof,nmpc,&id));
		    if(id>0){
			if(ikmpc[id-1]==idof)continue;
		    }
		    (*nmpc)++;
		    if(*nmpc>*nmpc_){
			if(*nmpc_<11)*nmpc_=11;
			*nmpc_=(int)(1.1**nmpc_);
			RENEW(ipompc,int,*nmpc_);
			RENEW(labmpc,char,20**nmpc_+1);
			RENEW(ikmpc,int,*nmpc_);
			RENEW(ilmpc,int,*nmpc_);
			RENEW(fmpc,double,*nmpc_);
		    }
		    ipompc[*nmpc-1]=*mpcfree;
		    strcpy1(&labmpc[20*(*nmpc-1)],"CONTACT            ",20);
		    for(k=*nmpc-1;k>id;k--){
			ikmpc[k]=ikmpc[k-1];
			ilmpc[k]=ilmpc[k-1];
		    }
		    ikmpc[id]=idof;
		    ilmpc[id]=*nmpc;
		    
		    /* first term */
		    
		    nodempc[3**mpcfree-3]=node;
		    nodempc[3**mpcfree-2]=idir;
		    coefmpc[*mpcfree-1]=2.;
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
		    
		    /* second term */
		    
		    nodempc[3**mpcfree-3]=node1;
		    nodempc[3**mpcfree-2]=idir;
		    coefmpc[*mpcfree-1]=-1.;
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
		    
		    /* third term */
		    
		    nodempc[3**mpcfree-3]=node2;
		    nodempc[3**mpcfree-2]=idir;
		    coefmpc[*mpcfree-1]=-1.;
		    index=*mpcfree;
		    *mpcfree=nodempc[3**mpcfree-1];
		    nodempc[3*index-1]=0;
		    if(*mpcfree==0){
			*mpcfree=*memmpc_+1;
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
	
    }else if(strcmp1(&lakon[8*i+3],"15")==0){
	
	if(strcmp1(&lakon[8*i+6]," ")==0){
	    
	    /* genuine 15-node element */
	    
	    if(jface<3){
		jmin=3;jmax=6;
	    }else{
		jmin=4;jmax=8;
	    }
	    for(ij=jmin;ij<jmax;ij++){
		j=ifacew[(jface-1)*8+ij]-1;
		node=kon[indexe+j];
		
		/* create a MPC between node and the 
		   corresponding end nodes */
		
		node1=kon[indexe+nonei15[(j-6)*3+1]-1];
		node2=kon[indexe+nonei15[(j-6)*3+2]-1];
		
		/* create a MPC between node, node1 and node2 */
		
		for(idir=1;idir<4;idir++){
		    idof=8*(node-1)+idir;
		    FORTRAN(nident,(ikboun,&idof,nboun,&id));
		    if(id>0){
			if(ikboun[id-1]==idof)continue;
		    }
		    FORTRAN(nident,(ikmpc,&idof,nmpc,&id));
		    if(id>0){
			if(ikmpc[id-1]==idof)continue;
		    }
		    (*nmpc)++;
		    if(*nmpc>*nmpc_){
			if(*nmpc_<11)*nmpc_=11;
			*nmpc_=(int)(1.1**nmpc_);
			RENEW(ipompc,int,*nmpc_);
			RENEW(labmpc,char,20**nmpc_+1);
			RENEW(ikmpc,int,*nmpc_);
			RENEW(ilmpc,int,*nmpc_);
			RENEW(fmpc,double,*nmpc_);
		    }
		    ipompc[*nmpc-1]=*mpcfree;
		    strcpy1(&labmpc[20*(*nmpc-1)],"CONTACT             ",20);
		    for(k=*nmpc-1;k>id;k--){
			ikmpc[k]=ikmpc[k-1];
			ilmpc[k]=ilmpc[k-1];
		    }
		    ikmpc[id]=idof;
		    ilmpc[id]=*nmpc;
		    
		    /* first term */
		    
		    nodempc[3**mpcfree-3]=node;
		    nodempc[3**mpcfree-2]=idir;
		    coefmpc[*mpcfree-1]=2.;
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
		    
		    /* second term */
		    
		    nodempc[3**mpcfree-3]=node1;
		    nodempc[3**mpcfree-2]=idir;
		    coefmpc[*mpcfree-1]=-1.;
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
		    
		    /* third term */
		    
		    nodempc[3**mpcfree-3]=node2;
		    nodempc[3**mpcfree-2]=idir;
		    coefmpc[*mpcfree-1]=-1.;
		    index=*mpcfree;
		    *mpcfree=nodempc[3**mpcfree-1];
		    nodempc[3*index-1]=0;
		    if(*mpcfree==0){
			*mpcfree=*memmpc_+1;
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
	    } /* j */
	    
	}else if(strcmp1(&lakon[8*i+6],"B")!=0){
	    
	    /* plane strain, plane stress, axisymmetric elements
	       or shell elements */
	    
	    if(jface<3){
		jmin=3;jmax=6;
	    }else{
		jmin=4;jmax=5;
	    }
	    for(ij=jmin;ij<jmax;ij++){
		j=ifacew[(jface-1)*8+ij]-1;
		if(jface>2){
		    node=(kon[indexe+j]+kon[indexe+j+3])/2;
		    node1=kon[indexe+nonei6[(j-6)*3+1]-1];
		    node2=kon[indexe+nonei6[(j-6)*3+2]-1];
		}else{
		    node=kon[indexe+j];
		    node1=kon[indexe+nonei15[(j-6)*3+1]-1];
		    node2=kon[indexe+nonei15[(j-6)*3+2]-1];
		}
		
		
		/* create a MPC between node, node1 and node2 */
		
		for(idir=1;idir<3;idir++){
		    idof=8*(node-1)+idir;
		    FORTRAN(nident,(ikboun,&idof,nboun,&id));
		    if(id>0){
			if(ikboun[id-1]==idof)continue;
		    }
		    FORTRAN(nident,(ikmpc,&idof,nmpc,&id));
		    if(id>0){
			if(ikmpc[id-1]==idof)continue;
		    }
		    (*nmpc)++;
		    if(*nmpc>*nmpc_){
			if(*nmpc_<11)*nmpc_=11;
			*nmpc_=(int)(1.1**nmpc_);
			RENEW(ipompc,int,*nmpc_);
			RENEW(labmpc,char,20**nmpc_+1);
			RENEW(ikmpc,int,*nmpc_);
			RENEW(ilmpc,int,*nmpc_);
			RENEW(fmpc,double,*nmpc_);
		    }
		    ipompc[*nmpc-1]=*mpcfree;
		    strcpy1(&labmpc[20*(*nmpc-1)],"CONTACT             ",20);
		    for(k=*nmpc-1;k>id;k--){
			ikmpc[k]=ikmpc[k-1];
			ilmpc[k]=ilmpc[k-1];
		    }
		    ikmpc[id]=idof;
		    ilmpc[id]=*nmpc;
		    
		    /* first term */
		    
		    nodempc[3**mpcfree-3]=node;
		    nodempc[3**mpcfree-2]=idir;
		    coefmpc[*mpcfree-1]=2.;
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
		    
		    /* second term */
		    
		    nodempc[3**mpcfree-3]=node1;
		    nodempc[3**mpcfree-2]=idir;
		    coefmpc[*mpcfree-1]=-1.;
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
		    
		    /* third term */
		    
		    nodempc[3**mpcfree-3]=node2;
		    nodempc[3**mpcfree-2]=idir;
		    coefmpc[*mpcfree-1]=-1.;
		    index=*mpcfree;
		    *mpcfree=nodempc[3**mpcfree-1];
		    nodempc[3*index-1]=0;
		    if(*mpcfree==0){
			*mpcfree=*memmpc_+1;
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
	
        }else if(strcmp1(&lakon[8*i+3],"10")==0){

		/* genuine 10-node element */

	    for(ij=3;ij<6;ij++){
		j=ifacet[(jface-1)*6+ij]-1;
		node=kon[indexe+j];
			
			/* create a MPC between node and the 
			   corresponding end nodes */
			
			node1=kon[indexe+nonei10[(j-4)*3+1]-1];
			node2=kon[indexe+nonei10[(j-4)*3+2]-1];
			
			/* create a MPC between node, node1 and node2 */

			for(idir=1;idir<4;idir++){
			    idof=8*(node-1)+idir;
			    FORTRAN(nident,(ikboun,&idof,nboun,&id));
			    if(id>0){
			      if(ikboun[id-1]==idof)continue;
			    }
			    FORTRAN(nident,(ikmpc,&idof,nmpc,&id));
			    if(id>0){
				if(ikmpc[id-1]==idof)continue;
			    }
			    (*nmpc)++;
			    if(*nmpc>*nmpc_){
				if(*nmpc_<11)*nmpc_=11;
				*nmpc_=(int)(1.1**nmpc_);
				RENEW(ipompc,int,*nmpc_);
				RENEW(labmpc,char,20**nmpc_+1);
				RENEW(ikmpc,int,*nmpc_);
				RENEW(ilmpc,int,*nmpc_);
				RENEW(fmpc,double,*nmpc_);
			    }
			    ipompc[*nmpc-1]=*mpcfree;
			    strcpy1(&labmpc[20*(*nmpc-1)],"CONTACT             ",20);
			    for(k=*nmpc-1;k>id;k--){
				ikmpc[k]=ikmpc[k-1];
				ilmpc[k]=ilmpc[k-1];
			    }
			    ikmpc[id]=idof;
			    ilmpc[id]=*nmpc;
			    
			    /* first term */
			    
			    nodempc[3**mpcfree-3]=node;
			    nodempc[3**mpcfree-2]=idir;
			    coefmpc[*mpcfree-1]=2.;
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
			    
			    /* second term */
			    
			    nodempc[3**mpcfree-3]=node1;
			    nodempc[3**mpcfree-2]=idir;
			    coefmpc[*mpcfree-1]=-1.;
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
			    
			    /* third term */
			    
			    nodempc[3**mpcfree-3]=node2;
			    nodempc[3**mpcfree-2]=idir;
			    coefmpc[*mpcfree-1]=-1.;
			    index=*mpcfree;
			    *mpcfree=nodempc[3**mpcfree-1];
			    nodempc[3*index-1]=0;
			    if(*mpcfree==0){
				*mpcfree=*memmpc_+1;
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
	    } /* j */

	}

	/*   RENEW(ipompc,int,*nmpc);
    RENEW(labmpc,char,20**nmpc+1);
    RENEW(ikmpc,int,*nmpc);
    RENEW(ilmpc,int,*nmpc);
    RENEW(fmpc,double,*nmpc);*/

    *ipompcp=ipompc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
    *fmpcp=fmpc;*nodempcp=nodempc;*coefmpcp=coefmpc;
    
    return;

}
	  
