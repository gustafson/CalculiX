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
#include <fcntl.h>
#include <ctype.h>

#include "CalculiX.h"
#include "readfrd.h"

void getglobalresults (char *jobnamec,int **integerglobp,double **doubleglobp,
                       int *nboun,int *iamboun,double *xboun, int *nload,
                       char *sideload,int *iamload, int *iglob)
{
 
    char  datin[MAX_LINE_LENGTH],text[13]="            ";
    Summen    anz[1]; 
    Nodes     *node=NULL;
    Elements  *elem=NULL;
    Datasets *lcase=NULL;
    
    int *kontet=NULL,*ifatet=NULL,*inodfa=NULL,*ipofa=NULL,type,n1,n2,n3,n4,
	*nnx=NULL,*nny=NULL,*nnz=NULL,*kon=NULL,*ipkon=NULL,*kontyp=NULL,
	*iparent=NULL,ifreefa=1,kflag=2,ne,netet,numnodes,nkon,
	indexe,istep,loadcase,nfaces,netet_,nktet=0,nfield,j,nodes[4],i,
	read_mode=1,nodenr,*integerglob=NULL,*ielemnr=NULL,istep_global;
    
    int i1[24]={3,7,8,6,4,3,8,1,3,8,5,6,3,5,8,1,2,3,5,6,2,5,3,1};
    int i2[12]={1,2,3,5,1,5,3,4,4,5,3,6};
    int i4[88]={5,20,17,13,20,8,19,16,19,7,18,15,18,6,17,14,
		1,9,12,13,12,11,4,16,11,10,3,15,9,2,10,14,
		9,13,11,12,11,13,16,12,13,17,19,20,13,19,16,20,
		17,14,19,18,19,14,15,18,9,11,14,10,11,15,14,10,
		11,19,16,13,11,15,19,14,11,14,19,17,11,19,13,17,
              11,14,17,9,11,17,13,9};
    int i5[56]={1,7,9,10,7,2,8,11,8,3,9,12,5,13,14,11,
              13,4,15,10,14,15,6,12,11,12,7,10,7,12,9,10,
		11,7,12,8,12,7,9,8,13,15,11,10,11,15,12,10,
		13,11,15,14,15,11,12,14};
    int i6[32]={8,9,10,4,1,5,7,8,7,6,3,10,9,8,10,7,
		8,9,5,7,9,10,6,7,5,6,7,9,5,2,6,9};
    
    double *planfa=NULL,*cotet=NULL,*cgtet=NULL,*field=NULL,
	*x=NULL,*y=NULL,*z=NULL,*xo=NULL,*yo=NULL,*zo=NULL,
	*doubleglob=NULL;
    
    integerglob=*integerglobp;doubleglob=*doubleglobp;
    
    /* reading the global coordinates and the topology from file
       (if any, else return) */
    
    if(strcmp1(&jobnamec[396]," ")==0)return;
    strcpy1(datin,&jobnamec[396],132); 
    for(i=0;i<MAX_LINE_LENGTH;i++){
	if(strcmp1(&datin[i]," ")==0){
	    datin[i]='\0';
	    break;
	}
    }
    
    /* determining the appropriate step number: scanning the SPC
       boundary conditions and distribed facial loads 
       if no global data is needed return*/
    
    istep=0;
    for(i=0;i<*nboun;i++){
	if((xboun[i]<1.9232931375)&&(xboun[i]>1.9232931373)){
	    istep=iamboun[i];
	    break;
	}
    }
    if(istep==0){
	for(i=0;i<*nload;i++){
	    if(strcmp1(&sideload[20*i+2],"SM")==0){
		istep=iamload[2*i];
		break;
	    }
	}
    }
    if(istep==0){
	return;
    }else{
	*iglob=1;
    }
    
    /* initialization of the size of fields used in readfrd.c */

    anz->orign=0;
    anz->n=0;
    anz->e=0;
    anz->f=0;
    anz->g=0;
    anz->t=0;
    anz->l=0;
    anz->olc=0;
    anz->orignmax=0;
    anz->nmax=0;
    anz->nmin=MAX_INTEGER;
    anz->emax=0;
    anz->emin=MAX_INTEGER;
    anz->sets=0;
    anz->mats=0;
    anz->amps=0;
    anz->noffs=0;
    anz->eoffs=0;
    
    readfrd( datin, anz, &node, &elem, &lcase, read_mode);
    
    /* calculation of the highest node number */
    
    nktet=0;
    for(i=0;i<anz[0].n;i++){
	if(node[i].nr>nktet) nktet=node[i].nr;
    }
    cotet=NNEW(double,3*nktet);
    
    /* storing the global coordinates */
    
    for (i=0;i<anz[0].n;i++){
	nodenr=node[i].nr;
	j=nodenr-1;
	cotet[3*j]=node[nodenr].nx;
	cotet[3*j+1]=node[nodenr].ny;
	cotet[3*j+2]=node[nodenr].nz;
    }
    
    /* number of elements (not highest element number; this number
       is not needed) */
    
    ne=anz[0].e;
    
    /* storing the topology */
    
    indexe=0;
    ielemnr=NNEW(int,ne);
    kontyp=NNEW(int,ne);
    ipkon=NNEW(int,ne);
    kon=NNEW(int,20*ne);
    for(i=0;i<anz[0].e;i++){
	ielemnr[i]=elem[i].nr;
	kontyp[i]=elem[i].type;
	ipkon[i]=indexe;
	if(kontyp[i]==1){
	    numnodes=8;
	}else if(kontyp[i]==2){
	    numnodes=6;
	}else if(kontyp[i]==3){
	    numnodes=4;
	}else if(kontyp[i]==4){
	    numnodes=20;
	}else if(kontyp[i]==5){
	    numnodes=15;
	}else if(kontyp[i]==6){
	    numnodes=10;
	}else{
	    printf("*WARNING in getglobalresults.c: element in global\n");
	    printf("         mesh not recognized; cgx element type=%d\n",kontyp[i]);
	    continue;
	}
	for(j=0;j<numnodes;j++){
	    kon[indexe++]=elem[i].nod[j];
	}
    }
    nkon=indexe;
    RENEW(kon,int,nkon);
    
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
	type=kontyp[i];
	indexe=ipkon[i]-1;
	if(type==1){
	    
	    /* C3D8* */
	    
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
	else if(type==2){
	    
	    /* C3D6 */
	    
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
	else if(type==3){
	    
	    /* C3D4 */
	    
	    nodes[0]=kon[indexe+1];
	    nodes[1]=kon[indexe+2];
	    nodes[2]=kon[indexe+3];
	    nodes[3]=kon[indexe+4];
	    iparent[netet]=i+1;
	    netet++;
	    FORTRAN(generatetet,(kontet,ifatet,&netet,inodfa,
				 &ifreefa,planfa,ipofa,nodes,cotet));
	}
	else if(type==4){
	    
	    /* C3D20* */
	    
	    for(j=0;j<22;j++){
		nodes[0]=kon[indexe+i4[4*j]];
		nodes[1]=kon[indexe+i4[4*j+1]];
		nodes[2]=kon[indexe+i4[4*j+2]];
		nodes[3]=kon[indexe+i4[4*j+3]];
		iparent[netet]=i+1;
		netet++;
		FORTRAN(generatetet,(kontet,ifatet,&netet,inodfa,
				     &ifreefa,planfa,ipofa,nodes,cotet));
	    }
	}
	else if(type==5){
	    
	    /* C3D15 */
	    
	    for(j=0;j<14;j++){
		nodes[0]=kon[indexe+i5[4*j]];
		nodes[1]=kon[indexe+i5[4*j+1]];
		nodes[2]=kon[indexe+i5[4*j+2]];
		nodes[3]=kon[indexe+i5[4*j+3]];
		iparent[netet]=i+1;
		netet++;
		FORTRAN(generatetet,(kontet,ifatet,&netet,inodfa,
				     &ifreefa,planfa,ipofa,nodes,cotet));
	    }
	}
	else if(type==6){
	    
	    /* C3D10 */
	    
	    for(j=0;j<8;j++){
		nodes[0]=kon[indexe+i6[4*j]];
		nodes[1]=kon[indexe+i6[4*j+1]];
		nodes[2]=kon[indexe+i6[4*j+2]];
		nodes[3]=kon[indexe+i6[4*j+3]];
		iparent[netet]=i+1;
		netet++;
		FORTRAN(generatetet,(kontet,ifatet,&netet,inodfa,
				     &ifreefa,planfa,ipofa,nodes,cotet));
	    }
	}
    }
    free(ipofa);
    
    nfaces=ifreefa-1;
    
    RENEW(ifatet,int,4*netet);
    RENEW(iparent,int,netet);
    RENEW(planfa,double,4*nfaces);
    
    /* writing the tet mesh in frd format */
    
//    FORTRAN(writetetmesh,(kontet,&netet,cotet,&nktet,field,&nfield));
    
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
    
    /* loading the step data : NDTEMP (1 variable), DISP (3 variables) and
       STRESS (6 variables), if present */
    
    field=NNEW(double,10*nktet);
    
    /* reading the temperatures */
    /* 1. determining the last temperature loadcase in the step */
    
    loadcase=-1;
    for(i=0;i<anz[0].l;i++){
	for(j=0;j<lcase[i].npheader;j++){
	    if(strcmp1(&lcase[i].pheader[j][5],"PSTEP")==0){
		strcpy1(text,&lcase[i].pheader[j][48],12);
		istep_global=atoi(text);
		break;
	    }
	}
	if((istep_global==istep)&&
	   (strcmp1(lcase[i].name,"NDTEMP")==0)){
	    loadcase=i;
	}else if(istep_global>istep){
	    break;
	}
    }
    
    /* 2. reading the data */
    
    if(loadcase>-1){
//	if( readfrdblock(loadcase, anz, node, lcase )==-1) 
	if(!read_mode && readfrdblock(loadcase, anz, node, lcase )==-1) 
	{
	    printf("ERROR in getglobalresults: Could not read data for Dataset:%d\n", i+1); 
	    FORTRAN(stop,());
	}
	
    /* 3. storing the data */
    
	for(i=0;i<anz[0].n;i++){
	    nodenr=node[i].nr;
	    field[10*(nodenr-1)]=lcase[loadcase].dat[0][nodenr];
	}
    }
    
    /* reading the displacements */
    /* 1. determining the last displacement loadcase in the step */
    
    loadcase=-1;
    for(i=0;i<anz[0].l;i++){
	for(j=0;j<lcase[i].npheader;j++){
	    if(strcmp1(&lcase[i].pheader[j][5],"PSTEP")==0){
		strcpy1(text,&lcase[i].pheader[j][48],12);
		istep_global=atoi(text);
		break;
	    }
	}
	if((istep_global==istep)&&
	   (strcmp1(lcase[i].name,"DISPR")==0)){
	    loadcase=i;
	}else if(istep_global>istep){
	    break;
	}
    }
    
    /* 2. reading the data */
    
    if(loadcase>-1){
//	if( readfrdblock(loadcase, anz, node, lcase )==-1) 
	if(!read_mode && readfrdblock(loadcase, anz, node, lcase )==-1) 
	{
	    printf("ERROR in getglobalresults: Could not read data for Dataset:%d\n", i+1); 
	    FORTRAN(stop,());
	}
	
    /* 3. storing the data */
    
	for(i=0;i<anz[0].n;i++){
	    nodenr=node[i].nr;
	    field[10*(nodenr-1)+1]=lcase[loadcase].dat[0][nodenr];
	    field[10*(nodenr-1)+2]=lcase[loadcase].dat[1][nodenr];
	    field[10*(nodenr-1)+3]=lcase[loadcase].dat[2][nodenr];
	}
    }
    
    /* reading the stresses */
    /* 1. determining the last stress loadcase in the step */
    
    loadcase=-1;
    for(i=0;i<anz[0].l;i++){
	for(j=0;j<lcase[i].npheader;j++){
	    if(strcmp1(&lcase[i].pheader[j][5],"PSTEP")==0){
		strcpy1(text,&lcase[i].pheader[j][48],12);
		istep_global=atoi(text);
		break;
	    }
	}
	if((istep_global==istep)&&
	   (strcmp1(lcase[i].name,"STRESS")==0)){
	    loadcase=i;
	}else if(istep_global>istep){
	    break;
	}
    }
    
    /* 2. reading the data */
    
    if(loadcase>-1){
//	if( readfrdblock(loadcase, anz, node, lcase )==-1) 
	if(!read_mode && readfrdblock(loadcase, anz, node, lcase )==-1) 
	{
	    printf("ERROR in getglobalresults: Could not read data for Dataset:%d\n", i+1); 
	    FORTRAN(stop,());
	}
	
    /* 3. storing the data */
    
	for(i=0;i<anz[0].n;i++){
	    nodenr=node[i].nr;
	    field[10*(nodenr-1)+4]=lcase[loadcase].dat[0][nodenr];
	    field[10*(nodenr-1)+5]=lcase[loadcase].dat[1][nodenr];
	    field[10*(nodenr-1)+6]=lcase[loadcase].dat[2][nodenr];
	    field[10*(nodenr-1)+7]=lcase[loadcase].dat[3][nodenr];
	    field[10*(nodenr-1)+8]=lcase[loadcase].dat[4][nodenr];
	    field[10*(nodenr-1)+9]=lcase[loadcase].dat[5][nodenr];
	}
    }
    
    free(kontet);free(inodfa);
    free(node);free(elem);
    for(j=0;j<anz->l;j++){
      freeDatasets(lcase,j);
    }
    free(lcase);lcase=NULL;
    
    /* storing the global data in a common block */
    
    
    integerglob=NNEW(int,5+3*ne+nkon+8*netet);
    
    integerglob[0]=nktet;
    integerglob[1]=netet;
    integerglob[2]=ne;
    integerglob[3]=nkon;
    integerglob[4]=nfaces;
    memcpy(&integerglob[5],&nnx[0],sizeof(int)*netet);
    memcpy(&integerglob[netet+5],&nny[0],sizeof(int)*netet);
    memcpy(&integerglob[2*netet+5],&nnz[0],sizeof(int)*netet);
    memcpy(&integerglob[3*netet+5],&ifatet[0],sizeof(int)*4*netet);
    memcpy(&integerglob[7*netet+5],&kontyp[0],sizeof(int)*ne);
    memcpy(&integerglob[ne+7*netet+5],&ipkon[0],sizeof(int)*ne);
    memcpy(&integerglob[2*ne+7*netet+5],&kon[0],sizeof(int)*nkon);
    memcpy(&integerglob[nkon+2*ne+7*netet+5],&iparent[0],sizeof(int)*netet);
    memcpy(&integerglob[nkon+2*ne+8*netet+5],&ielemnr[0],sizeof(int)*ne);
    
    doubleglob=NNEW(double,13*nktet+4*nfaces+6*netet);
    
    memcpy(&doubleglob[0],&x[0],sizeof(double)*netet);
    memcpy(&doubleglob[netet],&y[0],sizeof(double)*netet);
    memcpy(&doubleglob[2*netet],&z[0],sizeof(double)*netet);
    memcpy(&doubleglob[3*netet],&xo[0],sizeof(double)*netet);
    memcpy(&doubleglob[4*netet],&yo[0],sizeof(double)*netet);
    memcpy(&doubleglob[5*netet],&zo[0],sizeof(double)*netet);
    memcpy(&doubleglob[6*netet],&planfa[0],sizeof(double)*4*nfaces);
    memcpy(&doubleglob[4*nfaces+6*netet],&field[0],sizeof(double)*10*nktet);
    memcpy(&doubleglob[10*nktet+4*nfaces+6*netet],&cotet[0],sizeof(double)*3*nktet);
    
    free(nnx);free(nny);free(nnz);free(ifatet);free(kontyp);free(ipkon);
    free(kon);free(iparent);free(ielemnr);

    free(x);free(y);free(z);free(xo);free(yo);free(zo);
    free(planfa);free(field);free(cotet);

    *integerglobp=integerglob;*doubleglobp=doubleglob;
    
    return;

}
