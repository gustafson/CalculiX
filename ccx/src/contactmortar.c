/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2013 Guido Dhondt                     */

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

/**
 * \brief function to change koordinate system needed for mortar contact
 *
 * author: Saskia Sitzmann
 * e-mail: saskia.ssitzmann@fau.de
 * 
 * \f$ A -> \hat{A}, \quad U -> \hat{U}  \f$
 * @todo: have a look at input variables!
 * @param [in] 		ncont
 * @param [in] 		ntie	# ties
 * @param [in]		tieset	
 * @param [in]		nset 
 * @param [in]		set
 * @param [in]		istartset
 * @param [in]		iendset
 * @param [in]		ialset
 * @param [in]		itietri	(1,i)pointer to node where trangulation starts for i (2,i) pointer to end
 * @param [in]		lakon	element label
 * @param [in]		ipkon 	pointer into field kon
 * @param [in]		kon 	field containing connectivity of elements in succesive order
 * @param [in]		koncont (1:3,i) nodes belonging to triangle i (4,i)element face to wich triangle belongs
 * @param [in]		ne
 * @param [in,out]	cg		center of gravity per triangle
 * @param [in]		straight
 * @param [in]		co		coordinates of nodes
 * @param [in]		vold	displacement of nodes
 * @param [in]		ielmat 
 * @param [in]		cs 
 * @param [in]		elcon
 * @param [in]		istep	current step
 * @param [in]		iinc	current increment 
 * @param [in]		iit	current iteration
 * @param [in]		ncmat_	
 * @param [in]		ntmat_
 * @param [in]		ne0
 * @param [out]		vini
 * @param [in]		nmethod
 * @param [in]		neq
 * @param [in]		nzs		size of au
 * @param [in]		nactdof		???
 * @param [in]		itiefac		pointer into field islavsurf: (1,i) beginning slave_i (2,i) end of slave_i
 * @param [in]		islavsurf	islavsurf(1,i) slaveface i islavsurf(2,i) # integration points generated before looking at face i  
 * @param [in]		islavnode	field containing nodes of slave surface
 * @param [in]		imastnode	field containing nodes of master side
 * @param [in] 		nslavnode 	(i) pointer into field islavnode for tie i
 * @param [in]		nmastnode
 * @param [in]		ncone
 * @param [in]		ad
 * @param [out]		aup		pointer to au
 * @param [out]		b
 * @param [out]		irowp		pointer to irow 
 * @param [out]		icol
 * @param [in]		jq
 * @param [in] 		imastop		field containing connectivity of master triangles
 * @param [in] 		iponoels	pointer for node i into field inoels
 * @param [in] 		inoels	wich stores 1D&2D elements belonging to node plus pointer to next entry
 * @param [in] 		nzsc
 * @param [out] 	aucp	pointer to auc
 * @param [out] 	adc	
 * @param [out] 	irowcp 	pointer to irowwc
 * @param [out] 	jqc
 * @param [in,out] 	islavact	active set
 * @param [out] 	gap	transformed gap at node i
 * @param [out]		bdd
 * @param [out] 	auqdtp		pointer to auqdt	
 * @param [out] 	irowqdtp	pointer to irowwqdt
 * @param [in] 		jqqdt
 * @param [in] 		nzsqdt
 * @param [out] 	nzlc
 * @param [out]		slavnor 	normals in nodes of slave face
 * @param [in] 		bhat	
 * @param [out] 	icolc 
 * @param [in] 		aubdp 
 * @param [in] 		irowbdp
 * @param [in] 		jqbd 
 * @param [in] 		mi
 * @param [in] 		ipe
 * @param [in] 		ime
 * @param [in]		tietol
 * @see multimortar
 */
void contactmortar(int *ncont, int *ntie, char *tieset, int *nset, char *set,
        int *istartset, int *iendset, int *ialset, int *itietri,
        char *lakon, int *ipkon, int *kon, int *koncont, int *ne,
        double *cg, double *straight, double *co,
        double *vold, int *ielmat, double *cs, double *elcon,
        int *istep,int *iinc,int *iit,int *ncmat_,int *ntmat_,
        int *ne0, double *vini,
        int *nmethod,int *neq, int *nzs, int *nactdof, int *itiefac,
        int *islavsurf, int *islavnode, int *imastnode,
        int *nslavnode, int *nmastnode, int *ncone, double *ad,
        double **aup, double *b, int **irowp, int *icol, int *jq, int *imastop,
        int *iponoels, int *inoels, int *nzsc, double **aucp,
        double *adc, int **irowcp, int *jqc, int *islavact,
        double *gap, double *bdd, double **auqdtp, int **irowqdtp,
        int *jqqdt, int *nzsqdt, int *nzlc,double *slavnor,double *slavtan, 
        double *bhat,
	int *icolc, double **aubdp, int **irowbdp, int *jqbd, int *mi,
        int *ipe, int *ime,double *tietol,int* iflagact,double *cstress,double *bp_old,int* iflag_fric,
        int *nk, 
        int *nboun,int *ndirboun,int *nodeboun,double *xboun,
        int *nmpc,int *ipompc,int *nodempc,double *coefmpc,
        int *ikboun,int *ilboun,int *ikmpc,int *ilmpc,
        int *nslavspc,int *islavspc,int *nsspc,int *nslavmpc,int *islavmpc,int *nsmpc,
        int *nmastspc,int *imastspc,int *nmspc,int *nmastmpc,int *imastmpc,int *nmmpc,
        int *islavborder,double *pslavdual,
        double **Bdp,double *Dd,int *jqb,int **irowbp, int *nzsbd2, int *islavactdof,
	double *dhinv,int *islavactdoftie,
	double *plicon,int *nplicon, int *npmat_,int *nelcon){
    
    int i,j,k,numb,ntrimax,*nx=NULL,*ny=NULL,*nz=NULL,nintpoint=0,
        nzsbd,*irowbd=NULL,l,nstart,kflag,ntri,ii,number,regmode,
        *irowc=NULL,*imastsurf=NULL,
        *irow=NULL,*irowqdt=NULL,* jqctemp=NULL, 
        debug,*irowb=NULL,nacti, ninacti, nnogap,nstick,nnolm;
    
    double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL,*aubd=NULL, 
	*auc=NULL, *pmastsurf=NULL,*auqdt=NULL,*gapmints=NULL,
	*au=NULL,*pslavsurf=NULL,
	*areaslav=NULL,*anull=NULL,*Bd=NULL,mu,fkinv,p0,beta;
    
    
    clock_t debut;
    clock_t fin;
    irow = *irowp; au=*aup; auc=*aucp; irowc=*irowcp; auqdt=*auqdtp;
    irowqdt=*irowqdtp; aubd=*aubdp; irowbd=*irowbdp;
    Bd=*Bdp; irowb=*irowbp;
    
    debug=0;
    printf("contactmortar: start\n");
    if(*iflagact==0){	
      /* update the location of the center of gravity of 
	   the master triangles and the coefficients of their
	   bounding planes */
      FORTRAN(updatecont,(koncont,ncont,co,vold,
			    cg,straight,mi));		
      /* determining the size of the auxiliary fields 
	   (needed for the master triangle search for any
	   given location on the slave faces */	
      ntrimax=0;	
      for(i=0;i<*ntie;i++){	    
	if(itietri[2*i+1]-itietri[2*i]+1>ntrimax)		
	  ntrimax=itietri[2*i+1]-itietri[2*i]+1;  	
      }
      if ((*iinc==1)&&(*iit==1)){	    
	debut=clock();	    
	xo=NNEW(double,ntrimax);	    
	yo=NNEW(double,ntrimax);	    
	zo=NNEW(double,ntrimax);	    
	x=NNEW(double,ntrimax);	    
	y=NNEW(double,ntrimax);	    
	z=NNEW(double,ntrimax);	   
	nx=NNEW(int,ntrimax);	   
	ny=NNEW(int,ntrimax);	    
	nz=NNEW(int,ntrimax);	    
	areaslav=NNEW(double,itiefac[2*(*ntie-1)+1]);	    
	int ifree=0;	    
	FORTRAN(genfirstactif,(tieset,ntie,itietri,ne,ipkon,kon,lakon,
				   cg,straight,koncont,
				   co,vold,xo,yo,zo,x,y,z,nx,ny,nz,ielmat,cs,elcon,istep,
				   iinc,iit,ncmat_,ntmat_,ne0,vini,nmethod,mi,
				   imastop,nslavnode,islavnode,islavsurf,itiefac,areaslav,iponoels,
				   inoels,set,nset,istartset,iendset,ialset,islavact,&ifree,
				   tietol));
	printf("Frist Active Set : %d nodes\n",ifree);	    
	fin= clock();	    
	printf("genfirstactiv : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);	    
	free(xo);free(yo);free(zo);free(x);free(y);free(z);free(nx);	    
	free(ny);free(nz);	    
	free(areaslav);	
      }
      fflush(stdout);	
      nacti=0; ninacti=0; nnogap=0;nstick=0;nnolm=0;	
      for (i=0;i<*ntie;i++){	    
	printf("tie %d type %c\n",i+1, tieset[i*(81*3)+80]);	    
	if(tieset[i*(81*3)+80]=='C'){	   	      
	  for(j=nslavnode[i];j<nslavnode[i+1];j++){	
	    if(islavact[j]<0){islavact[j]=-3;}
	    if(islavact[j]==2){nacti++;}
            if(islavact[j]==1){nstick++;}
	    if(islavact[j]==0){ninacti++;}
	    if(islavact[j]==-1){nnogap++;}
	    if(islavact[j]==-2){nnolm++;}				      
	  }	    
	}	
      }
      printf("contactmortar: N_Activ: %d\t N_stick: %d\tN_Inactiv: %d\t N_nogap: %d\t N_nolm %d\n",nacti,nstick,ninacti,nnogap,nnolm);		
      xo=NNEW(double,ntrimax);	
      yo=NNEW(double,ntrimax);	
      zo=NNEW(double,ntrimax);	
      x=NNEW(double,ntrimax);	
      y=NNEW(double,ntrimax);	
      z=NNEW(double,ntrimax);	
      nx=NNEW(int,ntrimax);	
      ny=NNEW(int,ntrimax);	
      nz=NNEW(int,ntrimax);	
      /* calculating the normals in the nodes of the slave
	   surface */	
      debut=clock();	
      FORTRAN(gencontrel,(tieset,ntie,itietri,ipkon,kon,
			    lakon,set,cg,straight,
			    koncont,co,vold,nset,
			    iinc,iit,
			    islavsurf,imastsurf,pmastsurf,itiefac,
			    islavnode,nslavnode,slavnor,slavtan,imastop,
			    mi,ncont,ipe,ime,pslavsurf,pslavdual));
	
      fin= clock();
      printf("gencontrel : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);
      /* Calculating the location of the matched slave/master
	   integration points */
		
      debut=clock();
      
//      int imastsuf[THREADS];
//      for ( int i =0; i < THREADS ; ++i )
//	NNEW
      
      imastsurf=NNEW(int,66);	
      gapmints=NNEW(double,66);	
      pmastsurf=NNEW(double,132);	
      pslavsurf=NNEW(double,198);	
      islavsurf[1]=0;	
//#pragma omp for
      for(i=0;i<*ntie;i++){	    
	ii=i+1;	    
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
	  for(l=itiefac[2*i];l<=itiefac[2*i+1];l++){		    
	    RENEW(imastsurf,int,nintpoint+ntri*66);		    
	    RENEW(gapmints,double,nintpoint+ntri*66);		    
	    RENEW(pmastsurf,double,2*(nintpoint+ntri*66));		    
	    RENEW(pslavsurf,double,3*(nintpoint+ntri*66));		    
	    FORTRAN(slavintmortar,(tieset,ntie,itietri,ipkon,kon,
	         	lakon,set,cg,straight,&nintpoint,
	        	koncont,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,nset,
	        	iinc,iit,
	        	islavsurf,imastsurf,pmastsurf,itiefac,
	        	islavnode,nslavnode,slavnor,slavtan,imastop,gapmints,
	        	islavact,mi,ncont,ipe,ime,pslavsurf,pslavdual,&ii,&l,&ntri));
//	    FORTRAN(stop,());
	  }	    
	}	
      }
      fin= clock();	
      printf("slavintmortar : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);	
      printf(" number of slave integration points = %d\n\n",nintpoint);	
      nacti=0; ninacti=0; nnogap=0;nstick=0;nnolm=0;	
      for (i=0;i<*ntie;i++){	
	if(tieset[i*(81*3)+80]=='C'){	    
	  for(j=nslavnode[i];j<nslavnode[i+1];j++){		
	    if(islavact[j]==2){nacti++;}
	    if(islavact[j]==1){nstick++;}		
	    if(islavact[j]==0){ninacti++;}
	    if(islavact[j]==-1){nnogap++;}
	    if(islavact[j]==-2){nnolm++;}	    
	  }	 
	}	 	
      }
      printf("contactmortar: N_Activ: %d\t N_stick: %d\tN_Inactiv: %d\t N_nogap: %d\t N_nolm: %d\n",nacti,nstick,ninacti,nnogap,nnolm);	
      if (nintpoint!=0){	    
	RENEW(imastsurf,int,nintpoint);	
      }else{	    
	RENEW(imastsurf,int,1);	
      }	
      if (nintpoint!=0){	    
	RENEW(gapmints,double,nintpoint);	
      }else{	    
	RENEW(gapmints,double,1);	
      }
      if (nintpoint!=0){    
	RENEW(pmastsurf,double,2*nintpoint);	
      }else{	    
	RENEW(pmastsurf,double,2);	
      }		
      if (nintpoint!=0){	    
	RENEW(pslavsurf,double,3*nintpoint);	
      }else{	    
	RENEW(pslavsurf,double,3);
      }
      free(xo);free(yo);free(zo);free(x);free(y);free(z);free(nx);     
      free(ny);free(nz);	
      /* check SPC's and MPC's on slave nodes for compability */
      debut=clock();
      FORTRAN(checkspcmpc,(lakon,ipkon,kon,ntie,tieset,nset,set,
          itiefac,islavsurf,islavnode,
          imastnode,nslavnode,nmastnode,
          slavnor,slavtan,islavact,
          nboun,ndirboun,nodeboun,xboun,
          nmpc,ipompc,nodempc,coefmpc,
          ikboun,ilboun,ikmpc,ilmpc,
          nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
          nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
          islavborder));
      fin= clock();
      printf("checkspcmpc : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);		
      
      /* calculating the coeffs of dual basis functions */		
      debut=clock();	
      FORTRAN(gendualcoeffs,(tieset,ntie,itietri,ipkon,kon,
			       lakon,set,cg,straight,
			       koncont,co,vold,nset,
			       iinc,iit,islavact,
			       islavsurf,imastsurf,pmastsurf,itiefac,
			       islavnode,nslavnode,imastop,
			       mi,ncont,ipe,ime,pslavsurf,pslavdual));
      fin= clock();	
      printf("gendualcoeffs : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);
	
      /* calculating the coupling matrices bd (described by aubd, irowbd, jqbd and nzsbd; 
	   nonsymmetric, zero diagonal) and dd (in bdd; diagonal matrix */
      nzsbd = nzs[1];	
      if(debug==1){	    
	number=1;	    
	FORTRAN(writematrix,(au,ad,irow,jq,&neq[1],&number));	
      }
      if(debug==1){	    
	number=1;	    
	FORTRAN(writevector,(b,&neq[1],&number));	
      }
      debut=clock();
      bdfill(&irowbd, jqbd, &aubd, bdd, &nzsbd, ntie,
	       ipkon, kon, lakon, nslavnode, nmastnode, imastnode, islavnode, 
	       islavsurf, imastsurf, pmastsurf, itiefac,tieset, neq, nactdof,co,vold,
	       iponoels, inoels,mi,gapmints,gap,pslavsurf,pslavdual,&nintpoint,slavnor,nk,
               nboun,ndirboun,nodeboun,xboun,
               nmpc,ipompc,nodempc,coefmpc,
               ikboun,ilboun,ikmpc,ilmpc,
               nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
               nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
	       &Bd,Dd,jqb,&irowb,nzsbd2,dhinv); 
      fin= clock();
      if(debug==1){	
	number=2;	
	FORTRAN(writematrix,(aubd,bdd,irowbd,jqbd,&neq[1],&number));	
      }
      printf("bdfill : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);	
      free(imastsurf);free(pmastsurf);free(gapmints);free(pslavsurf);
      nacti=0; ninacti=0; nnogap=0;nstick=0;nnolm=0;	
      for (i=0;i<*ntie;i++){	    
	if(tieset[i*(81*3)+80]=='C'){	     
	  for(j=nslavnode[i];j<nslavnode[i+1];j++){		
	    if(islavact[j]==2){nacti++;}
            if(islavact[j]==1){nstick++;}			
	    if(islavact[j]==0){ninacti++;}
	    if(islavact[j]==-1){nnogap++;}
	    if(islavact[j]==-2){nnolm++;}	     
	  }	     	    
	}	
      }
      printf("contactmortar: N_Activ: %d\t N_stick: %d\tN_Inactiv: %d\t N_nogap: %d\t N_nolm: %d \n",nacti,nstick,ninacti,nnogap,nnolm);        
      fflush(stdout);   
    }
    
    /* coupling the active slave degrees of freedom with the corresponding
       slave node */
    FORTRAN(genislavactdof,(ntie,tieset,neq,nactdof,nslavnode,islavact,islavactdof,
			    islavnode,mi));
	
    /* modifying the stiffnes matrix with the coupling matrices; the
	   modified (symmetric) matrix is described in asymmetric form by
	   the fields auc, adc, irowc, jqc and nzsc */         
    nzsbd=jqbd[neq[1]]-1;		
    *nzsc = nzs[1];	
//    RENEW(auc,double, *nzsc);	
//    RENEW(irowc,int, *nzsc);	
    *nzsqdt=nzsbd;	
    auqdt=NNEW(double,*nzsqdt);	
    irowqdt=NNEW(int,*nzsqdt);
    debut=clock();		
    multimortar(au, ad, irow, jq, nzs,
		    aubd, bdd, irowbd, jqbd, &nzsbd,
		    &auc, adc, &irowc, jqc, nzsc,
		    auqdt,irowqdt,jqqdt,nzsqdt,
		    neq,b,bhat,islavnode,imastnode,nactdof,nslavnode,nmastnode,mi,
		    ntie,
		    nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
                    nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,islavact,islavactdof,dhinv,tieset);	
     fin= clock();	
     printf("multimortar : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);	
     number=10;	
     /* sorting the row numbers within each column in irowc */	
     kflag = 2;		
     for (j=0;j<neq[1];j++){	    
       if(jqc[j+1]-jqc[j]>0){		
	 numb=(jqc[j+1]-jqc[j]);		
	 FORTRAN(isortid,(&irowc[jqc[j]-1],&auc[jqc[j]-1],&numb,&kflag));	    
       }	
     }

    /* copying auc,adc,irowc, jqc and bhat into
       au,ad,irow,jq and b */   
    RENEW(au,double,*nzsc);
    RENEW(irow,int,*nzsc);
    for(i=0;i<neq[1];i++){	
      ad[i]=adc[i];	
      jq[i]=jqc[i];	
      b[i]=bhat[i];
    }
    jq[neq[1]]=jqc[neq[1]];
    for(i=0;i<*nzsc;i++){	
      au[i]=auc[i];	
      irow[i]=irowc[i];
    }    
    nzs[1]=jq[neq[1]]-1;
    
    /* changing au due to N and T (normal and tangential
       direction at the slave surface */
    
    
    /** get uhat_k-1 for first increment and first iteration**/
    double *u_old=NULL;
    int mt=mi[1]+1,nodes,id,islavk2,dim,idof1,idof2,jj,kk;
    if((*iinc==1)&&(*iit==1)){     
      u_old=NNEW(double,3*nslavnode[*ntie]);        
      for (i=0;i<*ntie;i++){             
	if(tieset[i*(81*3)+80]=='C'){		
	  for(j=nslavnode[i];j<nslavnode[i+1];j++){	    	  
	    nodes=islavnode[j];	    	  
	    u_old[j*3]=vold[mt*(nodes)-3];	    	  
	    u_old[j*3+1]=vold[mt*(nodes)-2];	    	    	  
	    u_old[j*3+2]=vold[mt*(nodes)-1];	    		
	  }             
	}        
      }
      for(i=0;i<*nk; i++){                  
	for (j=jqb[i]-1;j<jqb[i+1]-1;j++){	
	  nodes=irowb[j];	
	  islavk2=0;	
	  for(jj=0;jj<*ntie;jj++){	 	  
	    if(tieset[jj*(81*3)+80]=='C'){         	    
	      dim=nslavnode[jj+1]-nslavnode[jj];	          	    
	      FORTRAN(nident,(&islavnode[nslavnode[jj]], &nodes,&dim, &id));	 	    
	      if(id>0 && islavnode[nslavnode[jj]+id-1]==nodes){islavk2=nslavnode[jj]+id-1;}	 	  	  
	    }		  	
	  }
	  for(k=0;k<3;k++){	       	  
	    idof1=nactdof[mt*(i+1)-3+k]-1;	       	  
	    idof2=nactdof[mt*nodes-3+k]-1;	       	  
	    jj=-1;	       	  
	    if(idof1>-1){	       	    
	      for(kk=jqqdt[idof1]-1;kk<jqqdt[idof1+1]-1;kk++){			      
		if(irowqdt[kk]-1==idof2){jj=kk;} 	       	   
	      }	       	  
	    }	       	      
	    if(jj>-1) u_old[(islavk2)*3+k]=u_old[(islavk2)*3+k]-auqdt[jj]*vold[mt*(i+1)-3+k];	               	
	  }            		
	}      
      }   
    }
    double u_t1old,u_t2old,nu_told;	
    nacti=0; ninacti=0; nnogap=0;nstick=0;nnolm=0;		
    for (i=0;i<*ntie;i++){	  
      if(tieset[i*(81*3)+80]=='C'){	    
	for(j=nslavnode[i];j<nslavnode[i+1];j++){	      
	  /** adjust active set for first iteration of first increment **/	        
	  if((*iinc==1)&&(*iit==1)){		 
	    u_t1old=u_old[(j)*3+0]*slavtan[(j*6)+0]+u_old[(j)*3+1]*slavtan[(j*6)+1]+u_old[(j)*3+2]*slavtan[(j*6)+2]; 	
	    u_t2old=u_old[(j)*3+0]*slavtan[(j*6)+3]+u_old[(j)*3+1]*slavtan[(j*6)+4]+u_old[(j)*3+2]*slavtan[(j*6)+5];		 
	    nu_told=u_t1old*u_t1old+u_t2old*u_t2old;		
	    //if(islavact[j]>0)printf("node %d nuold %e \n",islavnode[j],nu_told);	    
	    FORTRAN(getcontactparams,(&mu,&regmode,&fkinv,&p0,&beta,tietol,elcon,&i,ncmat_,ntmat_));
	    /// in case of friction and no tangetial displacement node is set stick 
	    if(islavact[j]==2 &&mu>1.e-10 && nu_told<1.e-10){islavact[j]=1;}
	    /// in case of NO friction node must set "slip"
            if(islavact[j]==1 && mu<1.e-10){islavact[j]=2;}
//		if(islavact[j]==2 && friccoeff[j]<1.e-10){islavact[j]=0;}		
	    if(islavact[j]==0 && gap[j]<1.e-9 && mu<1.e-10 ) islavact[j]=2;
	    if(islavact[j]==0 && gap[j]<1.e-9 && mu>1.e-10 && nu_told<1.e-10) islavact[j]=1;
	    if(islavact[j]==0 && gap[j]<1.e-9 && mu>1.e-10 && nu_told>1.e-10 ) islavact[j]=2;
	    if(islavact[j]>0 &&  gap[j]>1.e-9) islavact[j]=0;		
	  }  
	  if(islavact[j]==2){nacti++;}
	  if(islavact[j]==1){nstick++;}		
	  if(islavact[j]==0){ninacti++;}
	  if(islavact[j]==-1){nnogap++;}
	  if(islavact[j]==-2){nnolm++;}		    
	}	  
      }	
    }
    if(*iinc==1 && *iit==1 ){
      /// set initial value for bp_old in first iteration of first increment 
      for (i=0;i<*ntie;i++){	 
	if(tieset[i*(81*3)+80]=='C'){	 
	  for(j=nslavnode[i];j<nslavnode[i+1];j++){	  
	    bp_old[j]=1.0;  	 
	  }	 
	}       
      }
      free(u_old);       
    }
    printf("contactmortar: N_Activ: %d\t N_stick: %d\tN_Inactiv: %d\t N_nogap: %d\t N_nolm: %d \n",nacti,nstick,ninacti,nnogap,nnolm);           
    debut=clock();    
    trafoNTmortar_fric3(neq,nzs,islavactdof,islavact,nslavnode,nmastnode,ncone, ad,au,b,
			irow,
		 jq,nzsc,auc,adc,irowc,jqc,gap,bdd,auqdt,irowqdt,
                 jqqdt,nzsqdt,nzlc,slavnor,slavtan,bhat,aubd,irowbd,jqbd,vold,
		 cstress,bp_old,nactdof,islavnode,ntie,mi,nk,nboun,ndirboun,
		 nodeboun,xboun,nmpc,ipompc,nodempc,coefmpc,ikboun,ilboun,ikmpc,ilmpc,
                 nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,nmastspc,imastspc,nmspc,nmastmpc,
		 imastmpc,nmmpc,Bd,Dd,jqb,irowb,nzsbd2,tieset,islavactdoftie,
		 nelcon,elcon,tietol,ncmat_,ntmat_,plicon,nplicon,npmat_);
    
    fin= clock();
    printf("trafoNTmortar : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);
    fflush(stdout); 
    if(debug==1){	
      anull=NNEW(double,1);	
      number=3;		
      FORTRAN(writematrix,(auc,adc,irowc,jqc,&neq[1],&number));		
      number=4;		
      FORTRAN(writematrix,(au,ad,irow,jq,&neq[1],&number));		
      printf("\n");	
      number=5;		
      FORTRAN(writevector,(b,&neq[1],&number));		
      number=6;		
      free(anull);    
    }
 
    /* calculating icol and icolc (needed for SPOOLES) */   
    for(i=0; i<neq[1]; i++){	
      icol[i] = jq[i+1]-jq[i];
    }   
    free(jqctemp);
    
    /* nzlc is the number of the rightmost column with 
       nonzero off-diagonal terms */    
    number=10;    
    *nzlc=0;
    for(i=neq[1]-1;i>-1;i--){	
      if(icolc[i]>0){	    
	*nzlc=i+1;	    
	break;	
      }    
    }    
    *irowp = irow; *aup=au; *aucp=auc; *irowcp=irowc; *auqdtp=auqdt;
    *irowqdtp=irowqdt; *aubdp=aubd; *irowbdp=irowbd;
    *Bdp=Bd; *irowbp=irowb;
    fflush(stdout);    
    return;
}
