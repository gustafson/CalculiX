/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                     */

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

/** determining the contact stress vectors and updating the active
 *       and inactive sets and the Langrange Multipliers (LM) 
 *
 * @param [in,out] 	bhat		in:\f$ \hat{f} \f$ before NTtrafo out: differential displacement
 * @param [in] 		adc		\f$ \hat{k} \f$ before NTtrafo, diagonal terms
 * @param [in] 		auc		...nondiagonal terms
 * @param [in] 		jqc
 * @param [in] 		irowc
 * @param [in] 		neq
 * @param [in] 		gap		transformed gap per slave node
 * @param [in] 		bdd		matrix bdd
 * @param [in,out] 	b		in: differenzial displacement out:real displacement
 * @param [in,out] 	islavact	active set
 * @param [in] 		auqdt		\f$ Q_d^T \f$, non-diagonal elements, diagonalelemets are 1
 * @param [in] 		irowqdt
 * @param [in] 		jqqdt
 * @param [in] 		ntie		# contraints
 * @param [in] 		nslavnode	(i) pointer into field islavnode for tie i
 * @param [in] 		islavnode	field containing nodes of slave surfaces
 * @param [in] 		slavnor		field containing normals in slave nodes
 * @param [in] 		jqc
 * @param [in] 		nzlc
 * @param [in] 		nactdof 
 * @param [in,out] 	iflagact	flag indicating whether activ set has changed
 * @param [in,out] 	cstress		contact stress
 * @param [in] 		mi		field containing connectivity of edges
 * @param [out] 	cdisp		(6*j) dispnormal (6*j+3)stressnormal for current slave node j	
 * @param [in,out] 	alm_old		values of lagrange multipliers
 * @param [in] 		iit		current iteration
 *
*/

void contactstress_fric2(double *bhat, double *adc, double *auc,int *jqc, 
    int *irowc, int *neq, double *gap, double *bdd,double *aubd, 
    int *jqbd, int *irowbd, double *b, int *islavact,
    double *auqdt, int *irowqdt, int *jqqdt, int *ntie, int *nslavnode,
    int *islavnode, int *nmastnode, int *imastnode, double *slavnor,double *slavtan,
    int *icolc, int *nzlc, int *nactdof,int* iflagact,double *cstress, int *mi,
    double *cdisp, double *f_cs, double *f_cm, int *iit,int *iwatchactiv,
    double *vold,double* bp,int *nk,double *friccoeff,
    int *nboun,int *ndirboun,int *nodeboun,double *xboun,
    int *nmpc,int *ipompc,int *nodempc,double *coefmpc,
    int *ikboun,int *ilboun,int *ikmpc,int *ilmpc,
    int *nslavspc,int *islavspc,int *nsspc,int *nslavmpc,int *islavmpc,int *nsmpc,
    int *nmastspc,int *imastspc,int *nmspc,int *nmastmpc,int *imastmpc,int *nmmpc,
    double *pslavdual,int *ipkon,int *kon,char *lakon,
    int *islavsurf,int *itiefac,int *iponoels,int *inoels,int *islavactdof,
    double *dhinv, double *Dd, double *Bd,int *jqb,int *irowb,int *nzsbd2){
    
    int i,j,l,jj,k,kk,idof1,idof2,idof3,nodes,mt=mi[1]+1,nstick=0,nslip=0,ninacti=0,nnogap=0,nolm=0,
        max_node,debug,idof,node1,node2,dirind,dirdep,index,ist,idofs,islavk2,id,dim,idofm,nodem,
        keepset;
    
    double aux,stressnormal,ddispnormal,*unitmatrix=NULL,t1,t2,nt1,nt2,
        constant=1.E10,stresst1,stresst2,ddispt1,ddispt2,disp_totalnormal,disp_totalt1,
        disp_totalt2,uholdt1,uholdt2,bpold,nw_t=0.0,du[3],f_csn, f_cmn,dd[3], ch,dhs,coefdep,
        w_t[3],*alm_new=NULL,max=0.0,lm_res,*delta_lm=NULL,*bhat2=NULL,
         *f_cs_tot=NULL, *f_cm_tot=NULL, *vectornull=NULL,
        *cstress2=NULL,*f_res=NULL,*u_old=NULL,dm[3],help[3],ndm;
    
    clock_t debut;
    clock_t fin;
    
    debug=0;
    keepset=0;
//    if(*iflagact==0)keepset=1;
    *iflagact=0;
    cstress2=NNEW(double,neq[1]);
    f_res=NNEW(double,neq[1]);    
    bhat2=NNEW(double,neq[1]);
    debut = clock();
    FORTRAN(opnonsym,(&neq[1],&aux,b,cstress2,adc,auc,jqc,irowc));
    u_old=NNEW(double,3*nslavnode[*ntie]);
    
    /** get uhat_k-1 **/
    
    for (i=0;i<*ntie;i++){
	for(j=nslavnode[i];j<nslavnode[i+1];j++){
	    nodes=islavnode[j];
	    u_old[j*3]=vold[mt*(nodes)-3];
	    u_old[j*3+1]=vold[mt*(nodes)-2];
	    u_old[j*3+2]=vold[mt*(nodes)-1];	
	}
    }

    for(i=0;i<*nk; i++){
      for (j=jqb[i]-1;j<jqb[i+1]-1;j++){
	nodes=irowb[j];
	islavk2=0;
	for(jj=0;jj<*ntie;jj++){         
	  dim=nslavnode[jj+1]-nslavnode[jj];	          
	  FORTRAN(nident,(&islavnode[nslavnode[jj]], &nodes,&dim, &id));	 
	  if(id>0 && islavnode[nslavnode[jj]+id-1]==nodes){islavk2=nslavnode[jj]+id-1;}	
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
//        if(nodes==2654) printf("\tj %d nodem %d coeff %e uold %e %e %e\n",nodes,i+1,1/(Dd[nodes-1])*Bd[j],vold[mt*(i+1)-3+0],vold[mt*(i+1)-3+1], vold[mt*(i+1)-3+2] );	       
//        if(nodes==2654) printf("\tislavk2 %d node %d uold %e %e %e\n\n",islavk2,nodes, u_old[islavk2*3], u_old[islavk2*3+1], u_old[islavk2*3+2] );
      }
    }
    /** loop over slave nodes **/
    for (i=0;i<*ntie;i++){	
      for(j=nslavnode[i];j<nslavnode[i+1];j++){    
	nodes=islavnode[j];	    
	idof1=nactdof[mt*nodes-3]-1;
	idof2=nactdof[mt*nodes-2]-1;
	idof3=nactdof[mt*nodes-1]-1;
        for(k=0;k<3;k++){
	  idof=nactdof[mt*nodes-3+k]-1;
	  if(idof>-1) {	  	  
	    dd[k]=bdd[idof];	  
	    cstress[3*j+k]=0.0;	  
	    for(jj=0;jj<3;jj++){  	    
	      idofs=nactdof[mt*nodes-3+jj]-1;	    
	      if(idofs>-1){		
		cstress[3*j+k]=cstress[3*j+k]+(bhat[idofs]-cstress2[idofs])*dhinv[j*9+3*k+jj];	      	    
	      }	      	     	  
	    }     	  
	    du[k]=b[idof];	   	 
	  }else{	    	  
	    cstress[3*j+k]=0.0;du[k]=0.0;dd[k]=0.0; 	     	  
	    for(jj=nslavmpc[2*(j)];jj<nslavmpc[2*(j)+1];jj++){             	    
	      ist=islavmpc[2*jj];            	    
	      dirdep=nodempc[3*(ist-1)+1];            	    
	      coefdep=coefmpc[ist-1];            	    
	      index=nodempc[3*(ist-1)+2];	     	    
	      node1=nodempc[3*(ist-1)];	     	    
	      if(dirdep==k+1){               	      
		while(index!=0){               		
		  dirind=nodempc[3*(index-1)+1];	       		
		  node2=nodempc[3*(index-1)];	       		
		  ch=0.0;	       		
		  dhs=0.0;	       	       		
		  for(kk=0;kk<3;kk++){	        		  
		    idofs=nactdof[mt*node2-3+kk]-1;               		  
		    if(idofs>-1)islavk2 = floor(islavactdof[idofs]/10.)-1;					       		  
		    if(idofs>-1 && islavk2>-1){	     	         		    
		      ch=ch+(bhat[idofs]-cstress2[idofs])*dhinv[islavk2*9+3*(dirind-1)+kk];	        		  
		    }	       		
		  }	       	       	       
	          idof=nactdof[mt*node2-3+(dirind-1)]-1;	       
		  if(idof>-1){	       		 
		    dhs=b[idof];	       
		  } 	       
	          cstress[3*j+dirdep-1]=cstress[3*j+dirdep-1]-coefmpc[index-1]*ch/coefdep;	       
		  du[dirdep-1]=du[dirdep-1]-coefmpc[index-1]*dhs/coefdep;               
		  index=nodempc[3*(index-1)+2];	       	      
		}	      	    
	      }	     	  
	    }	     	     	  	   	
	  }	
	}

	stressnormal=cstress[3*j+0]*slavnor[3*j]+cstress[3*j+1]*slavnor[3*j+1]+cstress[3*j+2]*slavnor[3*j+2];
	stresst1=cstress[3*j+0]*slavtan[6*j]+cstress[3*j+1]*slavtan[6*j+1]+cstress[3*j+2]*slavtan[6*j+2];
	stresst2=cstress[3*j+0]*slavtan[6*j+3]+cstress[3*j+1]*slavtan[6*j+4]+cstress[3*j+2]*slavtan[6*j+5];
	ddispnormal=du[0]*slavnor[3*j]+du[1]*slavnor[3*j+1]+du[2]*slavnor[3*j+2];
	ddispt1=du[0]*slavtan[6*j]+du[1]*slavtan[6*j+1]+du[2]*slavtan[6*j+2];
	ddispt2=du[0]*slavtan[6*j+3]+du[1]*slavtan[6*j+4]+du[2]*slavtan[6*j+5];		
	disp_totalnormal=(du[0]+u_old[j*3])*slavnor[3*j]
		+(du[1]+u_old[j*3+1])*slavnor[3*j+1]
		+(du[2]+u_old[j*3+2])*slavnor[3*j+2];
	disp_totalt1=(du[0]+u_old[j*3])*slavtan[6*j]
		+(du[1]+u_old[j*3+1])*slavtan[6*j+1]
		+(du[2]+u_old[j*3+2])*slavtan[6*j+2];
	disp_totalt2=(du[0]+u_old[j*3])*slavtan[6*j+3]
		+(du[1]+u_old[j*3+1])*slavtan[6*j+4]
		+(du[2]+u_old[j*3+2])*slavtan[6*j+5];
	uholdt1=u_old[j*3]*slavtan[6*j]
		+u_old[j*3+1]*slavtan[6*j+1]
		+u_old[j*3+2]*slavtan[6*j+2];
	uholdt2=u_old[j*3]*slavtan[6*j+3]
		+u_old[j*3+1]*slavtan[6*j+4]
		+u_old[j*3+2]*slavtan[6*j+5];
	bpold=bp[j];
        if(friccoeff[j]>1.E-10){ 
	  bp[j]=friccoeff[j]*(stressnormal+constant*(ddispnormal-gap[j]));
	}else{
	  bp[j]=(stressnormal+constant*(ddispnormal-gap[j]));
	}
        w_t[0]=stresst1+constant*disp_totalt1;
	w_t[1]=stresst2+constant*disp_totalt2;
	
	nw_t=sqrt( w_t[0]* w_t[0]+w_t[1]*w_t[1]);	   
//	    if(debug==1 && gap[j]<-1.e-10){		
//        if(nodes==1624||nodes==1575||nodes==1526||nodes==1625||nodes==1527||nodes==1576||nodes==1573){
        if(debug==1 && ((islavact[j]==0 &&stressnormal+constant*(ddispnormal-gap[j])>-1E-10) ||(islavact[j]==2 &&stressnormal+constant*(ddispnormal-gap[j])<-1E-10)|| islavact[j]>0)){  
		printf("u(%d): %d %d %d act %d\n", nodes, idof1+1,idof2+1, idof3+1,islavact[j]);
//		printf("\t slavnorn %e %e %e \n",slavnor[3*j],slavnor[3*j+1],slavnor[3*j+2]);
//		printf("\t t1 %e %e %e \n",slavtan[6*j],slavtan[6*j+1],slavtan[6*j+2]);
//		printf("\t t2 %e %e %e \n",slavtan[6*j+3],slavtan[6*j+4],slavtan[6*j+5]);		
		printf("\t bdd(%d): %e %e %e \n",nodes, dd[0], dd[1],dd[2]);
//		printf("\t u_old(%d): %e %e %e \n", nodes, u_old[3*j+0],u_old[3*j+1], u_old[3*j+2]);	
		printf("\t u(%d) : %e %e %e \n", nodes, du[0],du[1], du[2]);
		printf("\t u2(%d): %e %e %e \n", nodes, ddispnormal,ddispt1,ddispt2);
//		printf("\t u_tot2(%d): %e %e %e \n", nodes, disp_totalnormal,disp_totalt1,disp_totalt2);
		printf("\t cstress(%d): %e %e %e, actif: %d  \n",nodes,cstress[3*j+0],cstress[3*j+1],cstress[3*j+2],islavact[j]);
		printf("\t lm(%d)     : %e %e %e \n",nodes, stressnormal, stresst1,stresst2);		
		printf("\t dispnormal(%d)= %e, gap(%d)=%e bp=%e\n",nodes,ddispnormal,nodes,gap[j],bp[j]);
		if(friccoeff[j]>1.E-10){
	        printf("\t uh_t1= %e uh_t2= %e nuh_t= %e \n",disp_totalt1,disp_totalt2,sqrt(disp_totalt1*disp_totalt1+disp_totalt2*disp_totalt2)) ;		
		printf("\t w_t1= %e w_t2= %e nw_t= %e \n",w_t[0],w_t[1],nw_t);
		}
                printf("seta neggap n %d \n",nodes);
		
	}
//	if(islavact[j]>0 && i==3) printf("\t lm(%d): %e %e %e, actif: %d  \n",nodes,stressnormal, stresst1,stresst2,islavact[j]);
	if(keepset==0){
	if(friccoeff[j]>1.E-10){ ///contact tie with friction	
	  if (islavact[j]>-1){	
	    /*** Active/Inactive set condition cf: Hueeber Phd p.72 ***/		    
	    if(bp[j]>-1E-10 && nw_t<bp[j]){			
	      nstick++;		
	      	if(islavact[j]>0 && i==3)printf("\t stick!\n");		
	      if (islavact[j]!=1) {*iflagact = 1;}				
	      islavact[j]=1;			
	      cdisp[6*j]=1;
//			cdisp[6*j+1]=disp_totalt1;			
	      cdisp[6*j+1]=(sqrt(stresst1*stresst1+stresst2*stresst2));			
	      cdisp[6*j+2]=bp[j];			
//			cdisp[6*j+2]=disp_totalt2;
	      cdisp[6*j+3]=stressnormal;			
	      cdisp[6*j+4]=stresst1;			
	      cdisp[6*j+5]=stresst2;	      
	    }else if(bp[j]>-1E-10 && nw_t>bp[j]){			
	      nslip++;			

	      	if(islavact[j]>0 && i==3)printf("\t slip!\n");			
	      if (islavact[j]!=2) {*iflagact = 1;}	
	      islavact[j]=2;			
	      cdisp[6*j]=2;
//			cdisp[6*j+1]=disp_totalt1;
	      cdisp[6*j+1]=(sqrt(stresst1*stresst1+stresst2*stresst2));			
//			cdisp[6*j+2]=disp_totalt2;
	      cdisp[6*j+2]=bp[j];						
	      cdisp[6*j+3]=stressnormal;			
	      cdisp[6*j+4]=stresst1;			
	      cdisp[6*j+5]=stresst2;			    
	    }else{			
	      if (islavact[j]>0){ *iflagact = 1;}	
	      ninacti++;			
	      if(islavact[j]>0 && i==3)printf("\t inactiv!\n");			
	      islavact[j]=0;			
	      cstress[3*j+0]=0.;			
	      cstress[3*j+1]=0.;			
	      cstress[3*j+2]=0.;		 
	      cdisp[6*j]=0.;		        
	      cdisp[6*j+1]=0.;		        
	      cdisp[6*j+2]=0.;		        
	      cdisp[6*j+3]=0.;		        
	      cdisp[6*j+4]=0.;		        
	      cdisp[6*j+5]=0.;		    
	    }		
	  }else{
	    if(islavact[j]==-1){
	      nnogap++;	
	    }else{
	      nolm++;
	    }
//	    if(debug==1)printf("\t nogap!\n");	            
	    cstress[3*j+0]=0.;		    
	    cstress[3*j+1]=0.;		    
	    cstress[3*j+2]=0.;		    
	    cdisp[6*j]=0.;		    
	    cdisp[6*j+1]=0.;		    
	    cdisp[6*j+2]=0.;		    
	    cdisp[6*j+3]=0.;		    
	    cdisp[6*j+4]=0.;		    
	    cdisp[6*j+5]=0.;		
	  }   
	}else{ ///no friction            
	  if (islavact[j]>-1){  		    
	    /*** Active/Inactive set condition cf: Hueeber Phd p.72 ***/		    
	    if(stressnormal+constant*(ddispnormal-gap[j])>-1E-10){
//	      if(islavact[j]>0 && i==3)printf("\t slip!\n");
	      nslip++;			
	      if (islavact[j]!=2) {*iflagact = 1;}
	      islavact[j]=2;
	      cdisp[6*j]=2;
//			cdisp[6*j+1]=disp_totalt1;
	      cdisp[6*j+1]=(sqrt(stresst1*stresst1+stresst2*stresst2));			
//			cdisp[6*j+2]=disp_totalt2;
	      cdisp[6*j+2]=bp[j];						
	      cdisp[6*j+3]=stressnormal;			
	      cdisp[6*j+4]=stresst1;			
	      cdisp[6*j+5]=stresst2;			    
	    }else{			
	      if (islavact[j]!=0){ *iflagact = 1;}
//	      if(islavact[j]>0 && i==3)printf("\t inactiv!\n");
	      ninacti++;			
	      islavact[j]=0;			
	      cstress[3*j+0]=0.;			
	      cstress[3*j+1]=0.;			
	      cstress[3*j+2]=0.;		        
	      cdisp[6*j]=0.;		        
	      cdisp[6*j+1]=0.;
//		        cdisp[6*j+2]=0.;
	      cdisp[6*j+2]=0;						        
	      cdisp[6*j+3]=0.;		        
	      cdisp[6*j+4]=0.;		        
	      cdisp[6*j+5]=0.;		    
	    }		
	  }else{		    
	    if(islavact[j]==-1){
	      nnogap++;	
	    }else{
	      nolm++;
	    }	            
	    cstress[3*j+0]=0.;		    
	    cstress[3*j+1]=0.;		    
	    cstress[3*j+2]=0.;		    
	    cdisp[6*j]=-1.;		    
	    cdisp[6*j+1]=0.;
//		    cdisp[6*j+2]=0.;
            cdisp[6*j+2]=0.0;			    		    
	    cdisp[6*j+3]=0.;		    
	    cdisp[6*j+4]=0.;		    
	    cdisp[6*j+5]=0.;		
	  }			    
	}
	}else{

	 if(islavact[j]>0){			
	      cdisp[6*j]=islavact[j];
//			cdisp[6*j+1]=disp_totalt1;			
	      cdisp[6*j+1]=(sqrt(stresst1*stresst1+stresst2*stresst2));			
	      cdisp[6*j+2]=bp[j];			
//			cdisp[6*j+2]=disp_totalt2;
	      cdisp[6*j+3]=stressnormal;			
	      cdisp[6*j+4]=stresst1;			
	      cdisp[6*j+5]=stresst2;		   
	 }else{
	      cstress[3*j+0]=0.;			
	      cstress[3*j+1]=0.;			
	      cstress[3*j+2]=0.;		        
	      cdisp[6*j]=0.;		        
	      cdisp[6*j+1]=0.;
//		        cdisp[6*j+2]=0.;
	      cdisp[6*j+2]=0;						        
	      cdisp[6*j+3]=0.;		        
	      cdisp[6*j+4]=0.;		        
	      cdisp[6*j+5]=0.;		   
	 }
	  
	}  
	gap[j]=gap[j]-ddispnormal;
//	FORTRAN(stop,());
      }    
    }
    
    if(keepset==1)printf("contactstress_fric2: keep active set!!!\n");
    
    /** TODO:   fill in values for noLM nodes, average over nodes with LM  **/


     unitmatrix=NNEW(double,neq[1]);
     for(j=0;j<neq[1];j++){
	  unitmatrix[j]=1.;
	  bhat2[j]=b[j];
      }  
       
      FORTRAN(opnonsym, (&neq[1], &aux, bhat2, b, unitmatrix, auqdt, jqqdt, irowqdt));
  
/** handle SPC's on master nodes **/
/** @todo: handle SPC's on slave nodes-> get master nodes, but still stored in islavnode **/
/*
  for( i=0; i<*ntie; i++){
        for(j=nmastnode[i]; j<nmastnode[i+1]; j++){
	   nodem = imastnode[j];
	   dm[0]=0.0;dm[1]=0.0;dm[2]=0.0;
//	   printf("nodem %d\n",nodem);
	   for(jj=nmastspc[2*(j)];jj<nmastspc[2*(j)+1];jj++){
	     ist=imastspc[2*jj];
             dirdep=ndirboun[ist-1];
	     node1=nodeboun[ist-1];
	     dm[dirdep-1]=xboun[ist-1];
	     if(debug==1){printf("jj %d ist %d dir %d node %d\n",jj,ist,dirdep,node1);}
	   }
	   for(jj=nmastmpc[2*(j)];jj<nmastmpc[2*(j)+1];jj++){
//	    printf("nodem %d mpc %d type %d\n",nodem,jj,imastmpc[2*jj+1]); 
	    if(imastmpc[2*jj+1]==3){                   
	     ist=imastmpc[2*jj];                                           
	     dirdep=nodempc[3*(ist-1)+1];                                           
	     coefdep=coefmpc[ist-1];                                           
	     index=nodempc[3*(ist-1)+2];
//	     printf("jj %d ist %d dir %d node %d\n",jj,ist,dirdep,nodem);                             
	     while(index!=0){				   				            
		 node1=nodempc[3*(index-1)];                                             
		 dirind=nodempc[3*(index-1)+1];	                                    	    			            
		 idofm = nactdof[mt*(node1-1)+(dirind-1)+1];
		 if(idofm==0){
		   idofm=8*(node1-1)+dirind;
		   FORTRAN(nident,(ikboun,&idofm,nboun,&id));
		   if(id>0 &&ikboun[id-1]==idofm){
		     dm[dirdep-1]=dm[dirdep-1]-coefmpc[index-1]*xboun[ilboun[id-1]]/coefdep;
//		     printf("xboun %e\n",xboun[ilboun[id-1]]);
		   }
		 }
		 index=nodempc[3*(index-1)+2];	                                    
	       }
	       printf("dm %e %e %e\n",dm[0],dm[1],dm[2]);
	     }
	    if(debug==1){printf("jj %d ist %d dir %d node %d\n",jj,ist,dirdep,node1);}
	   }
           ndm=sqrt(dm[0]*dm[0]+dm[1]*dm[1]+dm[2]*dm[2]);
	   if(ndm>1.e-12){
           for (jj=jqb[nodem-1]-1;jj<jqb[nodem]-1;jj++){
	      nodes=irowb[jj];
              for(k=0;k<3;k++){
	       idofs=nactdof[mt*nodes-3+k]-1;
               help[k]=1/(Dd[nodes-1])*Bd[jj]*dm[k];
	       if(idofs>-1)b[idofs]=b[idofs]-help[k];
              }
 
	    }
	   }
	}	
  }      
 */ 
/*       for (i=0;i<*ntie;i++){
	for(j=nslavnode[i];j<nslavnode[i+1];j++){
	    nodes=islavnode[j];
	    idof1=nactdof[mt*nodes-3]-1;
	    idof2=nactdof[mt*nodes-2]-1;
	    idof3=nactdof[mt*nodes-1]-1;
	    printf("u(%d)=%e %e %e \n",nodes,b[idof1],b[idof2],b[idof3]);    
	}
    }
*/    	
      free(unitmatrix);
//      free(auqdt);
    
    	      /* calculating the master contact forces 
		 Ph.D. Thesis Stefan Hartmann eqn. (6.26) */
  
     for (i=0;i<neq[1];i++){
       f_cs[i]=0.0;
       f_cm[i]=0.0;
       cstress2[i]=0.0;       
      }
      for( i=0; i<*ntie; i++){
        for(j=nslavnode[i]; j<nslavnode[i+1]; j++){
	    nodes=islavnode[j];
	    for(k=0;k<3;k++){
	      idof1=nactdof[mt*nodes-3+k]-1;	  
	      if(idof1>-1) {
		f_cs[idof1]=bdd[idof1]*cstress[j*3+k];
		cstress2[idof1]=cstress[j*3+k];
	      }
	    }
	}
      }
/*      for(l=0;l<nslavnode[*ntie];l++){
	nodes=islavnode[l];
	for(j=1;j<mt;j++){
          i=nactdof[mt*(nodes-1)+j];
	  if(i==0){continue;}else{i--;}
	  f_cs[i]=(bhat[i]-cstress2[i]);
	}
      }*/	
       

     
      vectornull=NNEW(double,neq[1]);
      
      FORTRAN(opnonsymt,(&neq[1],&aux,cstress2,f_cm,vectornull,aubd,jqbd,irowbd));
      free(vectornull);
      f_cs_tot=NNEW(double, 3);
      f_cm_tot=NNEW(double, 3);
      for(i=0;i<3;i++){
	f_cs_tot[i]=0.0;
	f_cm_tot[i]=0.0;
      }
      for( i=0; i<*ntie; i++){
        for(j=nslavnode[i]; j<nslavnode[i+1]; j++){
	  nodes = islavnode[j];
	  for(l=0;l<3;l++){
	    idof1=nactdof[mt*nodes-3+l]-1;
//	    if(l==1)b[idof1]=islavact[j];
	    if(idof1>-1)f_cs_tot[l]+=f_cs[idof1];
	    if(idof1>-1)f_cs_tot[l]+=f_cm[idof1];
//	    if((f_cs[idof1]>1.e-17 ||f_cs[idof1]<-1.e-17 ) || (f_cm[idof1]>1.e-17 ||f_cm[idof1]<-1.e-17 )){
//	      if(idof1>-1)printf("nodes %d dof %d \tf_c %e %e\n",nodes,l,f_cs[idof1],f_cm[idof1]);
//	    }
	  }
	}
      }	      
      for( i=0; i<*ntie; i++){
	for(j=nmastnode[i]; j<nmastnode[i+1]; j++){
	  nodes = imastnode[j];
	    for(l=0;l<3;l++){
		idof1=nactdof[mt*nodes-3+l]-1;
		if(idof1>-1){
		 f_cm_tot[l]+=f_cm[idof1];
		}
	    }
	}
 	
      }
 //    for (i=0;i<neq[1];i++){
 //      printf("dof %d \tb %e\n",i+1,b[i]);
 //     }      
      f_csn=sqrt(f_cs_tot[0]*f_cs_tot[0]+f_cs_tot[1]*f_cs_tot[1]+f_cs_tot[2]*f_cs_tot[2]);
      f_cmn=sqrt(f_cm_tot[0]*f_cm_tot[0]+f_cm_tot[1]*f_cm_tot[1]+f_cm_tot[2]*f_cm_tot[2]);
      printf("\nslave contact force : %e %e %e and norm: %e \n",f_cs_tot[0],f_cs_tot[1], f_cs_tot[2], f_csn);
      printf("master contact force: %e %e %e and norm: %e \n",f_cm_tot[0],f_cm_tot[1], f_cm_tot[2], f_cmn );
      free(f_cs_tot);free(f_cm_tot);
 
      if(*iit>10){*iflagact=0;}
      
    printf("\n contacstress : N_stick : %d\t N_slip : %d\tN_Inactiv : %d\t N_nogap : %d  N_nolm : %d\n Flag = %d\n",nstick,nslip,ninacti,nnogap,nolm,*iflagact);
    fin= clock();
    printf("contactstress : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);
//    free(alm_new);
    free(bhat2);free(cstress2);
    return;
}
