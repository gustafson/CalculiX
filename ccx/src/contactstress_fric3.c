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

/** determining the contact stress vectors and updating the active
 *       and inactive sets and the Langrange Multipliers (LM) 
 * 
 * author: Saskia Sitzmann
 * e-mail: saskia.sitzmann@fau.de
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

void contactstress_fric3(double *bhat, double *adc, double *auc,int *jqc, 
    int *irowc, int *neq, double *gap, double *bdd,double *aubd, 
    int *jqbd, int *irowbd, double *b, int *islavact,
    double *auqdt, int *irowqdt, int *jqqdt, int *ntie, int *nslavnode,
    int *islavnode, int *nmastnode, int *imastnode, double *slavnor,double *slavtan,
    int *icolc, int *nzlc, int *nactdof,int* iflagact,double *cstress, int *mi,
    double *cdisp, double *f_cs, double *f_cm, int *iit,int *iwatchactiv,
    double *vold,double* bp,int *nk,
    int *nboun,int *ndirboun,int *nodeboun,double *xboun,
    int *nmpc,int *ipompc,int *nodempc,double *coefmpc,
    int *ikboun,int *ilboun,int *ikmpc,int *ilmpc,
    int *nslavspc,int *islavspc,int *nsspc,int *nslavmpc,int *islavmpc,int *nsmpc,
    int *nmastspc,int *imastspc,int *nmspc,int *nmastmpc,int *imastmpc,int *nmmpc,
    double *pslavdual,int *ipkon,int *kon,char *lakon,
    int *islavsurf,int *itiefac,int *iponoels,int *inoels,int *islavactdof,
    double *dhinv, double *Dd, double *Bd,int *jqb,int *irowb,int *nzsbd2,char *tieset,
    double  *elcon, double *tietol,int *ncmat_,int *ntmat_,
    double *plicon,int *nplicon, int *npmat_,int *nelcon){
    
    int i,j,l,jj,k,kk,idof1,idof2,idof3,nodes,mt=mi[1]+1,nstick=0,nslip=0,ninacti=0,nnogap=0,nolm=0,
        debug,idof,node2,dirind,dirdep,index,ist,idofs,islavk2,id,dim,node1,
        keepset,derivmode,regmode,node_max_ncf,ndof;
    
    double aux,stressnormal,ddispnormal,*unitmatrix=NULL,
        constant=1.E10,stresst1,stresst2,ddispt1,ddispt2,disp_totalnormal,disp_totalt1,
        disp_totalt2,uholdt1,uholdt2,bpold,nw_t=0.0,du[3],f_csn, f_cmn,dd[3], ch,dhs,coefdep,
        w_t[3],*bhat2=NULL,
         *f_cs_tot=NULL, *f_cm_tot=NULL, *vectornull=NULL,
        *cstress2=NULL,*u_old=NULL,aninvloc,gnc,cold[3],
	ln_old,ncf_n, max_ncf_n,mu,p0,beta;
    
    clock_t debut;
    clock_t fin;
    
    debug=0;
    keepset=0;
//    if(*iflagact==0)keepset=1;
    *iflagact=0;
    cstress2=NNEW(double,neq[1]);  
    bhat2=NNEW(double,neq[1]);
    debut = clock();
    FORTRAN(opnonsym,(&neq[1],&aux,b,cstress2,adc,auc,jqc,irowc));
    u_old=NNEW(double,3*nslavnode[*ntie]);
    
    /** get uhat_k-1 **/
    
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
    /** loop over slave nodes **/
    max_ncf_n=0.0;
    for (i=0;i<*ntie;i++){
      if(tieset[i*(81*3)+80]=='C'){      
	for(j=nslavnode[i];j<nslavnode[i+1];j++){ 
	  cold[0]=cstress[3*j];
	  cold[1]=cstress[3*j+1];
	  cold[2]=cstress[3*j+2];
	  FORTRAN(getcontactparams,(&mu,&regmode,&aninvloc,&p0,&beta,tietol,elcon,&i,ncmat_,ntmat_));       
	  nodes=islavnode[j];	    	
	  idof1=nactdof[mt*nodes-3]-1;	
	  idof2=nactdof[mt*nodes-2]-1;	
	  idof3=nactdof[mt*nodes-1]-1; 
	  ndof=0;
	  if(idof1>-1){ndof++;}if(idof2>-1){ndof++;}if(idof3>-1){ndof++;}
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
          ln_old=cold[0]*slavnor[3*j]+cold[1]*slavnor[3*j+1]+cold[2]*slavnor[3*j+2];
	  bpold=bp[j];

	  double gnc_old,dgnc_old;
	  derivmode=0;
	  FORTRAN(regularization_gn_c,(&ln_old,&derivmode,&regmode,&gnc_old,&aninvloc,&p0,&beta,elcon,nelcon,&i,ntmat_,
				       plicon,nplicon,npmat_,ncmat_,tietol));
	  derivmode=1;
	  FORTRAN(regularization_gn_c,(&ln_old,&derivmode,&regmode,&dgnc_old,&aninvloc,&p0,&beta,elcon,nelcon,&i,ntmat_,
				       plicon,nplicon,npmat_,ncmat_,tietol));
	  
	  derivmode=0;
	  FORTRAN(regularization_gn_c,(&stressnormal,&derivmode,&regmode,&gnc,&aninvloc,&p0,&beta,elcon,nelcon,&i,ntmat_,
				       plicon,nplicon,npmat_,ncmat_,tietol));
	  
	  if(mu>1.E-10){ 	  
	    bp[j]=mu*(stressnormal+constant*(ddispnormal-gap[j]-gnc));	
	  }else{	  
	    bp[j]=(stressnormal+constant*(ddispnormal-gap[j]-gnc));	
	  }
          w_t[0]=stresst1+constant*disp_totalt1;	
	  w_t[1]=stresst2+constant*disp_totalt2;		
	  nw_t=sqrt( w_t[0]* w_t[0]+w_t[1]*w_t[1]);  
	  if(bp[j]>0.0){
	    ncf_n=stressnormal-bp[j];
	    if( ncf_n<0.0)ncf_n=-ncf_n;
	  }else{
	    ncf_n=stressnormal;
	    if( ncf_n<0.0)ncf_n=-ncf_n;
	  }
	  if(ncf_n>max_ncf_n && ndof>0){max_ncf_n=ncf_n;node_max_ncf=nodes;}
//	  if(islavact[j]>0||bp[j]>0){
//	    printf("u(%d): %d %d %d act %d\n", nodes, idof1+1,idof2+1, idof3+1,islavact[j]);	    
//	    printf("ln_old=%e gnc %e dgnc %e ln_new=%e, %e=%e ?\n",ln_old,gnc_old,dgnc_old,stressnormal,ddispnormal-dgnc_old*stressnormal,gap[j]+gnc_old-dgnc_old*ln_old);	    
//	    printf("\t dispnormal(%d)= %e, d-gap(%d)=%e bp=%e  gnc=%e\n",nodes,ddispnormal,nodes,ddispnormal-gap[j],bp[j],gnc);
//	    printf("\t ncf_n %e\n",ncf_n);
//	  }	  
	  if(debug==1 && islavact[j]>0){
//          if(nodes==70108||nodes==226941||nodes==70098||nodes==69858||nodes==60358||nodes==60128){
//          if(nodes==43729||nodes==41778||nodes==13282||nodes==65222){
//        if(debug==1 && ((islavact[j]==0 &&stressnormal+constant*(ddispnormal-gap[j])>-1E-10) ||(islavact[j]>0 &&stressnormal+constant*(ddispnormal-gap[j])<-1E-10)|| islavact[j]>0)){  
		printf("u(%d): %d %d %d act %d\n", nodes, idof1+1,idof2+1, idof3+1,islavact[j]);
		printf("\t slavnorn %e %e %e \n",slavnor[3*j],slavnor[3*j+1],slavnor[3*j+2]);
//		printf("\t t1 %e %e %e \n",slavtan[6*j],slavtan[6*j+1],slavtan[6*j+2]);
//		printf("\t t2 %e %e %e \n",slavtan[6*j+3],slavtan[6*j+4],slavtan[6*j+5]);		
//		printf("\t bdd(%d): %e %e %e \n",nodes, dd[0], dd[1],dd[2]);
//		printf("\t u_old(%d): %e %e %e \n", nodes, u_old[3*j+0],u_old[3*j+1], u_old[3*j+2]);	
		printf("\t u(%d) : %e %e %e \n", nodes, du[0],du[1], du[2]);
		printf("\t u2(%d): %e %e %e \n", nodes, ddispnormal,ddispt1,ddispt2);
//		printf("\t u_tot2(%d): %e %e %e \n", nodes, disp_totalnormal,disp_totalt1,disp_totalt2);
		printf("\t cstress(%d): %e %e %e, actif: %d  \n",nodes,cstress[3*j+0],cstress[3*j+1],cstress[3*j+2],islavact[j]);
		printf("\t lm(%d)     : %e %e %e \n",nodes, stressnormal, stresst1,stresst2);		
		printf("\t dispnormal(%d)= %e, gap(%d)=%e bp=%e\n aninv=%e",nodes,ddispnormal,nodes,gap[j],bp[j],gap[j]+gnc);
		if(mu>1.E-10){
	        printf("\t uh_t1= %e uh_t2= %e nuh_t= %e \n",disp_totalt1,disp_totalt2,sqrt(disp_totalt1*disp_totalt1+disp_totalt2*disp_totalt2)) ;		
		printf("\t w_t1= %e w_t2= %e nw_t= %e \n",w_t[0],w_t[1],nw_t);
		}			
	  }
	  
//	  if(*iit>8){keepset=1;*iflagact=0;} 
	  
	  if(keepset==0){
	    if(mu>1.E-10){ ///contact tie with friction	
//	      if(perturbed==1){
//		printf("contactstress: perturbed Lagrange with friction does not work yet!\n stop! \n");
//		FORTRAN(stop,());
//	      }
	      if (islavact[j]>-1){		     
		/*** Active/Inactive set condition cf: Hueeber Phd p.72 ***/		    	    
		if(bp[j]>-1E-10 && nw_t<bp[j]){				      
		  nstick++;			
	          if (islavact[j]!=1) {*iflagact = 1;}				
	          islavact[j]=1;
	          cdisp[6*j]=-(ddispnormal-gap[j]);		  
	          cdisp[6*j+1]=disp_totalt1;
	          cdisp[6*j+2]=disp_totalt2;							      
		  cdisp[6*j+3]=stressnormal;				      
		  cdisp[6*j+4]=stresst1;				      
		  cdisp[6*j+5]=stresst2;	      	      	    
		}else if(bp[j]>-1E-10 && nw_t>bp[j]){				      
		  nslip++;						      
		  if (islavact[j]!=2) {*iflagact = 1;}	
	          islavact[j]=2;
	          cdisp[6*j]=ddispnormal-gap[j];		  
	          cdisp[6*j+1]=disp_totalt1;
	          cdisp[6*j+2]=disp_totalt2;						
	          cdisp[6*j+3]=stressnormal;				      
		  cdisp[6*j+4]=stresst1;				      
		  cdisp[6*j+5]=stresst2;;							      			    	    
		}else{				      
		  if (islavact[j]>0){ *iflagact = 1;}	
	          ninacti++;				      		  			
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
		}else if(islavact[j]==-2){	      
		  nolm++;	    
		}	            
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
		if(stressnormal+constant*(ddispnormal-gnc-gap[j])>-1E-10){	      
		  nslip++;				      
		  if (islavact[j]!=2) {*iflagact = 1;}
	          islavact[j]=2;
	          cdisp[6*j]=-(ddispnormal-gap[j]);
	          cdisp[6*j+1]=disp_totalt1;
	          cdisp[6*j+2]=disp_totalt2;						
	          cdisp[6*j+3]=stressnormal;			
	          cdisp[6*j+4]=stresst1;			
	          cdisp[6*j+5]=stresst2;			    	    
		}else{				          
		  if (islavact[j]!=0){ *iflagact = 1;}
	          ninacti++;				      
		  islavact[j]=0;
		  cstress[3*j+0]=0.;		  				      
		  cstress[3*j+1]=0.;				      
		  cstress[3*j+2]=0.;		        	      
		  cdisp[6*j]=0.;		        	      
		  cdisp[6*j+1]=0.;	      
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
		cdisp[6*j]=0.;			    	    
		cdisp[6*j+1]=0.;            
		cdisp[6*j+2]=0.0;			    		    	    
		cdisp[6*j+3]=0.;		    	    
		cdisp[6*j+4]=0.;		    	    
		cdisp[6*j+5]=0.;			  
	      }			    	
	    }	
	  }else{	 	    
	    if(islavact[j]>0){				     
	      cdisp[6*j]=-(ddispnormal-gap[j]);
	      cdisp[6*j+1]=disp_totalt1;
	      cdisp[6*j+2]=disp_totalt2;	      			
	      cdisp[6*j+3]=stressnormal;			
	      cdisp[6*j+4]=stresst1;			
	      cdisp[6*j+5]=stresst2;		   
	 
	    }else{
	      cstress[3*j+0]=0.;			
	      cstress[3*j+1]=0.;			
	      cstress[3*j+2]=0.;		        
	      cdisp[6*j]=0.;		        
	      cdisp[6*j+1]=0.;
	      cdisp[6*j+2]=0;						        
	      cdisp[6*j+3]=0.;		        
	      cdisp[6*j+4]=0.;		        
	      cdisp[6*j+5]=0.;		   	 
	    }	  	
	  }  	
	gap[j]=gap[j]-ddispnormal;      
	}       
      }   
    }    
    if(keepset==1)printf("contactstress_fric2: keep active set!!!\n");     
    unitmatrix=NNEW(double,neq[1]);
  
    printf("\n max_ncf_n %e node %d\n",max_ncf_n,node_max_ncf);
    if(max_ncf_n>1.e-4 && regmode>1 ){*iflagact=1;}
    
    for(j=0;j<neq[1];j++){	  
      unitmatrix[j]=1.;	  
      bhat2[j]=b[j];      
    }  
    FORTRAN(opnonsym, (&neq[1], &aux, bhat2, b, unitmatrix, auqdt, jqqdt, irowqdt));
 
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
    /* calculating the master contact forces 
		 Ph.D. Thesis Stefan Hartmann eqn. (6.26) */
    for (i=0;i<neq[1];i++){       
      f_cs[i]=0.0;       
      f_cm[i]=0.0;       
      cstress2[i]=0.0;             
    }
    for( i=0; i<*ntie; i++){       
      if(tieset[i*(81*3)+80]=='C'){        
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
    }
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
      if(tieset[i*(81*3)+80]=='C'){	        
	for(j=nslavnode[i]; j<nslavnode[i+1]; j++){	  
	  nodes = islavnode[j];	  
	  for(l=0;l<3;l++){	    
	    idof1=nactdof[mt*nodes-3+l]-1;
	    if(idof1>-1)f_cs_tot[l]+=f_cs[idof1];
	    if(idof1>-1)f_cs_tot[l]+=f_cm[idof1];
//	    if(((f_cs[idof1]>1.e-17 ||f_cs[idof1]<-1.e-17 ) || (f_cm[idof1]>1.e-17 ||f_cm[idof1]<-1.e-17 ))&&l==0){
//	      if(idof1>-1)printf("nodes %d dof %d \tf_c %e %e tie %d activ %d\n",nodes,l,f_cs[idof1],f_cm[idof1],i+1,islavact[j]);
//	    }
	  }
	}	
      }      
    }	      
    for( i=0; i<*ntie; i++){      	
      if(tieset[i*(81*3)+80]=='C'){		
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
    }    
    f_csn=sqrt(f_cs_tot[0]*f_cs_tot[0]+f_cs_tot[1]*f_cs_tot[1]+f_cs_tot[2]*f_cs_tot[2]);      
    f_cmn=sqrt(f_cm_tot[0]*f_cm_tot[0]+f_cm_tot[1]*f_cm_tot[1]+f_cm_tot[2]*f_cm_tot[2]);      
    printf("\nslave contact force : %e %e %e and norm: %e \n",f_cs_tot[0],f_cs_tot[1], f_cs_tot[2], f_csn);      
    printf("master contact force: %e %e %e and norm: %e \n",f_cm_tot[0],f_cm_tot[1], f_cm_tot[2], f_cmn );      
    free(f_cs_tot);free(f_cm_tot);  
    
    if(*iit>20){*iflagact=0;}    
//    if(*iit==2){FORTRAN(stop,());}
    printf("\n contacstress : N_stick : %d\t N_slip : %d\tN_Inactiv : %d\t N_nogap : %d  N_nolm : %d\n Flag = %d\n",nstick,nslip,ninacti,nnogap,nolm,*iflagact);
    fin= clock();
    printf("contactstress : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);
    
    free(bhat2);free(cstress2);free(u_old);
    return;
}
