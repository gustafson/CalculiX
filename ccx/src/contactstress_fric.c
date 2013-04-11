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

void contactstress_fric(double *bhat, double *adc, double *auc,int *jqc, 
    int *irowc, int *neq, double *gap, double *bdd,double *aubd, 
    int *jqbd, int *irowbd, double *b, int *islavact,
    double *auqdt, int *irowqdt, int *jqqdt, int *ntie, int *nslavnode,
    int *islavnode, int *nmastnode, int *imastnode, double *slavnor,double *slavtan,
    int *icolc, int *nzlc, int *nactdof,int* iflagact,double* cstress, int *mi,
    double *cdisp, double *f_cs, double *f_cm, int *iit,int *iwatchactiv,double *vold,double* bp,int *nk,double *friccoeff){
    
    int i,j,l,idof1,idof2,idof3,nodes,mt=mi[1]+1,nstick=0,nslip=0,ninacti=0,nnogap=0,
        max_node,debug;
    
    double aux,stressnormal,ddispnormal,*unitmatrix=NULL,uholdt1,uholdt2,bpold,
        constant=1.E10,stresst1,stresst2,ddispt1,ddispt2,helpa,helpb,
        *alm_new=NULL,max=0.0,lm_res,*delta_lm=NULL,*bhat2=NULL,*cstressold=NULL,
        f_csn, f_cmn, *f_cs_tot=NULL, *f_cm_tot=NULL, *vectornull=NULL,*w_t=NULL,nw_t,
        disp_totalnormal,disp_totalt1,disp_totalt2,*auqdt2=NULL,*u_old=NULL;
    
    clock_t debut;
    clock_t fin;
    
    debug=1;
//    constant=7.E5;
    *iflagact=0;
    
    
    cstressold=NNEW(double,neq[1]);
    bhat2=NNEW(double,neq[1]);
    w_t=NNEW(double,2);
    debut = clock();
    for(j=0;j<neq[1];j++){cstressold[j]=cstress[j];}   
    FORTRAN(opnonsym,(&neq[1],&aux,b,cstress,adc,auc,jqc,irowc));
    alm_new=NNEW(double,neq[1]); // Contains the new current value of LM
    delta_lm=NNEW(double,neq[1]);
    
    bhat2=NNEW(double,neq[1]);
    auqdt2=NNEW(double,jqqdt[neq[1]]);
    unitmatrix=NNEW(double,neq[1]);
    u_old=NNEW(double,neq[1]);

    for(j=1;j<*nk+1;j++){
	    idof1=nactdof[mt*j-3]-1;
	    idof2=nactdof[mt*j-2]-1;
	    idof3=nactdof[mt*j-1]-1; 
	    if(idof1>-1) u_old[idof1]=vold[mt*j-3];
	    if(idof2>-1) u_old[idof2]=vold[mt*j-2];
	    if(idof3>-1) u_old[idof3]=vold[mt*j-1];
	    if(idof1>-1 && idof2>-1 && idof3>-1){
//	    printf("u(%d)=%e %e %e \n",j,u_old[idof1],u_old[idof2],u_old[idof3]);
	    }
    } 
    for(j=0;j<neq[1];j++){
	    unitmatrix[j]=1.;
    }
    for(j=0;j<jqqdt[neq[1]];j++){ 
           auqdt2[j]=-auqdt[j];
    }      
	FORTRAN(opnonsym, (&neq[1], &aux, u_old, bhat2, unitmatrix, auqdt2, jqqdt, irowqdt));
	
	free(unitmatrix);free(auqdt2);
    for (i=0;i<*ntie;i++){
	for(j=nslavnode[i];j<nslavnode[i+1];j++){
	    nodes=islavnode[j];
	    
	    /* mechanical degrees of freedom */
	    
	    idof1=nactdof[mt*nodes-3]-1;
	    idof2=nactdof[mt*nodes-2]-1;
	    idof3=nactdof[mt*nodes-1]-1;
//	    printf("u(%d)=%e %e %e \n",j,u_old[idof1],u_old[idof2],u_old[idof3]);
//	    printf("uhold(%d)= %e %e %e\n",nodes,bhat2[idof1],bhat2[idof2],bhat2[idof3]);
	}
    }	
    
    for (i=0;i<*ntie;i++){
	for(j=nslavnode[i];j<nslavnode[i+1];j++){
	    nodes=islavnode[j];
	    
	    /* mechanical degrees of freedom */
	    
	    idof1=nactdof[mt*nodes-3]-1;
	    idof2=nactdof[mt*nodes-2]-1;
	    idof3=nactdof[mt*nodes-1]-1;
	    
	    /* calculation of the Lagrange multiplier
	       (= contact pressure) */
	    
	    if(bdd[idof1]>0.0){
		cstress[idof1]=(bhat[idof1]-cstress[idof1])/bdd[idof1];
		cstress[idof2]=(bhat[idof2]-cstress[idof2])/bdd[idof2];
		cstress[idof3]=(bhat[idof3]-cstress[idof3])/bdd[idof3];
	    }else{
		cstress[idof1]=(bhat[idof1]-cstress[idof1]);
		cstress[idof2]=(bhat[idof2]-cstress[idof2]);
		cstress[idof3]=(bhat[idof3]-cstress[idof3]);
	    }

	    /* scaled normal stress */

	    stressnormal=cstress[idof1]*slavnor[3*j]+cstress[idof2]*slavnor[3*j+1]+
		cstress[idof3]*slavnor[3*j+2];
	    stresst1=cstress[idof1]*slavtan[6*j]+cstress[idof2]*slavtan[6*j+1]+
		cstress[idof3]*slavtan[6*j+2];
	    stresst2=cstress[idof1]*slavtan[6*j+3]+cstress[idof2]*slavtan[6*j+4]+
		cstress[idof3]*slavtan[6*j+5];	
		
		
	    
	    /* Division by bdd is anulated by multiplication by bdd needed for
	       stressnormal (scaled) */
	    
	    ddispnormal=b[idof1]*slavnor[3*j]
		+b[idof2]*slavnor[3*j+1]
		+b[idof3]*slavnor[3*j+2];
	    ddispt1=b[idof1]*slavtan[6*j]
		+b[idof2]*slavtan[6*j+1]
		+b[idof3]*slavtan[6*j+2];
	    ddispt2=b[idof1]*slavtan[6*j+3]
		+b[idof2]*slavtan[6*j+4]
		+b[idof3]*slavtan[6*j+5];
	   disp_totalnormal=(b[idof1]+bhat2[idof1])*slavnor[3*j]
		+(b[idof2]+bhat2[idof2])*slavnor[3*j+1]
		+(b[idof3]+bhat2[idof3])*slavnor[3*j+2];
	   disp_totalt1=(b[idof1]+bhat2[idof1])*slavtan[6*j]
		+(b[idof2]+bhat2[idof2])*slavtan[6*j+1]
		+(b[idof3]+bhat2[idof3])*slavtan[6*j+2];
	   disp_totalt2=(b[idof1]+bhat2[idof1])*slavtan[6*j+3]
		+(b[idof2]+bhat2[idof2])*slavtan[6*j+4]
		+(b[idof3]+bhat2[idof3])*slavtan[6*j+5];

	   uholdt1=(bhat2[idof1])*slavtan[6*j]
		+(bhat2[idof2])*slavtan[6*j+1]
		+(bhat2[idof3])*slavtan[6*j+2];
	   uholdt2=(bhat2[idof1])*slavtan[6*j+3]
		+(bhat2[idof2])*slavtan[6*j+4]
		+(bhat2[idof3])*slavtan[6*j+5];
	   bpold=bp[j];
//	   if(islavact[j]==1){
//	     helpa=disp_totalt1+friccoeff*uholdt1*stressnormal/bpold;
//	     helpb=disp_totalt2+friccoeff*uholdt2*stressnormal/bpold;
//	     printf("node %d h1 %e=%e h2 %e=%e?\n",nodes,helpa,uholdt1,helpb,uholdt2);
//	   }
//	   if(islavact[j]==2){
//	     helpa=sqrt(stresst1*stresst1+stresst2*stresst2);
//	     printf("node %d ||lm_t|| %e=%e lm_n ?\n",nodes,helpa,friccoeff*stressnormal);
//	   }	
	   bp[j]=friccoeff[j]*(stressnormal+constant*(disp_totalnormal-gap[j]));
	   w_t[0]=stresst1;
	   w_t[1]=stresst2;
	    nw_t=sqrt( w_t[0]* w_t[0]+w_t[1]*w_t[1]);
//           printf("\t lm_t1= %e lm_t2= %e nw_t= %e \n",w_t[0],w_t[1],nw_t);
           w_t[0]=stresst1+constant*disp_totalt1;
	   w_t[1]=stresst2+constant*disp_totalt2;
	   nw_t=sqrt( w_t[0]* w_t[0]+w_t[1]*w_t[1]);
	    if(debug==1 && islavact[j]>0){		
//	  if(debug==1 && ((islavact[j]==0 &&stressnormal+constant*(dispnormal-gap[j])>-1E-10) ||(islavact[j]==1 &&stressnormal+constant*(dispnormal-gap[j])<-1E-10)|| islavact[j]>0)){  
		printf("u(%d): %d %d %d \n", nodes, idof1,idof2, idof3);
		printf("\t uold_t1= %e uold_t2= %e \n",uholdt1,uholdt2);
//		printf("\t fhat(%d): %e, k*uhat: %e, bdd: %e \n",nodes, bhat[idof1], cstress[idof1],bdd[idof1]);
//		printf("\t fhat(%d): %e, k*uhat: %e, bdd: %e \n",nodes, bhat[idof2], cstress[idof2],bdd[idof2]);
//		printf("\t fhat(%d): %e, k*uhat: %e, bdd: %e \n",nodes, bhat[idof3], cstress[idof3],bdd[idof3]);
		printf("\t duhat(%d) : %e %e %e \n", nodes, b[idof1],b[idof2], b[idof3]);
		printf("\t duhat2(%d): %e %e %e \n", nodes, ddispnormal,ddispt1,ddispt2);
		printf("\t cstress(%d): %e %e %e, actif: %d  \n",nodes,cstress[idof1],cstress[idof2],cstress[idof3],islavact[j]);
		printf("\t lm(%d)     : %e %e %e \n",nodes, stressnormal, stresst1,stresst2);
		printf("\t dispnormal(%d)= %e, gap(%d)=%e bp=%e\n",nodes,disp_totalnormal,nodes,gap[j],bp[j]);
	        printf("\t uh_t1= %e uh_t2= %e nuh_t= %e \n",disp_totalt1,disp_totalt2,sqrt(disp_totalt1*disp_totalt1+disp_totalt2*disp_totalt2)) ;		
		printf("\t w_t1= %e w_t2= %e nw_t= %e \n",w_t[0],w_t[1],nw_t);
//           w_t[0]=stresst1+constant*disp_totalt1;
//	   w_t[1]=stresst2+constant*disp_totalt2;
//	   nw_t=sqrt( w_t[0]* w_t[0]+w_t[1]*w_t[1]);
//	   	printf("\t w_totalt1= %e w_totalt2= %e nw_totalt= %e \n",w_t[0],w_t[1],nw_t);
	    }
	    
	    if(*iit<=10){ //After 8 iterations the active and inactive sets are frozen
		if (islavact[j]!=-1){
		    
		    /*** Active/Inactive set condition cf: Hueeber Phd p.72 ***/
		    if(bp[j]>-1E-10 && nw_t<bp[j]){
			nstick++;
//			if(islavact[j]>0){bp[j]=0.0;}
			if(debug==1)printf("\t stick!\n");
			if (islavact[j]!=1) {*iflagact = 1;
			}
			islavact[j]=1;
			cdisp[6*j]=1;
//			cdisp[6*j+1]=disp_totalt1;
			cdisp[6*j+1]=(sqrt(stresst1*stresst1+stresst2*stresst2));
			cdisp[6*j+2]=bp[j];			
//			cdisp[6*j+2]=disp_totalt2;
			cdisp[6*j+3]=stressnormal;
			cdisp[6*j+4]=stresst1;
			cdisp[6*j+5]=stresst2;
//			if (lm_res>max){
//			    max=lm_res;
//			    max_node=nodes;  
//			} 
		    }else if(bp[j]>-1E-10 && nw_t>bp[j]){
			nslip++;
//			if(islavact[j]>0){bp[j]=0.0;}
			if(debug==1)printf("\t slip!\n");
			if (islavact[j]!=2) {*iflagact = 1;
			}
			islavact[j]=2;
			cdisp[6*j]=2;
//			cdisp[6*j+1]=disp_totalt1;
			cdisp[6*j+1]=(sqrt(stresst1*stresst1+stresst2*stresst2));			
//			cdisp[6*j+2]=disp_totalt2;
			cdisp[6*j+2]=bp[j];			
			cdisp[6*j+3]=stressnormal;
			cdisp[6*j+4]=stresst1;
			cdisp[6*j+5]=stresst2;
//			if (lm_res>max){
//			    max=lm_res;
//			    max_node=nodes;  
//			} 	
		    }else{
			if (islavact[j]>0){ *iflagact = 1;
			}
			ninacti++;
			if(debug==1)printf("\t inactiv!\n");
			islavact[j]=0;
			cstress[idof1]=0.;
			cstress[idof2]=0.;
			cstress[idof3]=0.;
		        cdisp[6*j]=0.;
		        cdisp[6*j+1]=0.;
		        cdisp[6*j+2]=0.;
		        cdisp[6*j+3]=0.;
		        cdisp[6*j+4]=0.;
		        cdisp[6*j+5]=0.;
		    }
		}else{
		    nnogap++;
		    if(debug==1)printf("\t nogap!\n");
		    cstress[idof1]=0.;
		    cstress[idof2]=0.;
		    cstress[idof3]=0.;
		    cdisp[6*j]=0.;
		    cdisp[6*j+1]=0.;
		    cdisp[6*j+2]=0.;
		    cdisp[6*j+3]=0.;
		    cdisp[6*j+4]=0.;
		    cdisp[6*j+5]=0.;
		}
		
	    }else{ //ACTIF and INACTIF remain the same just update datas of activ set
		switch(islavact[j]){
		    
		    case 0 : //Inactive
		    {
                        ninacti++;
  		    	cstress[idof1]=0.;
			cstress[idof2]=0.;
	    		cstress[idof3]=0.;
            		cdisp[6*j]=0.;
	    		cdisp[6*j+1]=0.;
            		cdisp[6*j+2]=0.;
	    		cdisp[6*j+3]=0.;
            		cdisp[6*j+4]=0.;
	    		cdisp[6*j+5]=0.;
            		break;
		    }
		    case 1 : //stick
		    {
		    	nstick++;			
            		cdisp[6*j]=disp_totalnormal;
	    		cdisp[6*j+1]=disp_totalt1;
            		cdisp[6*j+2]=disp_totalt2;
	    		cdisp[6*j+3]=stressnormal;
            		cdisp[6*j+4]=stresst1;
	    		cdisp[6*j+5]=stresst2;
			if (lm_res>max){
			    max=lm_res;
			    max_node=nodes;  
 		        } 		
			break;
		    }
		    case 2 : //slip
		    {
		    	nslip++;			
            		cdisp[6*j]=disp_totalnormal;
	    		cdisp[6*j+1]=disp_totalt1;
            		cdisp[6*j+2]=disp_totalt2;
	    		cdisp[6*j+3]=stressnormal;
            		cdisp[6*j+4]=stresst1;
	    		cdisp[6*j+5]=stresst2;
			if (lm_res>max){
			    max=lm_res;
			    max_node=nodes;  
 		        } 		
			break;
		    }		    
		    default : //No gap
		    {
			nnogap++;
	    		cstress[idof1]=0.;
	    		cstress[idof2]=0.;
	    		cstress[idof3]=0.;
            		cdisp[6*j]=0.;
	    		cdisp[6*j+1]=0.;
            		cdisp[6*j+2]=0.;
	    		cdisp[6*j+3]=0.;
            		cdisp[6*j+4]=0.;
	    		cdisp[6*j+5]=0.;
			break;
		    }
		}
	    }
//	    gap[j]=gap[j]-ddispnormal;    
	}
    }
//    iwatchactiv[(*iit-1)*2+0]=nacti;
//    iwatchactiv[(*iit-1)*2+1]=*iflagact;
    
    /* transforming the constrained displacements and the langrange multipliers 
       into the standard displacements */
    

     unitmatrix=NNEW(double,neq[1]);
     for(j=0;j<neq[1];j++){
	  unitmatrix[j]=1.;
	  bhat2[j]=b[j];
      }        
      FORTRAN(opnonsym, (&neq[1], &aux, bhat2, b, unitmatrix, auqdt, jqqdt, irowqdt));
//    if(debug==1){  
//    for (i=0;i<*ntie;i++){
//	for(j=nslavnode[i];j<nslavnode[i+1];j++){
//	    nodes=islavnode[j];
	    
	    /* mechanical degrees of freedom */
	    
//	    idof1=nactdof[mt*nodes-3]-1;
//	    idof2=nactdof[mt*nodes-2]-1;
//	    idof3=nactdof[mt*nodes-1]-1;
//	    printf("u(%d)=%e %e %e \n",nodes,b[idof1],b[idof2],b[idof3]);    
//	}
//    }
//    }
	
      free(unitmatrix);
    
    	      /* calculating the master contact forces 
		 Ph.D. Thesis Stefan Hartmann eqn. (6.26) */
  

      for (i=0;i<neq[1];i++){
	  f_cs[i]=cstress[i]*bdd[i];
	  f_cm[i]=0.0;
      }	
      vectornull=NNEW(double,neq[1]);
      
      FORTRAN(opnonsymt,(&neq[1],&aux,cstress,f_cm,vectornull,aubd,jqbd,irowbd));
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
	    f_cs_tot[l]+=f_cs[idof1];
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
      f_csn=sqrt(f_cs_tot[0]*f_cs_tot[0]+f_cs_tot[1]*f_cs_tot[1]+f_cs_tot[2]*f_cs_tot[2]);
      f_cmn=sqrt(f_cm_tot[0]*f_cm_tot[0]+f_cm_tot[1]*f_cm_tot[1]+f_cm_tot[2]*f_cm_tot[2]);
      printf("\nslave contact force : %e %e %e and norm: %e \n",f_cs_tot[0],f_cs_tot[1], f_cs_tot[2], f_csn);
      printf("master contact force: %e %e %e and norm: %e \n",f_cm_tot[0],f_cm_tot[1], f_cm_tot[2], f_cmn );
      free(f_cs_tot);free(f_cm_tot);
    
    printf("\n contacstress : N_stick : %d\t N_slip : %d\tN_Inactiv : %d\t N_nogap : %d\n Flag = %d\n",nstick,nslip,ninacti,nnogap,*iflagact);
    fin= clock();
    printf("contactstress : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);
    free(alm_new);free(bhat2);
 
    
    
    return;
}
