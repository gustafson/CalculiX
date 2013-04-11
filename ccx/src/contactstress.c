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

void contactstress(double *bhat, double *adc, double *auc,int *jqc, 
    int *irowc, int *neq, double *gap, double *bdd, double *b, int *islavact,
    double *auqdt, int *irowqdt, int *jqqdt, int *ntie, int *nslavnode,
    int *islavnode, double *slavnor, int *icolc, int *nzlc, int *nactdof, 
    int* iflagact,double* cstress, int *mi, double *cdisp, double *alm_old,
    int *iit){
    
    int i,j,idof1,idof2,idof3,nodes,mt=mi[1]+1,nacti=0,ninacti=0,nnogap=0,
        max_node;

    double aux,stressnormal,stressnormal2,dispnormal,*unitmatrix=NULL,
        constant=10000000000.E1,lamb_1,lamb_2,lamb_3,
        *alm_new=NULL,max=0.0,lm_res,*delta_lm=NULL;

	clock_t debut;
	clock_t fin;

	
    /* determining the contact stress vectors and updating the active
       and inactive sets and the Langrange Multipliers (LM) */

    int number=11;
    *iflagact=0;

    
    debut = clock();
    FORTRAN(op,(&neq[1],&aux,b,cstress,adc,auc,icolc,irowc,nzlc));

    alm_new=NNEW(double,neq[1]); // Contains the new current value of LM
    delta_lm=NNEW(double,neq[1]);
    
    for (i=0;i<*ntie;i++){
      for(j=nslavnode[i];j<nslavnode[i+1];j++){
	nodes=islavnode[j];

        /* mechanical degrees of freedom */

	idof1=nactdof[mt*nodes-3]-1;
	idof2=nactdof[mt*nodes-2]-1;
	idof3=nactdof[mt*nodes-1]-1;

	/* calculation of the Lagrange multiplier
           (= contact pressure) */
 
        //delta 
	cstress[idof1]=(bhat[idof1]-cstress[idof1])/bdd[idof1];
	cstress[idof2]=(bhat[idof2]-cstress[idof2])/bdd[idof2];
	cstress[idof3]=(bhat[idof3]-cstress[idof3])/bdd[idof3];

//	cstress[idof1]=(bhat[idof1])/bdd[idof1];
//	cstress[idof2]=(bhat[idof2])/bdd[idof2];
//	cstress[idof3]=(bhat[idof3])/bdd[idof3];


	//global
	alm_new[idof1]=alm_old[idof1]+cstress[idof1];
	alm_new[idof2]=alm_old[idof2]+cstress[idof2];
	alm_new[idof3]=alm_old[idof3]+cstress[idof3];

	/* calculation of the LM residual */
	lm_res=(cstress[idof1]*cstress[idof1]+cstress[idof2]*cstress[idof2]+
		cstress[idof3]*cstress[idof3])/
		(alm_new[idof1]*alm_new[idof1]+alm_new[idof2]*alm_new[idof2]+
		alm_new[idof3]*alm_new[idof3]);

	/* scaled normal stress */

	stressnormal=alm_new[idof1]*slavnor[3*j]+alm_new[idof2]*slavnor[3*j+1]+
	             alm_new[idof3]*slavnor[3*j+2];

	/* Division by bdd is anulated by multiplication by bdd needed for
	stressnormal (scaled) */
	
	dispnormal=b[idof1]*slavnor[3*j]
	          +b[idof2]*slavnor[3*j+1]
	          +b[idof3]*slavnor[3*j+2];
	
	if(*iit<=8){ //After 8 iterations the active and inactive sets are frozen
	if (islavact[j]!=-1){
        
        /*** Active/Inactive set condition cf: Hueeber Phd p.72 ***/
	if(stressnormal+constant*(dispnormal-gap[j])>-1E-10){
	    nacti++;
	    if (islavact[j]!=1) {*iflagact = 1;
	    }
	    islavact[j]=1;
            cdisp[6*j]=dispnormal;
	    cdisp[6*j+3]=stressnormal;
	    if (lm_res>max){
                 max=lm_res;
                 max_node=nodes;  
            } 	
	}else{
	    if (islavact[j]!=0){ *iflagact = 1;
	    }
	    ninacti++;
	    islavact[j]=-1;
	    cstress[idof1]=0.;
	    cstress[idof2]=0.;
	    cstress[idof3]=0.;
	    alm_old[idof1]=0.;
	    alm_old[idof2]=0.;
	    alm_old[idof3]=0.;
	    alm_new[idof1]=0.;
	    alm_new[idof2]=0.;
	    alm_new[idof3]=0.;
            cdisp[6*j]=0.;
	    cdisp[6*j+3]=0.;
	}
	}else{
            nnogap++;
	    cstress[idof1]=0.;
	    cstress[idof2]=0.;
	    cstress[idof3]=0.;
	    alm_old[idof1]=0.;
	    alm_old[idof2]=0.;
	    alm_old[idof3]=0.;
	    alm_new[idof1]=0.;
	    alm_new[idof2]=0.;
	    alm_new[idof3]=0.;
	    cdisp[6*j]=0.;
	    cdisp[6*j+3]=0.;
	}

      }else{ //ACTIF and INACTIF remain the same just update datas of activ set
	switch(islavact[j]){
	
        case 0 : //Inactive
                        ninacti++;
		    	islavact[j]=-1;
  		    	cstress[idof1]=0.;
			cstress[idof2]=0.;
	    		cstress[idof3]=0.;
	    		alm_old[idof1]=0.;
	    		alm_old[idof2]=0.;
	    		alm_old[idof3]=0.;
            		cdisp[6*j]=0.;
	    		cdisp[6*j+3]=0.;
            		break;
	case 1 : //Active
			{
		    	nacti++;
//		    	islavact[j]=1;
			alm_old[idof1]=alm_new[idof1];
			alm_old[idof2]=alm_new[idof2];
			alm_old[idof3]=alm_new[idof3];
	            	cdisp[6*j]=dispnormal;
		    	cdisp[6*j+3]=stressnormal;
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
	    		alm_old[idof1]=0.;
	    		alm_old[idof2]=0.;
	    		alm_old[idof3]=0.;
          		cdisp[6*j]=0.;
	    		cdisp[6*j+3]=0.;
			break;
			}

	}

	}
    }

	}

	if((*iflagact==0)&&(*iit<=8)){
//	if((*iit<=8)){
   for (i=0;i<*ntie;i++){
      for(j=nslavnode[i];j<nslavnode[i+1];j++){
	nodes=islavnode[j];

        /* mechanical degrees of freedom */

	idof1=nactdof[mt*nodes-3]-1;
	idof2=nactdof[mt*nodes-2]-1;
	idof3=nactdof[mt*nodes-1]-1; //Update LM if active set not changed
	alm_old[idof1]=alm_new[idof1];
	alm_old[idof2]=alm_new[idof2];
	alm_old[idof3]=alm_new[idof3];
	}
     }
    }	

    /* transforming the constrained displacements and the langrange multipliers 
       into the standard displacements */

    unitmatrix=NNEW(double,neq[1]);
    delta_lm=NNEW(double,neq[1]);
    for(j=0;j<neq[1];j++){
      unitmatrix[j]=1.;
      bhat[j]=b[j];
     }

    FORTRAN(opnonsym, (&neq[1], &aux, bhat, b, unitmatrix, auqdt, jqqdt, irowqdt));
   

    free(unitmatrix);
    free(delta_lm);
	printf("\ncontacstress : N_Activ : %d\tN_Inactiv : %d\tN_nogap : %d\n Flag = %d, RES= %e in node %d\n",nacti,ninacti,nnogap,*iflagact,max,max_node);
	fin= clock();
	printf("contactress : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);
    free(alm_new);
    return;
}
