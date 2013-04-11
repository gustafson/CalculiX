/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2007 Guido Dhondt                     */

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
    int* iflagact,double* cstress, int *mi, double *cdisp, double *alambda,
    double *alambdad,int *iit){
    
    int i,j,idof1,idof2,idof3,nodes,mt=mi[1]+1,nacti=0,ninacti=0,nnogap=0;

    double aux,stressnormal,stressnormal2,dispnormal,*unitmatrix=NULL,
        constant=10000000000.E1,lamb_1,lamb_2,lamb_3,lambd_1,lambd_2,lambd_3;

	clock_t debut;
	clock_t fin;

	
    /* determining the contact stress vectors and updating the active
       and inactive sets */

    int number=11;
    *iflagact=0;

//  FORTRAN(writematrix,(auc,b,irowc,jqc,&neq[1],&number));
	
//   for(j=0;j<neq[1];j++){

//      printf("contactstress bdd[%d]=%e\n",j,bdd[j]);
//      printf("contactstress uhat b[%d]=%e\n",j,b[j]);
//      printf("contactstress bhat b[%d]=%e\n",j,bhat[j]);
//      printf("%e\n",j,b[j]);
//    }
    
	debut = clock();
    FORTRAN(op,(&neq[1],&aux,b,cstress,adc,auc,icolc,irowc,nzlc));

//   for(j=0;j<neq[1];j++){

//   printf("contactstress cstress cstress[%d]=%e\n",j,cstress[j]);
//    }

    for (i=0;i<*ntie;i++){
      for(j=nslavnode[i];j<nslavnode[i+1];j++){
	nodes=islavnode[j];

        /* mechanical degrees of freedom */

	idof1=nactdof[mt*nodes-3]-1;
	idof2=nactdof[mt*nodes-2]-1;
	idof3=nactdof[mt*nodes-1]-1;

	/* calculation of the Lagrange multiplier
           (= contact pressure) */

	cstress[idof1]=(bhat[idof1]-cstress[idof1])/bdd[idof1];
	cstress[idof2]=(bhat[idof2]-cstress[idof2])/bdd[idof2];
	cstress[idof3]=(bhat[idof3]-cstress[idof3])/bdd[idof3];
	
/*
	alambda[idof1]+=cstress[idof1];
	alambda[idof2]+=cstress[idof2];
	alambda[idof3]+=cstress[idof3];

	alambdad[idof1]+=cstress[idof1]*bdd[idof1];
	alambdad[idof2]+=cstress[idof2]*bdd[idof2];
	alambdad[idof3]+=cstress[idof3]*bdd[idof3];
*/

	lamb_1=alambda[idof1]+cstress[idof1];
	lamb_2=alambda[idof2]+cstress[idof2];
	lamb_3=alambda[idof3]+cstress[idof3];

	lambd_1=alambdad[idof1]+cstress[idof1]*bdd[idof1];
	lambd_2=alambdad[idof2]+cstress[idof2]*bdd[idof2];
	lambd_3=alambdad[idof3]+cstress[idof3]*bdd[idof3];
	

	/* scaled normal stress */

//	stressnormal2=cstress[idof1]*slavnor[3*j]+cstress[idof2]*slavnor[3*j+1]+
//	             cstress[idof3]*slavnor[3*j+2];
/*
	stressnormal=alambda[idof1]*slavnor[3*j]+alambda[idof2]*slavnor[3*j+1]+
	             alambda[idof3]*slavnor[3*j+2];

	stressnormal2=alambdad[idof1]*slavnor[3*j]+alambdad[idof2]*slavnor[3*j+1]+
	             alambdad[idof3]*slavnor[3*j+2];
*/


	stressnormal=lamb_1*slavnor[3*j]+lamb_2*slavnor[3*j+1]+
	             lamb_3*slavnor[3*j+2];

	stressnormal2=lambd_1*slavnor[3*j]+lambd_2*slavnor[3*j+1]+
	             lambd_3*slavnor[3*j+2];



				 
//	printf("contact stress Lambda_n[%d]=%e, cstress_n[%d]=%e, gap[%d]=%e\n",
//              nodes,stressnormal,nodes,stressnormal2,nodes,gap[j]);

	/* Division by bdd is anulated by multiplication by bdd needed for
	stressnormal (scaled) */
	
	dispnormal=b[idof1]*slavnor[3*j]
	          +b[idof2]*slavnor[3*j+1]
	          +b[idof3]*slavnor[3*j+2];
	
//	if(stressnormal+constant*(dispnormal-gap[j])>1E-13){
	if(*iit<=8){
	if (islavact[j]!=-1){
//	printf("Contactstress Node %d val = %e\n",nodes,stressnormal+constant*(dispnormal-gap[j]));
	if(stressnormal+constant*(dispnormal-gap[j])>-1E-10){
	    nacti++;
	    if (islavact[j]!=1) {*iflagact = 1;
//				 printf("%d change\n",nodes);
	    }
//	    printf("%d ACTIF\n",nodes);
	    islavact[j]=1;
            cdisp[6*j]=dispnormal;
	    cdisp[6*j+3]=stressnormal;
//            printf("node : %d stress=%e\n",nodes,stressnormal);
	}else{
	    if (islavact[j]!=0){ *iflagact = 1;
//				 printf("%d INACTIF\n",nodes);
	    }
	    ninacti++;
	    islavact[j]=0;
//	    printf("%d INACTIF\n",nodes);
	    cstress[idof1]=0.;
	    cstress[idof2]=0.;
	    cstress[idof3]=0.;
	    alambda[idof1]=0.;
	    alambda[idof2]=0.;
	    alambda[idof3]=0.;
	    alambdad[idof1]=0.;
	    alambdad[idof2]=0.;
	    alambdad[idof3]=0.;
            cdisp[6*j]=0.;
	    cdisp[6*j+3]=0.;
	}
	}else{
//	    printf("%d NON_ACTIF\n",nodes);
            nnogap++;
	    islavact[j]=0;
	    cstress[idof1]=0.;
	    cstress[idof2]=0.;
	    cstress[idof3]=0.;
	    alambda[idof1]=0.;
	    alambda[idof2]=0.;
	    alambda[idof3]=0.;
	    alambdad[idof1]=0.;
	    alambdad[idof2]=0.;
	    alambdad[idof3]=0.;
	    cdisp[6*j]=0.;
	    cdisp[6*j+3]=0.;
	}

	if(*iflagact==0){ //Update Lambda if active set not changed
	alambda[idof1]=lamb_1;
	alambda[idof2]=lamb_2;
	alambda[idof3]=lamb_3;

	alambdad[idof1]=lambd_1;
	alambdad[idof2]=lambd_2;
	alambdad[idof3]=lambd_3;
	}


      }else{ //ACTIF and INACTIF remain the same just update datas of activ set
	switch(islavact[j]){
	
        case 0 : //active
                        ninacti++;
	case 1 : //active
			{
		    	nacti++;
//		    	islavact[j]=1;
	            	cdisp[6*j]=dispnormal;
		    	cdisp[6*j+3]=stressnormal;
			break;
			}
	default : //inactiv
			{
			nnogap++;
//	    		islavact[j]=0;
//		    	printf("%d INACTIF\n",nodes);
//	    		cstress[idof1]=0.;
//	    		cstress[idof2]=0.;
//	    		cstress[idof3]=0.;
//	    		alambda[idof1]=0.;
//	    		alambda[idof2]=0.;
//	    		alambda[idof3]=0.;
//	    		alambdad[idof1]=0.;
//	    		alambdad[idof2]=0.;
//	    		alambdad[idof3]=0.;
//          		cdisp[6*j]=0.;
//	    		cdisp[6*j+3]=0.;
			break;
			}

	}

	}
    }

	}

    /* transforming the constrained displacements into the standard displacements */

	
/*
   for(j=0;j<nslavnode[*ntie];j++){

   printf("contactstress islavact[%d]=%d\n",j,islavact[j]);
    } */

    unitmatrix=NNEW(double,neq[1]);
    for(j=0;j<neq[1];j++){
      unitmatrix[j]=1.;
      bhat[j]=b[j];
    }

    FORTRAN(opnonsym, (&neq[1], &aux, bhat, b, unitmatrix, auqdt, jqqdt, irowqdt));
    
/*    for(j=0;j<neq[1];j++){
      printf("contactstress u b[%d]=%e\n",j,b[j]);
    }
*/

  number=12;

//  FORTRAN(writematrix,(auc,b,irowc,jqc,&neq[1],&number));

    free(unitmatrix);
	printf("\ncontacstress : N_Activ : %d\tN_Inactiv : %d\tN_nogap : %d\n",nacti,ninacti,nnogap);
	fin= clock();
	printf("contactress : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);

    return;
}
