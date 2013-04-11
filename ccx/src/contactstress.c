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
#include "CalculiX.h"

void contactstress(double *bhat, double *adc, double *auc,int *jqc, 
    int *irowc, int *neq, double *gap, double *bdd, double *b, int *islavact,
    double *auqdt, int *irowqdt, int *jqqdt, int *ntie, int *nslavnode,
    int *islavnode, double *slavnor, int *icolc, int *nzlc, int *nactdof, 
    int* iflagact,double* cstress, int *mi, double *lambda){
    
    int i,j,idof1,idof2,idof3,nodes,mt=mi[1]+1;

    double aux,stressnormal,stressnormal2,dispnormal,*unitmatrix=NULL,
        constant=1.E1;

    /* determining the contact stress vectors and updating the active
       and inactive sets */

    int number=11;
    *iflagact=0;

 // FORTRAN(writematrix,(auc,b,irowc,jqc,&neq[1],&number));
	
/*   for(j=0;j<neq[1];j++){

//      printf("contactstress bdd[%d]=%e\n",j,bdd[j]);
//      printf("contactstress uhat b[%d]=%e\n",j,b[j]);
//      printf("%e\n",j,b[j]);
    }*/
    
    FORTRAN(op,(&neq[1],&aux,b,cstress,adc,auc,icolc,irowc,nzlc));

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
	
	lambda[idof1]+=cstress[idof1];
	lambda[idof2]+=cstress[idof2];
	lambda[idof3]+=cstress[idof3];

	/* scaled normal stress */

	stressnormal2=cstress[idof1]*slavnor[3*j]+cstress[idof2]*slavnor[3*j+1]+
	             cstress[idof3]*slavnor[3*j+2];

	stressnormal=lambda[idof1]*slavnor[3*j]+lambda[idof2]*slavnor[3*j+1]+
	             lambda[idof3]*slavnor[3*j+2];
				 
    printf("contact stress Lambda_n[%d]=%e, cstress_n[%d]=%e\n",j,stressnormal,j,stressnormal2);

	/* Division by bdd is anulated by multiplication by bdd needed for
	stressnormal (scaled) */
	
	dispnormal=b[idof1]*slavnor[3*j]
	          +b[idof2]*slavnor[3*j+1]
	          +b[idof3]*slavnor[3*j+2];
	
	if(stressnormal+constant*(dispnormal-gap[j])>1E-10){
		if (islavact[j]!=1) {*iflagact = 1;
		//	printf("1 node %d, value = %f\n", j,stressnormal+constant*(dispnormal-gap[j]));
		}
		islavact[j]=1;
	}else{
		if (islavact[j]!=0){ *iflagact = 1;
		//		printf("2 node %d, value = %f\n", j,stressnormal+constant*(dispnormal-gap[j]));
		}
		islavact[j]=0.;
		lambda[idof1]=0.;
		lambda[idof2]=0.;
		lambda[idof3]=0.;
		cstress[idof1]=0.;
		cstress[idof2]=0.;
		cstress[idof3]=0.;
	}
//	printf("node %d, status :%d\n",j,islavact[j]);
      }
    }

    /* transforming the constrained displacements into the standard displacements */

    unitmatrix=NNEW(double,neq[1]);
    for(j=0;j<neq[1];j++){
      unitmatrix[j]=1.;
      bhat[j]=b[j];
//      printf("contactstress dof=%d,cstress=%e\n",j,cstress[j]);
    }

    FORTRAN(opnonsym, (&neq[1], &aux, bhat, b, unitmatrix, auqdt, jqqdt, irowqdt));
    
    for(j=0;j<neq[1];j++){
//      printf("contactstress u b[%d]=%e\n",j,b[j]);
    }


  number=12;

//  FORTRAN(writematrix,(auc,b,irowc,jqc,&neq[1],&number));

    free(unitmatrix);
    return;
}
