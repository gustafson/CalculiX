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
    int *islavnode, double *slavnor, int *icolc, int *nzlc, int *nactdof){
    
    int i,j,idof1,idof2,idof3,nodes;

    double *cstress=NULL,aux,stressnormal,dispnormal,*unitmatrix=NULL,
        constant=1.;

    /* determining the contact stress vectors and updating the active
       and inactive sets */

  int number=11;

  FORTRAN(writematrix,(auc,b,irowc,jqc,&neq[1],&number));


    cstress=NNEW(double,neq[1]);
    
    FORTRAN(op,(&neq[1],&aux,b,cstress,adc,auc,icolc,irowc,nzlc));

    for (i=0;i<*ntie;i++){
      for(j=nslavnode[i];j<nslavnode[i+1];j++){
	nodes=islavnode[j];

        /* mechanical degrees of freedom */

	idof1=nactdof[4*nodes-3]-1;
	idof2=nactdof[4*nodes-2]-1;
	idof3=nactdof[4*nodes-1]-1;

	/* calculation of the Lagrange multiplier
           (= contact pressure) */

	cstress[idof1]=(bhat[idof1]-cstress[idof1])/bdd[idof1];
	cstress[idof2]=(bhat[idof2]-cstress[idof2])/bdd[idof2];
	cstress[idof3]=(bhat[idof3]-cstress[idof3])/bdd[idof3];

	/* scaled normal stress */

	stressnormal=cstress[idof1]*bdd[idof1]+cstress[idof2]*bdd[idof2]+
	             cstress[idof3]*bdd[idof3];

	/* Division by bdd is anulated by multiplication by bdd needed for
	stressnormal (scaled) */
	
	//	stressnormal=(bhat[idof1]-cstress[idof1])/bdd[idof1]*slavnor[3*j]
	//          +(bhat[idof2]-cstress[idof2])/bdd[idof2]*slavnor[3*j+1]
	//          +(bhat[idof3]-cstress[idof3])/bdd[idof3]*slavnor[3*j+2];
	
	//stressnormal=(bhat[idof1]-cstress[idof1])*slavnor[3*j]
	//            +(bhat[idof2]-cstress[idof2])*slavnor[3*j+1]
	//            +(bhat[idof3]-cstress[idof3])*slavnor[3*j+2];
	
	dispnormal=b[idof1]*slavnor[3*j]
	          +b[idof2]*slavnor[3*j+1]
	          +b[idof3]*slavnor[3*j+2];
	
	if(stressnormal+constant*(dispnormal-gap[j])>0.){
	  islavact[j]=1;
	}else{
	  islavact[j]=0.;
	}
      }
    }

    /* transforming the constrained displacements into the standard displacements */

    unitmatrix=NNEW(double,neq[1]);
    for(j=0;j<neq[1];j++){
      unitmatrix[j]=1.;
      bhat[j]=b[j];
      printf("contactstress dof=%d,cstress=%e\n",j,cstress[j]);
    }

    FORTRAN(opnonsym, (&neq[1], &aux, bhat, b, unitmatrix, auqdt, jqqdt, irowqdt));
    

  number=12;

  FORTRAN(writematrix,(auc,b,irowc,jqc,&neq[1],&number));

    free(unitmatrix);

    return;
}
