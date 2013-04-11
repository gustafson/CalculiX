/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2007 Guido Dhondt                     */

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

/*
 *Calculate the entries of Bd and Dd, and insert them into the data structure
*/

void bdfill(int **irowbdp, int *jqbd,
        double **aubdp, double *bdd,int *nzsbd, int *ntie, int *ipkon, int *kon, 
        char *lakon, int *nslavnode, int *nmastnode, int *imastnode,
        int *islavnode, int *islavsurf, int *imastsurf, double *pmastsurf, 
        int *itiefac, int *neq, int *nactdof, double *co, double *vold,
	int *iponoels, int *inoels, int *mi, double *gapmints, double *gap,
        double* pslavsurf,double* pslavdual){
		
  int i, j, k,l,m, idof1,idofs,idofm, nodes, nodem, kflag, index, indexold,
      *mast1=NULL, *irowbd=NULL,ifree,*ipointer=NULL,mt=mi[1]+1,contint;
  
  double contribution=0.0, *aubd=NULL,gapcont=0.0;
  
  irowbd = *irowbdp; aubd=*aubdp;
  
  ifree = 0;
  ipointer=NNEW(int,neq[1]);
  mast1=NNEW(int,*nzsbd);
  
  /* calculating the off-diagonal terms and storing them in aubd */

  /* meaning of the fields in FORTRAN notation:
     ipointer(i): points to an element in field aubd belonging to column i
                  aubd(ipointer(i)): value of that element
                  irowbd(ipointer(i)): row to which that element belongs
     mast1(ipointer(i)): points to another element in field aubd belonging
                         to column i, unless zero.
  */
  
  for( i=0; i<*ntie; i++){
    for(j=nslavnode[i]; j<nslavnode[i+1]; j++){
      nodes = islavnode[j];
      for(k=nmastnode[i]; k<nmastnode[i+1]; k++){
	nodem = imastnode[k];
	FORTRAN(createbdentry, ( &i,ipkon,kon,lakon,&nodem,&nodes,
	     islavsurf,imastsurf,pmastsurf,itiefac,&contribution,co,vold,
	     iponoels,inoels,mi,pslavsurf,pslavdual));
	contribution=-contribution;
//	printf("S : %d, M: %d, %e\n",nodes,nodem,contribution);
	for(l=0; l<3; l++){
	  idofs = nactdof[mt*(nodes-1)+l+1];
	  //	  for(m=0;m<3;m++){	
	    idofm = nactdof[mt*(nodem-1)+l+1];						
	    if ((idofs>0)&&(idofm>0)){ //insertion for active dofs
	      insertas(ipointer, &irowbd, &mast1, &idofs, &idofm, &ifree, nzsbd,
		       &contribution, &aubd);
	     // printf("idofs =%d , ifofm =%d, cont = %e\n",idofs,idofm,contribution);
	    }
	    //}	
	}
      }
    }
  }
  
  *nzsbd=ifree;
  RENEW(irowbd,int,*nzsbd);
  RENEW(aubd,double,*nzsbd);
  
  //replace mast1 by the column numbers
  
  for(i=0; i<neq[1]; i++){
    if(ipointer[i]==0)
      continue;
    else{
      indexold = ipointer[i];
      while(1){
	index = mast1[indexold-1];
	mast1[indexold-1] = i+1;
	if(index == 0){
	  break;}
	else{
	  indexold = index;
	}
      }
    }
  }
  
  
  /* Sort mast1, irowbd and aubd; 
     Outcome: the values in field aubd are sorted, column by
     column; no sorting is done within the columns */

  kflag = 2;
  FORTRAN(isortiid, (mast1, irowbd, aubd, nzsbd, &kflag));
  
  
  /*  fill in jqbd
      jqbd(i): first element in field aubd belonging to column i  */

  j = 0;
  for(i=0; i<neq[1]; i++){
    if(j == *nzsbd){
      for(k=i; k<neq[1]; k++) 
	jqbd[k] = *nzsbd+1;
      break;
    }
    
    
    if(mast1[j] != i+1){
      jqbd[i] = j+1;
      continue;
    }
    
    jqbd[i] = j+1;
    
    while(1){
      j++;
      if(j == *nzsbd) break;
      if(mast1[j] != i+1) break;
    }
  }
  
  jqbd[neq[1]] = *nzsbd + 1;

  free(ipointer);free(mast1);
  
  /* determining the diagonal entries and storing them in bdd */

  for(i=0;i<neq[1];i++){bdd[i]=0.;}
  
  for(i=0; i<*ntie; i++){
    for(j=nslavnode[i]; j<nslavnode[i+1]; j++){
      nodes = islavnode[j];
      FORTRAN(createddentry,(&i,ipkon,kon,&nodes,
      	      lakon,islavsurf,itiefac,&contribution,co,vold,
	      iponoels,inoels,mi,&gapcont,gapmints,pslavdual));
      gap[j]=gapcont/contribution;  
//      printf(" gap[%d]= %e\n",j,gap[j]);
//        gap[j]=0.0;
      for(l=0; l<3; l++){
	idof1 = nactdof[mt*(nodes-1)+l+1];
	if (idof1>0)
	  bdd[idof1-1]+=contribution;
      }
    }
  }
  
  *irowbdp = irowbd; *aubdp=aubd;
  
  return;
}
