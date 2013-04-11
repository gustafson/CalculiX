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

/**
 *Multiplication of (Bd)^T*A*Bd
*/

void multimortar(double *au, double *ad, int *irow, int *jq, int *nzs,
	   double *aubd, double *bdd, int *irowbd, int *jqbd, int *nzsbd,
	   double **aucp, double *adc, int **irowcp, int *jqc, int *nzsc,
           double *auqdt,int *irowqdt,int *jqqdt,int *nzsqdt,
	   int *neq,double *b, double *bhat){
	
  int i, j, k, l,m,n, icol, kflag,index,indexold, *irowi=NULL,ifree, 
  *irows=NULL, *irowc=NULL, *jqi=NULL, *jqs=NULL,nzsi, *irowbdt=NULL,
  *mast1=NULL,*ipointer=NULL,*irowqd=NULL,*jqqd=NULL,nzsqd;

  double *adi=NULL, *aui=NULL, *aus=NULL, value, *aux=NULL, *auc=NULL,
        *auqd=NULL,*unitmatrix=NULL;
 			 
  irowc = *irowcp; auc=*aucp;

    //Reference: Stefan Hueber's thesis page 28-33
    
	
  /* copy aubd into auqdt (transpose of qd) */

  for (j=0;j<*nzsbd;j++){
    auqdt[j]=aubd[j];
    irowqdt[j]=irowbd[j];
  }
  for(j=0;j<neq[1]+1;j++){
    jqqdt[j]=jqbd[j];
  }
  *nzsqdt=*nzsbd;
    
/*	Premultiply Bd by -Dd^(-1)
           before multiplication: Bd is stored in aubd,
                                  Dd is stored in bdd
	    after multiplication: -Dd^(-1).Bd is stored in auqdt
                                  The diagonal contains only 1's                      
*/

  for(j=0; j<neq[1]; j++){
    //    if (bdd[j]==0) continue; wrong!
    for(i=jqqdt[j]-1; i<jqqdt[j+1]-1; i++){
      icol=irowqdt[i]-1;
      auqdt[i] /= -bdd[icol];
      //printf("aubd[%d]=%e\n",i,aubd[i]);
    }
  }

  //   int number=1;

  //   FORTRAN(writematrix,(auqdt,bdd,irowqdt,jqqdt,&neq[1],&number));

  

  //   number=2;

  //  FORTRAN(writematrix,(au,ad,irow,jq,&neq[1],&number));
	
  /* determining the symmetry part of au and storing it
     into aus */

  ifree = 0;
  mast1=NNEW(int,nzs[1]);

  aus = NNEW(double,nzs[1]);
  irows = NNEW(int,nzs[1]);
  jqs = NNEW(int,neq[1]+1);
  
  for(j=0; j<neq[1]; j++){
    for(i= jq[j]-1; i< jq[j+1]-1; i++){
      mast1[i] = irow[i];
      irows[i] = j+1;
      aus[i] = au[i];
    }
  }    
    	
  // Order the matrix aus and determine jqs
  
  kflag = 2;
  FORTRAN(isortiid, (mast1, irows, aus, &nzs[1], &kflag));
  
  j=0;
  for(i=0; i<neq[1]; i++){
    if(j == nzs[1]){
      for(k=i; k<neq[1]; k++) 
	jqs[k] = nzs[1]+1;
      break;
    }
    
    if(mast1[j] != i+1){
      jqs[i] = j+1;
      continue;
    }
    
    jqs[i] = j+1;
    
    while(1){
      j++;
      if(j == nzs[1]) break;
      if(mast1[j] != i+1) break;
    }
  }
  jqs[neq[1]] = nzs[1] + 1;
  
  free(mast1);

  //  number=3;

  //FORTRAN(writematrix,(aus,ad,irows,jqs,&neq[1],&number));
	
  /* Calculate ai = A*Qdt
  **A = aus + au + ad
  **Qdt = auqdt + [unit]
  */
  
  ifree = 0;
  nzsi=nzs[1];
  mast1 = NNEW(int,nzsi);
  ipointer=NNEW(int,neq[1]);
  
  irowi = NNEW(int,nzsi);
  aui = NNEW(double,nzsi);
  adi = NNEW(double,neq[1]);
  jqi = NNEW(int,neq[1]+1);

  for(j=0; j<neq[1]; j++){
    m=j+1; // column number of the right matrix (FORTRAN)
    
    for(i= jqqdt[j]-1; i< jqqdt[j+1]-1; i++){
      
      icol = irowqdt[i]-1; //column number of the left matrix for multiplication (C)

      //au*auqdt

      for(k=jq[icol]-1; k<jq[icol+1]-1; k++){ 
	l = irow[k];
	value = au[k]*auqdt[i];
	if(l==m) adi[j]+=value;
	else{
	  insertas(ipointer, &irowi, &mast1, &l, &m, 
		   &ifree, &nzsi, &value, &aui);
	}
      }

      //ad*auqdt

      value = ad[icol]*auqdt[i];
      n=icol+1;
      if(n==m) printf("anormal !!!\n");
      else{
	insertas(ipointer, &irowi, &mast1, &n, &m, 
		 &ifree, &nzsi, &value, &aui);
      }

      //aus*auqdt

      for(k=jqs[icol]-1; k<(jqs[icol+1]-1); k++){ 
	l = irows[k];
	value = aus[k]*auqdt[i];
	if(l==m) adi[j]+=value;
	else{
	  insertas(ipointer, &irowi, &mast1, &l, &m,
		   &ifree, &nzsi, &value, &aui);}
      }
      
    }
    
    //au*[unit]
    
    for(i=jq[j]-1; i<jq[j+1]-1; i++){
      value = au[i];
      n=irow[i];
      if(n==m) {
	printf("au add actually impossible\n");
	adi[j]+=value;}
      else{
	insertas(ipointer, &irowi, &mast1, &n, &m, 
		 &ifree, &nzsi, &value, &aui);
      }	
    }
    
    //aus*[unit]
    
    for(i=jqs[j]-1; i<jqs[j+1]-1; i++){ 
      value = aus[i];
      n=irows[i];
      if(n==m) {
	printf("aus bdd actually impossible\n");
	adi[j]+=value;}
      else{
	insertas(ipointer, &irowi, &mast1, &n, &m, 
		 &ifree, &nzsi, &value, &aui);
      }	
    }
    
    //ad*[unit]
    
    adi[j] += ad[j];
    
  }
  nzsi=ifree;
	
  // Order the matrix aui and determine jqi
  
  //replace mast1 by the column numbers
  
  for(i=0; i<neq[1]; i++){
    if(ipointer[i]==0)
      continue;
    else{
      indexold = ipointer[i];
      while(1){
	index = mast1[indexold-1];
	mast1[indexold-1] = i+1;
	if(index == 0){break;}
	else{
	  indexold = index;
	}
      }
    }
  }    
  
  kflag = 2;
  FORTRAN(isortiid, (mast1, irowi, aui, &nzsi, &kflag));
  
  j=0;
  for(i=0; i<neq[1]; i++){
    if(j == nzsi){
      for(k=i; k<neq[1]; k++) 
	jqi[k] = nzsi+1;
      break;
    }
    
    if(mast1[j] != i+1){
      jqi[i] = j+1;
      continue;
    }
    
    jqi[i] = j+1;
    
    while(1){
      j++;
      if(j == nzsi) break;
      if(mast1[j] != i+1) break;
    }
    
  }
  jqi[neq[1]] = nzsi + 1;	

   //int  number=4;

   //FORTRAN(writematrix,(aui,adi,irowi,jqi,&neq[1],&number));
  
  free(mast1);free(ipointer);

  /* transpose Qdt and storing it in Qd */ 
  
  nzsqd=*nzsqdt;
  auqd=NNEW(double,nzsqd);
  irowqd=NNEW(int,nzsqd);
  mast1=NNEW(int,nzsqd);
  jqqd=NNEW(int,neq[1]+1);
  
  for(j=0;j<nzsqd;j++){auqd[j]=auqdt[j];}

  /* mast1 contains the rows of the original matrix,
     irowbt the columns */
  
  for(j=0;j<neq[1];j++){
    for(i=jqqdt[j]-1;i<jqqdt[j+1]-1;i++){
      mast1[i]=irowqdt[i];
      irowqd[i]=j+1;
    }
  }
  
  kflag = 2;
  FORTRAN(isortiid, (mast1, irowqd, auqd, &nzsqd, &kflag));
  
  j=0;
  for(i=0; i<neq[1]; i++){
    
    if(j == nzsqd){
      for(k=i; k<neq[1]; k++) //Guido fiche <= 21.04.09
	jqqd[k] = nzsqd+1;
      break;
    }
    
    if(mast1[j] != i+1){
      jqqd[i] = j+1;
      continue;
    }
    
    jqqd[i] = j+1;
    
    while(1){
      j++;
      if(j == nzsqd) break;
      if(mast1[j] != i+1) break;
    }
  }
  jqqd[neq[1]] = nzsqd + 1;	
  
  free(mast1);
	
  // calculate bhat = Qd*b, the new rhs
  
  unitmatrix=NNEW(double,neq[1]);
  for(j=0;j<neq[1];j++){unitmatrix[j]=1.;}

     //number=5;

   //FORTRAN(writematrix,(auqd,unitmatrix,irowqd,jqqd,&neq[1],&number));
  
  FORTRAN(opnonsym, (&neq[1], aux, b, bhat, unitmatrix, auqd, jqqd, irowqd));

  free(unitmatrix);

  /* Calculating Ac = Qd.Ai */

  ifree = 0;
  mast1 = NNEW(int,*nzsc);
  ipointer = NNEW(int,neq[1]);

  // auqd * aui
  
  for(j=0; j<neq[1]; j++){
    m=j+1;
    for(i= jqi[j]-1; i< jqi[j+1]-1; i++){
      icol = irowi[i]-1;
      for(k=jqqd[icol]-1; k<jqqd[icol+1]-1; k++){
	l = irowqd[k];
	value = auqd[k]*aui[i];
	if (l==m) adc[j]+=value;
	else{
	  insertas(ipointer, &irowc, &mast1, &l, &m, 
		   &ifree, nzsc, &value, &auc);
	}
      }
      // [unit] * aui
      
      l=icol+1;
      value = aui[i]; 
      if (l==m) adc[j]+=value;
      else{
	insertas(ipointer, &irowc, &mast1, &l, &m, 
		 &ifree, nzsc, &value, &auc);
      }
    }
    
    //auqd*adi
    
    for(i=jqqd[j]-1; i<jqqd[j+1]-1; i++){
      value = auqd[i]*adi[j];
      n=irowqd[i];
      if (n==m) adc[j]+=value;
      else{
	insertas(ipointer, &irowc, &mast1, &n, &m, 
		 &ifree, nzsc, &value, &auc);
      }
    }
    
    //[unit]*adi
    
    adc[j] += adi[j]; 
    
  }
  *nzsc=ifree;

  free(auqd);free(irowqd);free(jqqd);

  // Order the matrix auc and determine jqc
  
  //replace mast1 by the column numbers
  
  for(i=0; i<neq[1]; i++){
    if(ipointer[i]==0)
      continue;
    else{
      indexold = ipointer[i];
      while(1){
	index = mast1[indexold-1];
	mast1[indexold-1] = i+1;
	if(index == 0){break;}
	else{
	  indexold = index;
	}
      }
    }
  }
  
  /* ordering the column numbers in mast1 */

  kflag = 2;
  FORTRAN(isortiid, (mast1, irowc, auc,nzsc, &kflag));
  
  j=0;
  for(i=0; i<neq[1]; i++){
    if(j == *nzsc){
      for(k=i; k<neq[1]; k++) //Guido fiche <= 21.04.09
	jqc[k] = *nzsc+1;
      break;
    }
    
    if(mast1[j] != i+1){
      jqc[i] = j+1;
      continue;
    }
    
    jqc[i] = j+1;
    
    while(1){
      j++;
      if(j == *nzsc) break;
      if(mast1[j] != i+1) break;
    }
  }
  jqc[neq[1]] = *nzsc + 1;

     //number=6;

   //FORTRAN(writematrix,(auc,adc,irowc,jqc,&neq[1],&number));

  free(mast1);free(ipointer);
 
  /* free intermediate matrix */

  free(aui);free(adi);free(irowi);free(jqi);
 
  /* free symmetric part of au */

  free(aus);free(irows);free(jqs);

  *irowcp = irowc; *aucp=auc;

  return;
}
