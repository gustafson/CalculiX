/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2013 Guido Dhondt                     */

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
#include <time.h>
#include "CalculiX.h"

void matrixsort(double *au,int *mast1,int *irow, int *jq, 
	int *nzs,int *dim){
		
  int i, j,jj, k,kk,l,m, ll, kflag,numb;
  
  kflag = 2;
  FORTRAN(isortiid, (mast1, irow, au, nzs, &kflag));
  /*  fill in jqbd
      jqbd(i): first element in field aubd belonging to column i  */

  j = 0;
  for(i=0; i<*dim; i++){
    if(j == *nzs){
      for(k=i; k<*dim; k++) 
	jq[k] = *nzs+1;
      break;
    }
    
    
    if(mast1[j] != i+1){
      jq[i] = j+1;
      continue;
    }
    
    jq[i] = j+1;
    
    while(1){
      j++;
      if(j == *nzs) break;
      if(mast1[j] != i+1) break;
    }
  }
  
  jq[*dim] = *nzs + 1;

  /* Sorting of the rows*/
 for (i=0;i<*dim;i++){
    if(jq[i+1]-jq[i]>0){
   numb=jq[i+1]-jq[i]; 
   FORTRAN(isortid,(&irow[jq[i]-1],&au[jq[i]-1],&numb,&kflag));
    }
  }

  return;
}