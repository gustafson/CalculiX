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

void insertas(int *ipointer, int **irowp, int **mast1p, int *i1,
	 int *i2, int *ifree, int *nzs_, double *contribution, double **bdp){

  /*   inserts a new nonzero matrix position into the data structure
       the structure is not assumed to be symmetric 
       i1: row number (FORTRAN convention) 
       i2: column number (FORTRAN convention) */

  int idof1,idof2,istart,*irow=NULL,*mast1=NULL;
  double *bd=NULL;

  irow=*irowp;   
  mast1=*mast1p;
  bd=*bdp;

  idof1 = *i1;
  idof2 = *i2;
    
  if(ipointer[idof2-1]==0){
    ++*ifree;
    if(*ifree>*nzs_){
      *nzs_=(int)(1.1**nzs_);
      RENEW(irow,int,*nzs_);
      RENEW(mast1,int,*nzs_);
      RENEW(bd,double,*nzs_);
    }
    ipointer[idof2-1]=*ifree;
    irow[*ifree-1]=idof1;
    mast1[*ifree-1]=0;
    bd[*ifree-1]=*contribution;
  }
  else{
    istart=ipointer[idof2-1];
    while(1){
      if(irow[istart-1]==idof1){ //Former row number -> update value
	bd[istart-1]+=*contribution;
	break;
      }
      if(mast1[istart-1]==0){ // new row number value
	++*ifree;
	if(*ifree>*nzs_){
	  *nzs_=(int)(1.1**nzs_);
	  RENEW(bd,double,*nzs_);
	  RENEW(mast1,int,*nzs_);
	  RENEW(irow,int,*nzs_);
	}
	mast1[istart-1]=*ifree;
	irow[*ifree-1]=idof1;
	mast1[*ifree-1]=0;
	bd[*ifree-1]=*contribution;  //first value
	
	break;
      }
      else{
	istart=mast1[istart-1];
      }
    }
  }
  
  *irowp=irow;
  *mast1p=mast1;
  *bdp=bd;

  return;
  
}
