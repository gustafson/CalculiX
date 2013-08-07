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
#include "CalculiX.h"

void insertas(int **irowp, int **mast1p, int *i1,
	 int *i2, int *ifree, int *nzs_, double *contribution, double **bdp){

  /*   inserts a new nonzero matrix position into the data structure
       the structure is not assumed to be symmetric 
       i1: row number (FORTRAN convention) 
       i2: column number (FORTRAN convention) */

  int idof1,idof2,*irow=NULL,*mast1=NULL;
  double *bd=NULL;

  irow=*irowp;   
  mast1=*mast1p;
  bd=*bdp;

  idof1 = *i1;
  idof2 = *i2;

    if(*ifree>*nzs_){
//      printf("Insertas RENEW ifree = %d,nzs = %d\n",*ifree,*nzs_);
//      *nzs_=(int)(1.1**nzs_);
      *nzs_=(int)(1.5**nzs_);
      RENEW(irow,int,*nzs_);
      if (irow==NULL) printf("WARNING !!!!\n");
      RENEW(mast1,int,*nzs_);
      RENEW(bd,double,*nzs_);
    }
    mast1[*ifree-1]=idof2;
    irow[*ifree-1]=idof1;
    bd[*ifree-1]=*contribution;
    ++*ifree;
//    printf("ifree %d\n",*ifree);
  
  *irowp=irow;
  *mast1p=mast1;
  *bdp=bd;

  return;
  
}
