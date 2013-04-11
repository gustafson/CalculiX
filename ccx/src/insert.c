/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2007 Guido Dhondt                          */

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

void insert(int *ipointer, int **mast1p, int **mast2p, int *i1,
	    int *i2, int *ifree, int *nzs_){

  /*   inserts a new nonzero matrix position into the data structure */

  int idof1,idof2,istart,*mast1=NULL,*mast2=NULL;

  mast1=*mast1p;
  mast2=*mast2p;

  if(*i1<*i2){
    idof1=*i1;
    idof2=*i2;
  }
  else{
    idof1=*i2;
    idof2=*i1;
  }

  if(ipointer[idof2-1]==0){
    ++*ifree;
    if(*ifree>*nzs_){
      *nzs_=(int)(1.1**nzs_);
      RENEW(mast1,int,*nzs_);
      RENEW(mast2,int,*nzs_);
      /*     printf(" reallocation: nzs_=%d\n\n",*nzs_);*/
    }
    ipointer[idof2-1]=*ifree;
/*    printf("idof1=%d,idof2=%d,ifree=%d\n",idof1,idof2,*ifree);*/
    mast1[*ifree-1]=idof1;
    mast2[*ifree-1]=0;
  }
  else{
    istart=ipointer[idof2-1];
    while(1){
      if(mast1[istart-1]==idof1) break;
      if(mast2[istart-1]==0){
	++*ifree;
	if(*ifree>*nzs_){
	  *nzs_=(int)(1.1**nzs_);
	  RENEW(mast1,int,*nzs_);
	  RENEW(mast2,int,*nzs_);
/*	  printf(" reallocation: nzs_=%d\n\n",*nzs_);*/
	}
	mast2[istart-1]=*ifree;
	mast1[*ifree-1]=idof1;
	mast2[*ifree-1]=0;
	break;
      }
      else{
	istart=mast2[istart-1];
      }
    }
  }

  *mast1p=mast1;
  *mast2p=mast2;
  
  return;

}
	  
/*

Here starts the original FORTRAN code, which was transferred to the
C-code above in order to allow automatic reallocation

      subroutine insert(ipointer,mast1,mast2,i1,i2,ifree,nzs_)
!
!     inserts a new nonzero matrix position into the data structure
!
      implicit none
!
      integer ipointer(*),mast1(*),mast2(*),i1,i2,ifree,nzs_,idof1,
     &  idof2,istart
!
      if(i1.lt.i2) then
        idof1=i1
        idof2=i2
      else
        idof1=i2
        idof2=i1
      endif
!
      if(ipointer(idof2).eq.0) then
        ifree=ifree+1
        if(ifree.gt.nzs_) then
           write(*,*) '*ERROR in insert: increase nzs_'
           stop
        endif
        ipointer(idof2)=ifree
        mast1(ifree)=idof1
        mast2(ifree)=0
      else
        istart=ipointer(idof2)
        do
          if(mast1(istart).eq.idof1) exit
          if(mast2(istart).eq.0) then
            ifree=ifree+1
            if(ifree.gt.nzs_) then
               write(*,*) '*ERROR in insert: increase nzs_'
               stop
            endif
            mast2(istart)=ifree
            mast1(ifree)=idof1
            mast2(ifree)=0
            exit
          else
            istart=mast2(istart)
          endif
        enddo
      endif
!
      return
      end

      */
