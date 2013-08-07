!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2013 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine smooth(adb,aub,adl,sol,aux,icol,irow,jq,
     &  neq,nzl,csmooth)
!
!     smoothing the finite element solution
!
!     Ref: The Finite Element Method for Fluid Dynamics,
!          O.C. Zienkiewicz, R.L. Taylor & P. Nithiarasu
!          6th edition (2006) ISBN 0 7506 6322 7
!          p. 61
!
      implicit none
!
      integer icol(*),irow(*),jq(*),neq,nzl,i,j,k
!
      real*8 adb(*),aub(*),adl(*),sol(*),aux(*),p,csmooth,c1,c2
!
!     multiplying the original matrix with zero
!     diagonal with the actual solution 
!
      call op(neq,p,sol,aux,adb,aub,icol,irow,nzl)
c      call op(neq,p,sol,aux,adl,aub,icol,irow,nzl)
!
!     lumping the matrix set (adb,aux) and storing the resulting
!     diagonal terms in adl
!
      do i=1,neq
         adl(i)=adb(i)
      enddo
!
      do j=1,neq
         do k=jq(j),jq(j+1)-1
            i=irow(k)
            adl(i)=adl(i)+aub(k)
            adl(j)=adl(j)+aub(k)
         enddo
      enddo
!
!     determining the multiplicative constants
!
      c2=1.d0+csmooth/2.d0
      c1=1.d0/c2
      c2=csmooth/c2
!
!     smoothing the solution
!
      do i=1,neq
         sol(i)=c1*sol(i)+c2*aux(i)/adl(i)/2.d0
c         sol(i)=c1*sol(i)+c2*aux(i)/adl(i)
      enddo
!
      return
      end
