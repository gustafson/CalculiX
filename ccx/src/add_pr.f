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
      subroutine add_pr(au,ad,icol,jq,i,j,value,i0,i1)
!
!     stores coefficient (i,j) in the stiffness matrix stored in
!     profile format
!
      implicit none
!
      integer icol(*),jq(*),i,j,ii,jj,ipointer,i0,i1
      real*8 ad(*),au(*),value
!
      if(i.eq.j) then
         if(i0.eq.i1) then
            ad(i)=ad(i)+value
         else
            ad(i)=ad(i)+2.d0*value
         endif
         return
      elseif(i.gt.j) then
         ii=j
         jj=i
      else
         ii=i
         jj=j
      endif
!
      if(ii.lt.jq(jj)) then
         write(*,*) '*ERROR in add_pr: coefficient should be 0'
         stop
      else
         ipointer=icol(jj)-jj+ii+1
         au(ipointer)=au(ipointer)+value
      endif
!
      return
      end













