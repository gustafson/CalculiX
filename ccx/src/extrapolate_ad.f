!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2014 Guido Dhondt
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
      subroutine extrapolate_ad(nface,ielfa,xrlfa,ad,adfa,nactdoh)
!
!     inter/extrapolation of ad at the center of the elements
!     to the center of the faces in the form 1/ad
!
      implicit none
!
      integer nface,ielfa(3,*),ipo1,iel2,ipo3,i,nactdoh(*)
!
      real*8 xrlfa(3,*),xl1,adfa(*),ad(*)
!
      do i=1,nface
         ipo1=nactdoh(ielfa(1,i))
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
            adfa(i)=xl1/ad(ipo1)+xrlfa(2,i)/ad(nactdoh(iel2))
         else
            ipo3=nactdoh(ielfa(3,i))
            adfa(i)=xl1/ad(ipo1)+xrlfa(3,i)/ad(ipo3)
         endif
      enddo
!            
      return
      end
