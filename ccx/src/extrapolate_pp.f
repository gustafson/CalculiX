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
      subroutine extrapolate_pp(nface,ielfa,xrlfa,ppel,ppfa,
     &  ifabou)
!
!     inter/extrapolation of p' element values to the faces
!
      implicit none
!
      integer nface,ielfa(3,*),ifabou(*),i,iel1,iel2
!
      real*8 xrlfa(3,*),ppel(*),ppfa(*),xl1
!     
      do i=1,nface
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
            ppfa(i)=xl1*ppel(iel1)+xrlfa(2,i)*ppel(iel2)
         elseif(iel2.lt.0) then
            if(ifabou(-iel2)+4.ne.0) then
!
!              p given
!               
               ppfa(i)=0.d0
            else
!
!              extrapolation
!
               ppfa(i)=xl1*ppel(iel1)+xrlfa(3,i)*ppel(ielfa(3,i))
            endif
         else
!
!           extrapolation
!
            ppfa(i)=xl1*ppel(iel1)+xrlfa(3,i)*ppel(ielfa(3,i))
         endif
      enddo
!            
      return
      end
