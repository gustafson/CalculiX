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
      subroutine extrapolate_gradv(nface,ielfa,xrlfa,gradv,gradvfa)
!
!     inter/extrapolation of p element values to the faces
!
      implicit none
!
      integer nface,ielfa(3,*),i,iel1,iel2,k,l
!
      real*8 xrlfa(3,*),gradv(3,3,*),gradvfa(3,3,*),xl1
!     
      do i=1,nface
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.ne.0) then
            do k=1,3
               do l=1,3
                  gradvfa(k,l,i)=xl1*gradv(k,l,iel1)+
     &                           xrlfa(2,i)*gradv(k,l,iel2)
               enddo
            enddo
         else
            do k=1,3
               do l=1,3
                  gradvfa(k,l,i)=xl1*gradv(k,l,iel1)+
     &                           xrlfa(3,i)*gradv(k,l,ielfa(3,i))
               enddo
            enddo
         endif
      enddo
!            
      return
      end
