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
      subroutine extrapolate_p(nface,ielfa,xrlfa,vel,vfa,ipbou,
     &  ifabou,xfabou)
!
!     inter/extrapolation of p element values to the faces
!
      implicit none
!
      integer nface,ielfa(3,*),ipbou(*),ifabou(*),i,iel1,iel2
!
      real*8 xrlfa(3,*),vel(0:5,*),vfa(0:5,*),xfabou(*),xl1
!     
      do i=1,nface
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.ne.0) then
            vfa(4,i)=xl1*vel(4,iel1)+xrlfa(2,i)*vel(4,iel2)
         elseif(ipbou(i).ne.0) then
            if(ifabou(ipbou(i)+5).ne.0) then
!
!              p given
!               
               vfa(4,i)=xfabou(ipbou(i)+5)
            else
!
!              extrapolation
!
               vfa(4,i)=xl1*vel(4,iel1)+xrlfa(3,i)*vel(4,ielfa(3,i))
            endif
         else
!
!           extrapolation
!
            vfa(4,i)=xl1*vel(4,iel1)+xrlfa(3,i)*vel(4,ielfa(3,i))
         endif
      enddo
!            
      return
      end
