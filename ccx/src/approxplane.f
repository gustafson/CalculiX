!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998 Guido Dhondt
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
      subroutine approxplane(col,straight,xn)
!
!     calculate the equation of the planes through the
!     edges of a quadrilateral and parallel to the vector xn together
!     with a plane perpendicular to xn and through the center of gravity 
!     of the four corner nodes of the quadrilateral  
!     (so-called mean quadrilateral plane) with 
!     (col(1,1),col(2,1),col(3,1)),(col(1,2),col(2,2),col(3,2)),
!     (col(1,3),col(2,3),col(3,3)),(col(1,4),col(2,4),col(3,4))
!     as vertices. The equation of the planes through the first edge 
!     (connecting the first and the second node) is of the form
!     straight(1)*x+straight(2)*y+straight(3)*z+straight(4)=0, such that the
!     vector (straight(1),straight(2),straight(3)) points outwards (replace 
!     (1) by (5),(9) and (13) for the second, third and fourth edge,
!     similar offset for (2),(3) and (4); 
!     The equation of the mean quadrilateral plane is
!     straight(17)*x+straight(18)*y+straight(19)*z+straight(20)=0 such
!     that the quadrilateral is numbered clockwise when looking in the 
!     direction of vector (straight(17),straight(18),straight(19)).
!
      implicit none
!
      integer i
!
      real*8 col(3,4),straight(20),p12(3),p23(3),p34(3),p41(3),dd,xn(3)
!
!     sides of the quadrilateral
!
      do i=1,3
         p12(i)=col(i,2)-col(i,1)
         p23(i)=col(i,3)-col(i,2)
         p34(i)=col(i,4)-col(i,3)
         p41(i)=col(i,1)-col(i,4)
      enddo
!
!     mean normal to the quadrilateral (given)
!
      do i=17,19
         straight(i)=xn(i)
      enddo
!
!     p12 x xn
!
      straight(1)=p12(2)*xn(3)-p12(3)*xn(2)
      straight(2)=p12(3)*xn(1)-p12(1)*xn(3)
      straight(3)=p12(1)*xn(2)-p12(2)*xn(1)
      dd=dsqrt(straight(1)*straight(1)+straight(2)*straight(2)+
     &         straight(3)*straight(3))
      do i=1,3
         straight(i)=straight(i)/dd
      enddo
!
!     p23 x xn
!
      straight(5)=p23(2)*xn(3)-p23(3)*xn(2)
      straight(6)=p23(3)*xn(1)-p23(1)*xn(3)
      straight(7)=p23(1)*xn(2)-p23(2)*xn(1)
      dd=dsqrt(straight(5)*straight(5)+straight(6)*straight(6)+
     &         straight(7)*straight(7))
      do i=5,7
         straight(i)=straight(i)/dd
      enddo
!
!     p34 x xn
!
      straight(9)=p34(2)*xn(3)-p34(3)*xn(2)
      straight(10)=p34(3)*xn(1)-p34(1)*xn(3)
      straight(11)=p34(1)*xn(2)-p34(2)*xn(1)
      dd=dsqrt(straight(9)*straight(9)+straight(10)*straight(10)+
     &         straight(11)*straight(11))
      do i=9,11
         straight(i)=straight(i)/dd
      enddo
!
!     p41 x xn
!
      straight(13)=p41(2)*xn(3)-p41(3)*xn(2)
      straight(14)=p41(3)*xn(1)-p41(1)*xn(3)
      straight(15)=p41(1)*xn(2)-p41(2)*xn(1)
      dd=dsqrt(straight(13)*straight(13)+straight(14)*straight(14)+
     &         straight(15)*straight(15))
      do i=13,15
         straight(i)=straight(i)/dd
      enddo
!
!     determining the inhomogeneous terms
!
      straight(4)=-straight(1)*col(1,1)-straight(2)*col(2,1)-
     &             straight(3)*col(3,1)
      straight(8)=-straight(5)*col(1,2)-straight(6)*col(2,2)-
     &             straight(7)*col(3,2)
      straight(12)=-straight(9)*col(1,3)-straight(10)*col(2,3)-
     &             straight(11)*col(3,3)
      straight(16)=-straight(13)*col(1,4)-straight(14)*col(2,4)-
     &             straight(15)*col(3,4)
      straight(20)=-xn(1)*(col(1,1)+col(1,2)+col(1,3)+col(1,4))/4.d0
     &             -xn(2)*(col(2,1)+col(2,2)+col(2,3)+col(2,4))/4.d0
     &             -xn(3)*(col(3,1)+col(3,2)+col(3,3)+col(3,4))/4.d0
!
      return
      end

