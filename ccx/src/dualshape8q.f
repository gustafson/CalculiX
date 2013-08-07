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
      subroutine dualshape8q(xi,et,xl,xsj,xs,shp,iflag)
!
!     shape functions and derivatives for a 8-node quadratic
!     isoparametric quadrilateral element. -1<=xi,et<=1 
!
!     iflag=2: calculate the value of the shape functions,
!              their derivatives w.r.t. the local coordinates
!              and the Jacobian vector (local normal to the
!              surface)
!     iflag=3: calculate the value of the shape functions, the
!              value of their derivatives w.r.t. the global
!              coordinates and the Jacobian vector (local normal
!              to the surface)
!
      implicit none
!
      integer i,j,k,iflag
!
      real*8 shp(4,8),xs(3,2),xsi(2,3),xl(0:3,8),sh(3),xsj(3)
!
      real*8 xi,et
!
!     shape functions and their glocal derivatives for an element
!     described with two local parameters and three global ones.
!
!     local derivatives of the shape functions: xi-derivative
!
      shp(1,1)=(1.d0-et)*(2.d0*xi+et)/4.d0
      shp(1,2)=(1.d0-et)*(2.d0*xi-et)/4.d0
      shp(1,3)=(1.d0+et)*(2.d0*xi+et)/4.d0
      shp(1,4)=(1.d0+et)*(2.d0*xi-et)/4.d0
      shp(1,5)=-xi*(1.d0-et)
      shp(1,6)=(1.d0-et*et)/2.d0
      shp(1,7)=-xi*(1.d0+et)
      shp(1,8)=-(1.d0-et*et)/2.d0
!
!     local derivatives of the shape functions: eta-derivative
!
      shp(2,1)=(1.d0-xi)*(2.d0*et+xi)/4.d0
      shp(2,2)=(1.d0+xi)*(2.d0*et-xi)/4.d0
      shp(2,3)=(1.d0+xi)*(2.d0*et+xi)/4.d0
      shp(2,4)=(1.d0-xi)*(2.d0*et-xi)/4.d0
      shp(2,5)=-(1.d0-xi*xi)/2.d0
      shp(2,6)=-et*(1.d0+xi)
      shp(2,7)=(1.d0-xi*xi)/2.d0
      shp(2,8)=-et*(1.d0-xi)
!
!     standard shape functions
!
      shp(3,1)=(1.d0-xi)*(1.d0-et)*(-xi-et-1.d0)/4.d0
      shp(3,2)=(1.d0+xi)*(1.d0-et)*(xi-et-1.d0)/4.d0
      shp(3,3)=(1.d0+xi)*(1.d0+et)*(xi+et-1.d0)/4.d0
      shp(3,4)=(1.d0-xi)*(1.d0+et)*(-xi+et-1.d0)/4.d0
      shp(3,5)=(1.d0-xi*xi)*(1.d0-et)/2.d0
      shp(3,6)=(1.d0+xi)*(1.d0-et*et)/2.d0
      shp(3,7)=(1.d0-xi*xi)*(1.d0+et)/2.d0
      shp(3,8)=(1.d0-xi)*(1.d0-et*et)/2.d0
!
!     Dual shape functions
!
      shp(4,1)=shp(3,1)+(shp(3,5)+shp(3,8))/5.d0
      shp(4,2)=shp(3,2)+(shp(3,5)+shp(3,6))/5.d0
      shp(4,3)=shp(3,3)+(shp(3,6)+shp(3,7))/5.d0
      shp(4,4)=shp(3,4)+(shp(3,7)+shp(3,8))/5.d0
!
!     computation of the local derivative of the global coordinates
!     (xs)
!
      do i=1,3
        do j=1,2
          xs(i,j)=0.d0
          do k=1,8
            xs(i,j)=xs(i,j)+xl(i,k)*shp(j,k)
          enddo
        enddo
      enddo
!
!     computation of the jacobian vector
!
      xsj(1)=xs(2,1)*xs(3,2)-xs(3,1)*xs(2,2)
      xsj(2)=xs(1,2)*xs(3,1)-xs(3,2)*xs(1,1)
      xsj(3)=xs(1,1)*xs(2,2)-xs(2,1)*xs(1,2)
!
      if(iflag.eq.2) return
!
!     computation of the global derivative of the local coordinates
!     (xsi) (inversion of xs)
!
c      xsi(1,1)=xs(2,2)/xsj(3)
c      xsi(2,1)=-xs(2,1)/xsj(3)
c      xsi(1,2)=-xs(1,2)/xsj(3)
c      xsi(2,2)=xs(1,1)/xsj(3)
c      xsi(1,3)=-xs(2,2)/xsj(1)
c      xsi(2,3)=xs(2,1)/xsj(1)
      if(dabs(xsj(3)).gt.1.d-10) then
         xsi(1,1)=xs(2,2)/xsj(3)
         xsi(2,2)=xs(1,1)/xsj(3)
         xsi(1,2)=-xs(1,2)/xsj(3)
         xsi(2,1)=-xs(2,1)/xsj(3)
         if(dabs(xsj(2)).gt.1.d-10) then
            xsi(2,3)=xs(1,1)/xsj(2)
            xsi(1,3)=-xs(1,2)/xsj(2)
         elseif(dabs(xsj(1)).gt.1.d-10) then
            xsi(2,3)=xs(2,1)/xsj(1)
            xsi(1,3)=-xs(2,2)/xsj(1)
         else
            xsi(2,3)=0.d0
            xsi(1,3)=0.d0
         endif
      elseif(dabs(xsj(2)).gt.1.d-10) then
         xsi(1,1)=xs(3,2)/xsj(2)
         xsi(2,3)=xs(1,1)/xsj(2)
         xsi(1,3)=-xs(1,2)/xsj(2)
         xsi(2,1)=-xs(3,1)/xsj(2)
         if(dabs(xsj(1)).gt.1.d-10) then
            xsi(1,2)=xs(3,2)/xsj(1)
            xsi(2,2)=-xs(3,1)/xsj(1)
         else
            xsi(1,2)=0.d0
            xsi(2,2)=0.d0
         endif
      else
         xsi(1,2)=xs(3,2)/xsj(1)
         xsi(2,3)=xs(2,1)/xsj(1)
         xsi(1,3)=-xs(2,2)/xsj(1)
         xsi(2,2)=-xs(3,1)/xsj(1)
         xsi(1,1)=0.d0
         xsi(2,1)=0.d0
      endif
!
!     computation of the global derivatives of the shape functions
!
      do k=1,8
        do j=1,3
          sh(j)=shp(1,k)*xsi(1,j)+shp(2,k)*xsi(2,j)
        enddo
        do j=1,3
          shp(j,k)=sh(j)
        enddo
      enddo
!
      return
      end
