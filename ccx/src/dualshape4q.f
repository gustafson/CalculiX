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
      subroutine dualshape4q(xi,et,xl,xsj,xs,shp,ns,pslavdual,iflag)
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
      integer i,j,k,iflag,ns
!
      real*8 shp(4,8),xs(3,2),xsi(2,3),xl(3,8),sh(3),xsj(3)
!
      real*8 xi,et,pslavdual(16,*)
!
!     shape functions and their glocal derivatives for an element
!     described with two local parameters and three global ones.
!
!     local derivatives of the shape functions: xi-derivative
!      
      shp(1,1)=-(1.d0-et)/4.d0
      shp(1,2)=(1.d0-et)/4.d0
      shp(1,3)=(1.d0+et)/4.d0
      shp(1,4)=-(1.d0+et)/4.d0
!
!     local derivatives of the shape functions: eta-derivative
!
      shp(2,1)=-(1.d0-xi)/4.d0
      shp(2,2)=-(1.d0+xi)/4.d0
      shp(2,3)=(1.d0+xi)/4.d0
      shp(2,4)=(1.d0-xi)/4.d0
!
!     standard shape functions
!
      shp(3,1)=(1.d0-xi)*(1.d0-et)/4.d0
      shp(3,2)=(1.d0+xi)*(1.d0-et)/4.d0
      shp(3,3)=(1.d0+xi)*(1.d0+et)/4.d0
      shp(3,4)=(1.d0-xi)*(1.d0+et)/4.d0
!
!     Dual shape functions
!
c      shp(4,1)=4.d0*shp(3,1)-2.d0*shp(3,2)+shp(3,3)-2.d0*shp(3,4)
c      shp(4,2)=4.d0*shp(3,2)-2.d0*shp(3,1)+shp(3,4)-2.d0*shp(3,3)
c      shp(4,3)=4.d0*shp(3,3)-2.d0*shp(3,2)+shp(3,1)-2.d0*shp(3,4)
c      shp(4,4)=4.d0*shp(3,4)-2.d0*shp(3,1)+shp(3,2)-2.d0*shp(3,3)
!
!    with Mass Matrix pslavdual
!
      shp(4,1)=pslavdual(1,ns)*shp(3,1)+pslavdual(2,ns)*shp(3,2)+
     & pslavdual(3,ns)*shp(3,3)+pslavdual(4,ns)*shp(3,4)
      shp(4,2)=pslavdual(5,ns)*shp(3,1)+pslavdual(6,ns)*shp(3,2)+
     & pslavdual(7,ns)*shp(3,3)+pslavdual(8,ns)*shp(3,4)
      shp(4,3)=pslavdual(9,ns)*shp(3,1)+pslavdual(10,ns)*shp(3,2)+
     & pslavdual(11,ns)*shp(3,3)+pslavdual(12,ns)*shp(3,4)
      shp(4,4)=pslavdual(13,ns)*shp(3,1)+pslavdual(14,ns)*shp(3,2)+
     & pslavdual(15,ns)*shp(3,3)+pslavdual(16,ns)*shp(3,4)
c       if(ns.eq.214)then
c        write(*,*)'xi',xi,'et',et 
c        write(*,*)pslavdual(13,ns),pslavdual(14,ns),
c     &  pslavdual(15,ns),pslavdual(16,ns)
c        write(*,*)shp(3,1),shp(3,2),shp(3,3),shp(3,4)
c       endif
!
!     computation of the local derivative of the global coordinates
!     (xs)
!
      do i=1,3
        do j=1,2
          xs(i,j)=0.d0
          do k=1,4
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
      xsi(1,1)=xs(2,2)/xsj(3)
      xsi(2,1)=-xs(2,1)/xsj(3)
      xsi(1,2)=-xs(1,2)/xsj(3)
      xsi(2,2)=xs(1,1)/xsj(3)
      xsi(1,3)=-xs(2,2)/xsj(1)
      xsi(2,3)=xs(2,1)/xsj(1)
!
!     computation of the global derivatives of the shape functions
!
      do k=1,4
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
