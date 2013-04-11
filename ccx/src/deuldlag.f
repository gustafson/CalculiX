!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine deuldlag(xi,et,ze,xlag,xeul,xj,xs)
!
!     calculation of the coefficients of the linearization
!     of J:=det(dx/dX)=1 at (xi,et,ze) for a 20-node quadratic
!     isoparametric brick element. Watch out: 0<=xi<=2, -1<=et,ze<=1 !!!
!     X are Lagrangian coordinates, x are Eulerian coordinates
!
      implicit none
!
      integer i,j,k
!
      real*8 shp(3,20),xs(3,3),xlag(3,20),xeul(3,20)
!
      real*8 xi,et,ze,xj
!
!     shape functions and their glocal derivatives
!
!     local derivatives of the shape functions: xi-derivative
!
      shp(1, 1)=-0.125+et*(0.250+et*(-0.125-ze*0.125)
     &  +ze*(0.375+ze*0.125))+ze*(-0.250-ze*0.125)
     &  +xi*(0.250+et*(-0.250-0.250*ze)+ze*0.250)
      shp(1, 2)=-0.375+et*(0.250+et*(0.125+ze*0.125)
     &  +ze*(0.125-ze*0.125))+ze*(-0.250+ze*0.125)
     &  +xi*(0.250+et*(-0.250-0.250*ze)+ze*0.250)
      shp(1, 3)=-0.375+et*(0.250+et*(0.125-ze*0.125)
     &  +ze*(-0.125-ze*0.125))+ze*(0.250+ze*0.125)
     &  +xi*(0.250+et*(-0.250+ze*0.250)-0.250*ze)
      shp(1, 4)=-0.125+et*(0.250+et*(-0.125+ze*0.125)
     &  +ze*(-0.375+ze*0.125))+ze*(0.250-ze*0.125)
     &  +xi*(0.250+et*(-0.250+ze*0.250)-0.250*ze)
      shp(1, 5)=-0.125+et*(-0.250+et*(-0.125-ze*0.125)
     &  +ze*(-0.375-ze*0.125))+ze*(-0.250-ze*0.125)
     &  +xi*(0.250+et*(0.250+ze*0.250)+ze*0.250)
      shp(1, 6)=-0.375+et*(-0.250+et*(0.125+ze*0.125)
     &  +ze*(-0.125+ze*0.125))+ze*(-0.250+ze*0.125)
     &  +xi*(0.250+et*(0.250+ze*0.250)+ze*0.250)
      shp(1, 7)=-0.375+et*(-0.250+et*(0.125-ze*0.125)
     &  +ze*(0.125+ze*0.125))+ze*(0.250+ze*0.125)
     &  +xi*(0.250+et*(0.250-0.250*ze)-0.250*ze)
      shp(1, 8)=-0.125+et*(-0.250+et*(-0.125+ze*0.125)
     &  +ze*(0.375-ze*0.125))+ze*(0.250-ze*0.125)
     &  +xi*(0.250+et*(0.250-0.250*ze)-0.250*ze)
      shp(1, 9)= 0.500+et*(-0.500-0.500*ze)+ze*0.500
     &  +xi*(-0.500+et*(0.500+ze*0.500)-0.500*ze)
      shp(1,10)= 0.250+et*(-0.250+ze*(ze*0.250))+ze*(-ze*0.250)
      shp(1,11)= 0.500+et*(-0.500+ze*0.500)-0.500*ze
     &  +xi*(-0.500+et*(0.500-0.500*ze)+ze*0.500)
      shp(1,12)=-0.250+et*(0.250+ze*(-ze*0.250))+ze*(ze*0.250)
      shp(1,13)= 0.500+et*(0.500+ze*0.500)+ze*0.500
     &  +xi*(-0.500+et*(-0.500-0.500*ze)-0.500*ze)
      shp(1,14)= 0.250+et*(0.250+ze*(-ze*0.250))+ze*(-ze*0.250)
      shp(1,15)= 0.500+et*(0.500-0.500*ze)-0.500*ze
     &  +xi*(-0.500+et*(-0.500+ze*0.500)+ze*0.500)
      shp(1,16)=-0.250+et*(-0.250+ze*(ze*0.250))+ze*(ze*0.250)
      shp(1,17)=-0.250+et*(et*(0.250+ze*0.250))-0.250*ze
      shp(1,18)= 0.250+et*(et*(-0.250-ze*0.250))+ze*0.250
      shp(1,19)= 0.250+et*(et*(-0.250+ze*0.250))-0.250*ze
      shp(1,20)=-0.250+et*(et*(0.250-ze*0.250))+ze*0.250
!
!     local derivatives of the shape functions: eta-derivative
!
      shp(2, 1)=et*(0.500
     &  +ze*0.500)+ze*(-0.250-ze*0.250)+xi*(0.250+et*(-0.250
     &  -0.250*ze)+ze*(0.375+ze*0.125)+xi*(-0.125-0.125*ze))
      shp(2, 2)=xi*(0.250+et*(0.250+ze*0.250)+ze*(0.125-ze*0.125)
     &  +xi*(-0.125-0.125*ze))
      shp(2, 3)=xi*(0.250+et*(0.250-0.250*ze)+ze*(-0.125-ze*0.125)
     &  +xi*(-0.125+ze*0.125))
      shp(2, 4)=et*(0.500-0.500*ze)+ze*(0.250-ze*0.250)
     &  +xi*(0.250+et*(-0.250+ze*0.250)+ze*(-0.375+ze*0.125)
     &  +xi*(-0.125+ze*0.125))
      shp(2, 5)=et*(0.500+ze*0.500)+ze*(0.250+ze*0.250)
     &  +xi*(-0.250+et*(-0.250-0.250*ze)+ze*(-0.375-ze*0.125)
     &  +xi*(0.125+ze*0.125))
      shp(2, 6)=xi*(-0.250+et*(0.250+ze*0.250)+ze*(-0.125+ze*0.125)
     &  +xi*(0.125+ze*0.125))
      shp(2, 7)=xi*(-0.250+et*(0.250-0.250*ze)+ze*(0.125+ze*0.125)
     &  +xi*(0.125-0.125*ze))
      shp(2, 8)=et*(0.500-0.500*ze)+ze*(-0.250+ze*0.250)
     &  +xi*(-0.250+et*(-0.250+ze*0.250)+ze*(0.375-ze*0.125)
     &  +xi*(0.125-0.125*ze))
      shp(2, 9)=xi*(-0.500-0.500*ze+xi*(0.250+ze*0.250))
      shp(2,10)=xi*(-0.250+ze*(ze*0.250))
      shp(2,11)=xi*(-0.500+ze*0.500+xi*(0.250-0.250*ze))
      shp(2,12)=-0.500+ze*(ze*0.500)+xi*(0.250+ze*(-ze*0.250))
      shp(2,13)=xi*(0.500+ze*0.500+xi*(-0.250-0.250*ze))
      shp(2,14)=xi*(0.250+ze*(-ze*0.250))
      shp(2,15)=xi*(0.500-0.500*ze+xi*(-0.250+ze*0.250))
      shp(2,16)= 0.500+ze*(-ze*0.500)+xi*(-0.250+ze*(ze*0.250))
      shp(2,17)=et*(-1.000+ze*(-1.000))+xi*(et*(0.500+ze*0.500))
      shp(2,18)=xi*(et*(-0.500-0.500*ze))
      shp(2,19)=xi*(et*(-0.500+ze*0.500))
      shp(2,20)=et*(-1.000+ze*( 1.000))+xi*(et*(0.500-0.500*ze))
!
!     local derivatives of the shape functions: zeta-derivative
!
      shp(3, 1)=et*(-0.250+et*0.250-0.500*ze)+ze*0.500
     &  +xi*(-0.250+et*(0.375-0.125*et+ze*0.250)-0.250*ze
     &  +xi*(0.125+et*(-0.125)))
      shp(3, 2)=xi*(-0.250+et*(0.125+et*0.125-0.250*ze)+ze*0.250
     &  +xi*(0.125+et*(-0.125)))
      shp(3, 3)=xi*(0.250+et*(-0.125-0.125*et-0.250*ze)+ze*0.250
     &  +xi*(-0.125+et*(0.125)))
      shp(3, 4)=et*(0.250-0.250*et-0.500*ze)+ze*0.500
     &  +xi*(0.250+et*(-0.375+et*0.125+ze*0.250)-0.250*ze
     &  +xi*(-0.125+et*(0.125)))
      shp(3, 5)=et*(0.250+et*0.250+ze*0.500)+ze*0.500
     &  +xi*(-0.250+et*(-0.375-0.125*et-0.250*ze)-0.250*ze
     &  +xi*(0.125+et*(0.125)))
      shp(3, 6)=xi*(-0.250+et*(-0.125+et*0.125+ze*0.250)+ze*0.250
     &  +xi*(0.125+et*(0.125)))
      shp(3, 7)=xi*(0.250+et*(0.125-0.125*et+ze*0.250)+ze*0.250
     &  +xi*(-0.125+et*(-0.125)))
      shp(3, 8)=et*(-0.250-0.250*et+ze*0.500)+ze*0.500
     &  +xi*(0.250+et*(0.375+et*0.125-0.250*ze)-0.250*ze
     &  +xi*(-0.125+et*(-0.125)))
      shp(3, 9)=xi*(0.500+et*(-0.500)+xi*(-0.250+et*(0.250)))
      shp(3,10)=xi*(et*(ze*0.500)-0.500*ze)
      shp(3,11)=xi*(-0.500+et*(0.500)+xi*(0.250+et*(-0.250)))
      shp(3,12)=et*(ze*( 1.000))+ze*(-1.000)+xi*(et*(-0.500*ze)+
     &  ze*0.500)
      shp(3,13)=xi*(0.500+et*(0.500)+xi*(-0.250+et*(-0.250)))
      shp(3,14)=xi*(et*(-0.500*ze)-0.500*ze)
      shp(3,15)=xi*(-0.500+et*(-0.500)+xi*(0.250+et*(0.250)))
      shp(3,16)=et*(ze*(-1.000))+ze*(-1.000)+xi*(et*(ze*0.500)+ze*0.500)
      shp(3,17)= 0.500+et*(et*(-0.500))+xi*(-0.250+et*(et*0.250))
      shp(3,18)=xi*(0.250+et*(et*(-0.250)))
      shp(3,19)=xi*(-0.250+et*(et*0.250))
      shp(3,20)=-0.500+et*(et*0.500)+xi*(0.250+et*(et*(-0.250)))
!
!     computation of the local derivative of the global coordinates
!     (xs)
!
      do i=1,3
        do j=1,3
          xs(i,j)=0.d0
          do k=1,20
            xs(i,j)=xs(i,j)+xlag(i,k)*shp(j,k)
          enddo
        enddo
      enddo
!
!     computation of the jacobian determinant
!
      xj=xs(1,1)*(xs(2,2)*xs(3,3)-xs(2,3)*xs(3,2))
     &   -xs(1,2)*(xs(2,1)*xs(3,3)-xs(2,3)*xs(3,1))
     &   +xs(1,3)*(xs(2,1)*xs(3,2)-xs(2,2)*xs(3,1))
      xj=1.d0/xj
!
!     computation of the local derivative of the global coordinates
!     (xs)
!
      do i=1,3
        do j=1,3
          xs(i,j)=0.d0
          do k=1,20
            xs(i,j)=xs(i,j)+xeul(i,k)*shp(j,k)
          enddo
        enddo
      enddo
!
      return
      end
