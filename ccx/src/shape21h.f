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
      subroutine shape21h(xi,et,ze,xl,xsj,shp,iflag)
!
!     shape functions and derivatives for a 21-node quadratic
!     isoparametric brick element. -1<=xi,et,ze<=1
!
!     This is a 27-node Lagrange element in which the central
!     node and all nodes in the middle of the faces except in
!     face one were condensed
!
!     iflag=1: calculate only the value of the shape functions
!     iflag=2: calculate the value of the shape functions and
!              the Jacobian determinant
!     iflag=3: calculate the value of the shape functions, the
!              value of their derivatives w.r.t. the global
!              coordinates and the Jacobian determinant
!
      implicit none
!
      integer i,j,k,iflag
!
      real*8 shp(4,21),xs(3,3),xsi(3,3),xl(3,20),shpe(4,21),dd,
     &  dd1,dd2,dd3,shp27(4,27)
!
      real*8 xi,et,ze,xsj,fxi1,fxi2,fxi3,fet1,fet2,fet3,
     &  fze1,fze2,fze3,dfxi1,dfxi2,dfxi3,dfet1,dfet2,dfet3,
     &  dfze1,dfze2,dfze3
!
!     shape functions in one dimension
!
      fxi1=xi*(xi-1.d0)/2.d0
      fxi2=(1.d0-xi)*(1.d0+xi)
      fxi3=xi*(xi+1.d0)/2.d0
!
      fet1=et*(et-1.d0)/2.d0
      fet2=(1.d0-et)*(1.d0+et)
      fet3=et*(et+1.d0)/2.d0
!
      fze1=ze*(ze-1.d0)/2.d0
      fze2=(1.d0-ze)*(1.d0+ze)
      fze3=ze*(ze+1.d0)/2.d0
!
!     shape functions of the 27-node element
!
      shp27(4,1)=fxi1*fet1*fze1
      shp27(4,2)=fxi3*fet1*fze1
      shp27(4,3)=fxi3*fet3*fze1
      shp27(4,4)=fxi1*fet3*fze1
      shp27(4,5)=fxi1*fet1*fze3
      shp27(4,6)=fxi3*fet1*fze3
      shp27(4,7)=fxi3*fet3*fze3
      shp27(4,8)=fxi1*fet3*fze3
      shp27(4,9)=fxi2*fet1*fze1
      shp27(4,10)=fxi3*fet2*fze1
      shp27(4,11)=fxi2*fet3*fze1
      shp27(4,12)=fxi1*fet2*fze1
      shp27(4,13)=fxi2*fet1*fze3
      shp27(4,14)=fxi3*fet2*fze3
      shp27(4,15)=fxi2*fet3*fze3
      shp27(4,16)=fxi1*fet2*fze3
      shp27(4,17)=fxi1*fet1*fze2
      shp27(4,18)=fxi3*fet1*fze2
      shp27(4,19)=fxi3*fet3*fze2
      shp27(4,20)=fxi1*fet3*fze2
      shp27(4,21)=fxi2*fet2*fze2
      shp27(4,22)=fxi2*fet2*fze1
      shp27(4,23)=fxi2*fet2*fze3
      shp27(4,24)=fxi2*fet1*fze2
      shp27(4,25)=fxi3*fet2*fze2
      shp27(4,26)=fxi2*fet3*fze2
      shp27(4,27)=fxi1*fet2*fze2
!
!     shape functions of the 21-node element
!
      shp(4,1)=shp27(4,1)-(shp27(4,21)+shp27(4,24)+shp27(4,28))/4.d0
      shp(4,2)=shp27(4,2)-(shp27(4,21)+shp27(4,24)+shp27(4,25))/4.d0
      shp(4,3)=shp27(4,3)-(shp27(4,21)+shp27(4,25)+shp27(4,26))/4.d0
      shp(4,4)=shp27(4,4)-(shp27(4,21)+shp27(4,26)+shp27(4,27))/4.d0
      shp(4,5)=shp27(4,5)-(shp27(4,21)+shp27(4,23)+shp27(4,24)+
     &         shp27(4,27))/4.d0
      shp(4,6)=shp27(4,6)-(shp27(4,21)+shp27(4,23)+shp27(4,24)+
     &         shp27(4,25))/4.d0
      shp(4,7)=shp27(4,7)-(shp27(4,21)+shp27(4,23)+shp27(4,25)+
     &         shp27(4,26))/4.d0
      shp(4,8)=shp27(4,8)-(shp27(4,21)+shp27(4,23)+shp27(4,26)+
     &         shp27(4,27))/4.d0
      shp(4,9)=shp27(4,9)+shp27(4,21)/2.d0+shp27(4,24)/2.d0
      shp(4,10)=shp27(4,10)+shp27(4,21)/2.d0+shp27(4,25)/2.d0
      shp(4,11)=shp27(4,11)+shp27(4,21)/2.d0+shp27(4,26)/2.d0
      shp(4,12)=shp27(4,12)+shp27(4,21)/2.d0+shp27(4,27)/2.d0
      shp(4,13)=shp27(4,13)+shp27(4,21)/2.d0+(shp27(4,23)+
     &          shp27(4,24))/2.d0
      shp(4,14)=shp27(4,14)+shp27(4,21)/2.d0+(shp27(4,23)+
     &          shp27(4,25))/2.d0
      shp(4,15)=shp27(4,15)+shp27(4,21)/2.d0+(shp27(4,23)+
     &          shp27(4,26))/2.d0
      shp(4,16)=shp27(4,16)+shp27(4,21)/2.d0+(shp27(4,23)+
     &          shp27(4,27))/2.d0
      shp(4,17)=shp27(4,17)+shp27(4,21)/2.d0+(shp27(4,24)+
     &          shp27(4,27))/2.d0
      shp(4,18)=shp27(4,18)+shp27(4,21)/2.d0+(shp27(4,24)+
     &          shp27(4,25))/2.d0
      shp(4,19)=shp27(4,19)+shp27(4,21)/2.d0+(shp27(4,25)+
     &          shp27(4,26))/2.d0
      shp(4,20)=shp27(4,20)+shp27(4,21)/2.d0+(shp27(4,26)+
     &          shp27(4,27))/2.d0
      shp(4,21)=shp27(4,22)
!
      if(iflag.eq.1) return
!
!     derivative of the shape functions in one dimension
!
      dfxi1=(2.d0*xi-1.d0)/2.d0
      dfxi2=-2.d0*xi
      dfxi3=(2.d0*xi+1.d0)/2.d0
!
      dfet1=(2.d0*et-1.d0)/2.d0
      dfet2=-2.d0*et
      dfet3=(2.d0*et+1.d0)/2.d0
!
      dfze1=(2.d0*ze-1.d0)/2.d0
      dfze2=-2.d0*ze
      dfze3=(2.d0*ze+1.d0)/2.d0
!
!     local derivatives of the shape functions of the 27-node element: 
!     xi-derivative
!
      shp27(1,1)=dfxi1*fet1*fze1
      shp27(1,2)=dfxi3*fet1*fze1
      shp27(1,3)=dfxi3*fet3*fze1
      shp27(1,4)=dfxi1*fet3*fze1
      shp27(1,5)=dfxi1*fet1*fze3
      shp27(1,6)=dfxi3*fet1*fze3
      shp27(1,7)=dfxi3*fet3*fze3
      shp27(1,8)=dfxi1*fet3*fze3
      shp27(1,9)=dfxi2*fet1*fze1
      shp27(1,10)=dfxi3*fet2*fze1
      shp27(1,11)=dfxi2*fet3*fze1
      shp27(1,12)=dfxi1*fet2*fze1
      shp27(1,13)=dfxi2*fet1*fze3
      shp27(1,14)=dfxi3*fet2*fze3
      shp27(1,15)=dfxi2*fet3*fze3
      shp27(1,16)=dfxi1*fet2*fze3
      shp27(1,17)=dfxi1*fet1*fze2
      shp27(1,18)=dfxi3*fet1*fze2
      shp27(1,19)=dfxi3*fet3*fze2
      shp27(1,20)=dfxi1*fet3*fze2
      shp27(1,21)=dfxi2*fet2*fze2
      shp27(1,22)=dfxi2*fet2*fze1
      shp27(1,23)=dfxi2*fet2*fze3
      shp27(1,24)=dfxi2*fet1*fze2
      shp27(1,25)=dfxi3*fet2*fze2
      shp27(1,26)=dfxi2*fet3*fze2
      shp27(1,27)=dfxi1*fet2*fze2
!
!     local derivatives of the shape functions of the 27-node element: 
!     eta-derivative
!
      shp27(2,1)=fxi1*dfet1*fze1
      shp27(2,2)=fxi3*dfet1*fze1
      shp27(2,3)=fxi3*dfet3*fze1
      shp27(2,4)=fxi1*dfet3*fze1
      shp27(2,5)=fxi1*dfet1*fze3
      shp27(2,6)=fxi3*dfet1*fze3
      shp27(2,7)=fxi3*dfet3*fze3
      shp27(2,8)=fxi1*dfet3*fze3
      shp27(2,9)=fxi2*dfet1*fze1
      shp27(2,10)=fxi3*dfet2*fze1
      shp27(2,11)=fxi2*dfet3*fze1
      shp27(2,12)=fxi1*dfet2*fze1
      shp27(2,13)=fxi2*dfet1*fze3
      shp27(2,14)=fxi3*dfet2*fze3
      shp27(2,15)=fxi2*dfet3*fze3
      shp27(2,16)=fxi1*dfet2*fze3
      shp27(2,17)=fxi1*dfet1*fze2
      shp27(2,18)=fxi3*dfet1*fze2
      shp27(2,19)=fxi3*dfet3*fze2
      shp27(2,20)=fxi1*dfet3*fze2
      shp27(2,21)=fxi2*dfet2*fze2
      shp27(2,22)=fxi2*dfet2*fze1
      shp27(2,23)=fxi2*dfet2*fze3
      shp27(2,24)=fxi2*dfet1*fze2
      shp27(2,25)=fxi3*dfet2*fze2
      shp27(2,26)=fxi2*dfet3*fze2
      shp27(2,27)=fxi1*dfet2*fze2
!
!     local derivatives of the shape functions of the 27-node element: 
!     zeta-derivative
!
      shp27(3,1)=fxi1*fet1*dfze1
      shp27(3,2)=fxi3*fet1*dfze1
      shp27(3,3)=fxi3*fet3*dfze1
      shp27(3,4)=fxi1*fet3*dfze1
      shp27(3,5)=fxi1*fet1*dfze3
      shp27(3,6)=fxi3*fet1*dfze3
      shp27(3,7)=fxi3*fet3*dfze3
      shp27(3,8)=fxi1*fet3*dfze3
      shp27(3,9)=fxi2*fet1*dfze1
      shp27(3,10)=fxi3*fet2*dfze1
      shp27(3,11)=fxi2*fet3*dfze1
      shp27(3,12)=fxi1*fet2*dfze1
      shp27(3,13)=fxi2*fet1*dfze3
      shp27(3,14)=fxi3*fet2*dfze3
      shp27(3,15)=fxi2*fet3*dfze3
      shp27(3,16)=fxi1*fet2*dfze3
      shp27(3,17)=fxi1*fet1*dfze2
      shp27(3,18)=fxi3*fet1*dfze2
      shp27(3,19)=fxi3*fet3*dfze2
      shp27(3,20)=fxi1*fet3*dfze2
      shp27(3,21)=fxi2*fet2*dfze2
      shp27(3,22)=fxi2*fet2*dfze1
      shp27(3,23)=fxi2*fet2*dfze3
      shp27(3,24)=fxi2*fet1*dfze2
      shp27(3,25)=fxi3*fet2*dfze2
      shp27(3,26)=fxi2*fet3*dfze2
      shp27(3,27)=fxi1*fet2*dfze2
!
!     local derivatives of the shape functions of the 21-node element: 
!     xi-derivative
!
      shpe(1,1)=shp27(1,1)-(shp27(1,21)+shp27(1,24)+shp27(1,28))/4.d0
      shpe(1,2)=shp27(1,2)-(shp27(1,21)+shp27(1,24)+shp27(1,25))/4.d0
      shpe(1,3)=shp27(1,3)-(shp27(1,21)+shp27(1,25)+shp27(1,26))/4.d0
      shpe(1,4)=shp27(1,4)-(shp27(1,21)+shp27(1,26)+shp27(1,27))/4.d0
      shpe(1,5)=shp27(1,5)-(shp27(1,21)+shp27(1,23)+shp27(1,24)+
     &         shp27(1,27))/4.d0
      shpe(1,6)=shp27(1,6)-(shp27(1,21)+shp27(1,23)+shp27(1,24)+
     &         shp27(1,25))/4.d0
      shpe(1,7)=shp27(1,7)-(shp27(1,21)+shp27(1,23)+shp27(1,25)+
     &         shp27(1,26))/4.d0
      shpe(1,8)=shp27(1,8)-(shp27(1,21)+shp27(1,23)+shp27(1,26)+
     &         shp27(1,27))/4.d0
      shpe(1,9)=shp27(1,9)+shp27(1,21)/2.d0+shp27(1,24)/2.d0
      shpe(1,10)=shp27(1,10)+shp27(1,21)/2.d0+shp27(1,25)/2.d0
      shpe(1,11)=shp27(1,11)+shp27(1,21)/2.d0+shp27(1,26)/2.d0
      shpe(1,12)=shp27(1,12)+shp27(1,21)/2.d0+shp27(1,27)/2.d0
      shpe(1,13)=shp27(1,13)+shp27(1,21)/2.d0+(shp27(1,23)+
     &          shp27(1,24))/2.d0
      shpe(1,14)=shp27(1,14)+shp27(1,21)/2.d0+(shp27(1,23)+
     &          shp27(1,25))/2.d0
      shpe(1,15)=shp27(1,15)+shp27(1,21)/2.d0+(shp27(1,23)+
     &          shp27(1,26))/2.d0
      shpe(1,16)=shp27(1,16)+shp27(1,21)/2.d0+(shp27(1,23)+
     &          shp27(1,27))/2.d0
      shpe(1,17)=shp27(1,17)+shp27(1,21)/2.d0+(shp27(1,24)+
     &          shp27(1,27))/2.d0
      shpe(1,18)=shp27(1,18)+shp27(1,21)/2.d0+(shp27(1,24)+
     &          shp27(1,25))/2.d0
      shpe(1,19)=shp27(1,19)+shp27(1,21)/2.d0+(shp27(1,25)+
     &          shp27(1,26))/2.d0
      shpe(1,20)=shp27(1,20)+shp27(1,21)/2.d0+(shp27(1,26)+
     &          shp27(1,27))/2.d0
      shpe(1,21)=shp27(1,22)
!
!     local derivatives of the shape functions of the 21-node element: 
!     et-derivative
!
      shpe(2,1)=shp27(2,1)-(shp27(2,21)+shp27(2,24)+shp27(2,28))/4.d0
      shpe(2,2)=shp27(2,2)-(shp27(2,21)+shp27(2,24)+shp27(2,25))/4.d0
      shpe(2,3)=shp27(2,3)-(shp27(2,21)+shp27(2,25)+shp27(2,26))/4.d0
      shpe(2,4)=shp27(2,4)-(shp27(2,21)+shp27(2,26)+shp27(2,27))/4.d0
      shpe(2,5)=shp27(2,5)-(shp27(2,21)+shp27(2,23)+shp27(2,24)+
     &         shp27(2,27))/4.d0
      shpe(2,6)=shp27(2,6)-(shp27(2,21)+shp27(2,23)+shp27(2,24)+
     &         shp27(2,25))/4.d0
      shpe(2,7)=shp27(2,7)-(shp27(2,21)+shp27(2,23)+shp27(2,25)+
     &         shp27(2,26))/4.d0
      shpe(2,8)=shp27(2,8)-(shp27(2,21)+shp27(2,23)+shp27(2,26)+
     &         shp27(2,27))/4.d0
      shpe(2,9)=shp27(2,9)+shp27(2,21)/2.d0+shp27(2,24)/2.d0
      shpe(2,10)=shp27(2,10)+shp27(2,21)/2.d0+shp27(2,25)/2.d0
      shpe(2,11)=shp27(2,11)+shp27(2,21)/2.d0+shp27(2,26)/2.d0
      shpe(2,12)=shp27(2,12)+shp27(2,21)/2.d0+shp27(2,27)/2.d0
      shpe(2,13)=shp27(2,13)+shp27(2,21)/2.d0+(shp27(2,23)+
     &          shp27(2,24))/2.d0
      shpe(2,14)=shp27(2,14)+shp27(2,21)/2.d0+(shp27(2,23)+
     &          shp27(2,25))/2.d0
      shpe(2,15)=shp27(2,15)+shp27(2,21)/2.d0+(shp27(2,23)+
     &          shp27(2,26))/2.d0
      shpe(2,16)=shp27(2,16)+shp27(2,21)/2.d0+(shp27(2,23)+
     &          shp27(2,27))/2.d0
      shpe(2,17)=shp27(2,17)+shp27(2,21)/2.d0+(shp27(2,24)+
     &          shp27(2,27))/2.d0
      shpe(2,18)=shp27(2,18)+shp27(2,21)/2.d0+(shp27(2,24)+
     &          shp27(2,25))/2.d0
      shpe(2,19)=shp27(2,19)+shp27(2,21)/2.d0+(shp27(2,25)+
     &          shp27(2,26))/2.d0
      shpe(2,20)=shp27(2,20)+shp27(2,21)/2.d0+(shp27(2,26)+
     &          shp27(2,27))/2.d0
      shpe(2,21)=shp27(2,22)
!
!     local derivatives of the shape functions of the 21-node element: 
!     ze-derivative
!
      shpe(3,1)=shp27(3,1)-(shp27(3,21)+shp27(3,24)+shp27(3,28))/4.d0
      shpe(3,2)=shp27(3,2)-(shp27(3,21)+shp27(3,24)+shp27(3,25))/4.d0
      shpe(3,3)=shp27(3,3)-(shp27(3,21)+shp27(3,25)+shp27(3,26))/4.d0
      shpe(3,4)=shp27(3,4)-(shp27(3,21)+shp27(3,26)+shp27(3,27))/4.d0
      shpe(3,5)=shp27(3,5)-(shp27(3,21)+shp27(3,23)+shp27(3,24)+
     &         shp27(3,27))/4.d0
      shpe(3,6)=shp27(3,6)-(shp27(3,21)+shp27(3,23)+shp27(3,24)+
     &         shp27(3,25))/4.d0
      shpe(3,7)=shp27(3,7)-(shp27(3,21)+shp27(3,23)+shp27(3,25)+
     &         shp27(3,26))/4.d0
      shpe(3,8)=shp27(3,8)-(shp27(3,21)+shp27(3,23)+shp27(3,26)+
     &         shp27(3,27))/4.d0
      shpe(3,9)=shp27(3,9)+shp27(3,21)/2.d0+shp27(3,24)/2.d0
      shpe(3,10)=shp27(3,10)+shp27(3,21)/2.d0+shp27(3,25)/2.d0
      shpe(3,11)=shp27(3,11)+shp27(3,21)/2.d0+shp27(3,26)/2.d0
      shpe(3,12)=shp27(3,12)+shp27(3,21)/2.d0+shp27(3,27)/2.d0
      shpe(3,13)=shp27(3,13)+shp27(3,21)/2.d0+(shp27(3,23)+
     &          shp27(3,24))/2.d0
      shpe(3,14)=shp27(3,14)+shp27(3,21)/2.d0+(shp27(3,23)+
     &          shp27(3,25))/2.d0
      shpe(3,15)=shp27(3,15)+shp27(3,21)/2.d0+(shp27(3,23)+
     &          shp27(3,26))/2.d0
      shpe(3,16)=shp27(3,16)+shp27(3,21)/2.d0+(shp27(3,23)+
     &          shp27(3,27))/2.d0
      shpe(3,17)=shp27(3,17)+shp27(3,21)/2.d0+(shp27(3,24)+
     &          shp27(3,27))/2.d0
      shpe(3,18)=shp27(3,18)+shp27(3,21)/2.d0+(shp27(3,24)+
     &          shp27(3,25))/2.d0
      shpe(3,19)=shp27(3,19)+shp27(3,21)/2.d0+(shp27(3,25)+
     &          shp27(3,26))/2.d0
      shpe(3,20)=shp27(3,20)+shp27(3,21)/2.d0+(shp27(3,26)+
     &          shp27(3,27))/2.d0
      shpe(3,21)=shp27(3,22)
!
!     computation of the local derivative of the global coordinates
!     (xs)
!
      do i=1,3
        do j=1,3
          xs(i,j)=0.d0
          do k=1,20
            xs(i,j)=xs(i,j)+xl(i,k)*shpe(j,k)
          enddo
        enddo
      enddo
!
!     computation of the jacobian determinant
!
      dd1=xs(2,2)*xs(3,3)-xs(2,3)*xs(3,2)
      dd2=xs(2,3)*xs(3,1)-xs(2,1)*xs(3,3)
      dd3=xs(2,1)*xs(3,2)-xs(2,2)*xs(3,1)
      xsj=xs(1,1)*dd1+xs(1,2)*dd2+xs(1,3)*dd3
!
      if(iflag.eq.2) return
!
      dd=1.d0/xsj
!
!     computation of the global derivative of the local coordinates
!     (xsi) (inversion of xs)
!
      xsi(1,1)=dd1*dd
      xsi(1,2)=(xs(1,3)*xs(3,2)-xs(1,2)*xs(3,3))*dd
      xsi(1,3)=(xs(1,2)*xs(2,3)-xs(2,2)*xs(1,3))*dd
      xsi(2,1)=dd2*dd
      xsi(2,2)=(xs(1,1)*xs(3,3)-xs(3,1)*xs(1,3))*dd
      xsi(2,3)=(xs(1,3)*xs(2,1)-xs(1,1)*xs(2,3))*dd
      xsi(3,1)=dd3*dd
      xsi(3,2)=(xs(1,2)*xs(3,1)-xs(1,1)*xs(3,2))*dd
      xsi(3,3)=(xs(1,1)*xs(2,2)-xs(2,1)*xs(1,2))*dd
!
!     computation of the global derivatives of the shape functions
!
      do k=1,21
        do j=1,3
          shp(j,k)=shpe(1,k)*xsi(1,j)+shpe(2,k)*xsi(2,j)
     &          +shpe(3,k)*xsi(3,j)
        enddo
      enddo
!
      return
      end
