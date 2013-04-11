
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
      subroutine neartriangle(p,xn,xo,yo,zo,x,y,z,nx,ny,nz,n,neigh,
     &  kneigh,itietri,ntie,straight,imastop,itri,ctie)
!
!     check for a triangle such that a straight line
!     through p and with direction xn cuts this triangle
!     
      implicit none
!
      character*1 c
!
      integer nx(*),ny(*),nz(*),n,kneigh,neigh(*),itietri(2,ntie),k,m1,
     &  isol,itri,imastop(3,*),i,ntie,ctie,r1
!
      real*8 p(3),xn(3),xo(*),yo(*),zo(*),x(*),y(*),z(*),straight(16,*),
     &  al,al1,al2,rand,err
!     
!     determining the kneigh neighboring triangles
!     
!     
      err=0d0
      CALL time(r1)
c      CALL srand(REAL(MOD(r1,1000)))
!
      call near3d(xo,yo,zo,x,y,z,nx,ny,nz,p(1),p(2),p(3),
     &     n,neigh,kneigh)
!     
      isol=0
!     
      loop1: do k=1,kneigh
         itri=neigh(k)+itietri(1,ctie)-1
         loop2: do
            al=-(straight(16,itri)+straight(13,itri)*(p(1)+err*rand(1))
     &        +straight(14,itri)*(p(2)+err*rand(1))+
     &         straight(15,itri)*(p(3)+err*rand(1)))/
     &        (straight(13,itri)*xn(1)+straight(14,itri)*xn(2)
     &        +straight(15,itri)*xn(3))
!     
            do m1=1,3
               al1=straight(4*m1-3,itri)*(p(1)+err*rand(1))+
     &              straight(4*m1-2,itri)*(p(2)+err*rand(1))+
     &              straight(4*m1-1,itri)*(p(3)+err*rand(1))
               al2=straight(4*m1-3,itri)*xn(1)+
     &              straight(4*m1-2,itri)*xn(2)+
     &              straight(4*m1-1,itri)*xn(3)
               if(al1+al*al2+straight(4*m1,itri).gt.1.d-10)then
                  if(al.lt.1.d-10) cycle loop1
                  itri=imastop(m1,itri)
                  if(itri.eq.0) cycle loop1
                  cycle loop2
               endif
            enddo
!     
            isol=1
!     
            exit loop1
         enddo loop2
      enddo loop1
!
      if(isol.ne.1) itri=0
!     
      return
      end
      
