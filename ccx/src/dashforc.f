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
      subroutine dashforc(xl,konl,vl,imat,elcon,nelcon,
     &  elas,fn,ncmat_,ntmat_,nope,lakonl,t0l,t1l,kode,elconloc,
     &  plicon,nplicon,npmat_,vel)
!
!     calculates the force of the dashpot
!
      implicit none
!
      character*8 lakonl
!
      integer konl(20),i,j,imat,ncmat_,ntmat_,nope,
     &  kode,nelcon(2,*),nplicon(0:ntmat_,*),npmat_
!
      real*8 xl(3,20),elas(21),t0l,t1l,vl(0:3,20),plconloc(82),
     &  pl(0:3,9),xn(3),al,dd,fn(0:3,*),vel(1:3,20),
     &  elcon(0:ncmat_,ntmat_,*),elconloc(21),xk,fk,
     &  plicon(0:2*npmat_,ntmat_,*)
!
!     actual positions of the nodes belonging to the contact spring
!
      do i=1,nope
         do j=1,3
            pl(j,i)=xl(j,i)+vl(j,i)
         enddo
      enddo
!
      dd=dsqrt((pl(1,2)-pl(1,1))**2
     &     +(pl(2,2)-pl(2,1))**2
     &     +(pl(3,2)-pl(3,1))**2)
      do i=1,3
         xn(i)=(pl(i,2)-pl(i,1))/dd
      enddo
!
      al=(vel(1,2)-vel(1,1))*xn(1)+
     &   (vel(2,2)-vel(2,1))*xn(2)+
     &   (vel(3,2)-vel(3,1))*xn(3)
!     
!     interpolating the material data
!     
      call materialdata_sp(elcon,nelcon,imat,ntmat_,i,t0l,t1l,
     &     elconloc,kode,plicon,nplicon,npmat_,plconloc,ncmat_)
!     
!     calculating the spring force and the spring constant
!     
      xk=elconloc(1)
      fk=xk*al
!     
      do i=1,3
         fn(i,konl(1))=-fk*xn(i)
         fn(i,konl(2))=fk*xn(i)
      enddo
!
      return
      end
      
      
