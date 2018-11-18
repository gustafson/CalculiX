!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2018 Guido Dhondt
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
!     Calculating the residual in CFD-calculations
!
      subroutine calcresvfluid(n,a,b,au,ia,ja,nflnei,x,res)
!
      implicit none
!
      integer i,j,n,ia(*),ja(*),nflnei
!
      real*8 xmax,a(*),b(*),au(*),x(*),res,resi,vel
!
      xmax=0.d0
      res=0.d0
!
!     x
! 
      do i=1,n
         resi=a(ja(i)+1)*x(ia(ja(i)+1))
         do j=ja(i)+2,ja(i+1)
            resi=resi+a(j)*x(ia(j))
         enddo
         res=res+((resi-b(i))/au(nflnei+i))**2
      enddo
!
!     y
!
      do i=1,n
         resi=a(ja(i)+1)*x(n+ia(ja(i)+1))
         do j=ja(i)+2,ja(i+1)
            resi=resi+a(j)*x(n+ia(j))
         enddo
         res=res+((resi-b(n+i))/au(nflnei+i))**2
      enddo
!
!     z
!
      do i=1,n
         resi=a(ja(i)+1)*x(2*n+ia(ja(i)+1))
         do j=ja(i)+2,ja(i+1)
            resi=resi+a(j)*x(2*n+ia(j))
         enddo
         res=res+((resi-b(2*n+i))/au(nflnei+i))**2
      enddo
!
      do i=1,n
         vel=dsqrt(x(i)**2+x(n+i)**2+x(2*n+i)**2)
         if(vel.gt.xmax) xmax=vel
      enddo
!
      res=dsqrt(res/n)/xmax
!
      return
      end
