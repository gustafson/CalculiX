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
      subroutine calcvel(ne,nactdoh,vel,b,neq)
!
!     stores the velocities into field vel
!
      implicit none
!
      integer ne,neq,nactdoh(*),i,j
!
      real*8 vel(0:5,*),b(neq,3)
!
      do i=1,ne
         if(nactdoh(i).eq.0) cycle
         do j=1,3
            vel(j,i)=b(nactdoh(i),j)
         enddo
      enddo
!     
      return
      end
