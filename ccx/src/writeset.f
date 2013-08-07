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
      subroutine writeset(nset,set,istartset,iendset,ialset)
!
!     writes an MPC to standard output (for debugging purposes)
!
      implicit none
!
      character*81 set(*)
      integer nset,istartset(*),iendset(*),ialset(*),i,j
!
      do i=1,nset
         write(*,'(''SET '',i10,1x,a81)') i,set(i)
         do j=istartset(i),iendset(i)
            write(*,'(i10)') ialset(j)
         enddo
      enddo
!
      return
      end

