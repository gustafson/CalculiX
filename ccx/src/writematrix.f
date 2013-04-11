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
      subroutine writematrix(au,ad,irow,jq,neq,number)
!
!     writes an MPC to standard output (for debugging purposes)
!
      implicit none
!
      integer irow(*),jq(*),neq,i,j,number
      real*8 au(*),ad(*)
!
      write(*,*) 'matrix number ',number
!
      do i=1,neq
         write(*,*) 'row ',i,' value ',ad(i)
      enddo
!
      do i=1,neq
         do j=jq(i),jq(i+1)-1
            write(*,*) 'colomn ',i,' row ',irow(j),' value ',au(j)
         enddo
      enddo
!
      return
      end

