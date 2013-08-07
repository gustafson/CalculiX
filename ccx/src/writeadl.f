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
      subroutine writeadl(nk,nactdoh,adlt,adlv)
!
!     lumping the matrix stored in adb,aub and storing the result
!     in adl
!
      implicit none
!
      integer i,j,nk,nactdoh(0:4,*)
!
      real*8 adlt(*),adlv(*)
!
      do i=1,nk
         if(nactdoh(0,i).ne.0) then
            j=0
            write(*,*) i,j,adlt(nactdoh(0,i))
         endif
         do j=1,3
            if(nactdoh(j,i).ne.0) then
               write(*,*) i,j,adlv(nactdoh(j,i))
            endif
         enddo
      enddo
!
      return
      end
