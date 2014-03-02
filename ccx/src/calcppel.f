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
      subroutine calcppel(ne,nactdoh,ppel,b)
!
!     stores the first pressure correction p' into field ppel
!
      implicit none
!
      integer ne,nactdoh(*),i
!
      real*8 ppel(*),b(*)
!
      do i=1,ne
         if(nactdoh(i).eq.0) cycle
         ppel(i)=b(nactdoh(i))
      enddo
!     
      return
      end
