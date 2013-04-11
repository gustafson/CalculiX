!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine genran(iix,g,k)
!
!     used in the Lanczos routines (cf. netlib CD)
!     generates k random numbers between 0. and 1. from seed iix 
!     and stores them in g
!
      integer iix,k
      real*4 g(k)
!
      i1=nint(iix*1974./2546.)
      i2=nint(iix*235./2546.)
      i3=nint(iix*337./2546.)
!
!     initialisation of ranewr
!
      call iniran(i1,i2,i3)
!
!     repeatedly calling ranewr to generate k random numbers
!
      do i=1,k
         g(i)=ranewr()
      enddo
!
      return
      end
