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
      real*8 function finpro(n,v,j,w,k)
!
!     function used in the Lanczos routines
!     calculates the internal product of v and w
!
      integer n,j,k
      real*8 v(n),w(n)
!
      finpro=0.d0
!
      if((j.eq.1).and.(k.eq.1)) then
         do i=1,n
            finpro=finpro+v(i)*w(i)
         enddo
      else
         write(*,*) '*ERROR in finpro'
         write(*,*) '       j.ne.1 or k.ne.1'
         stop 
      endif
!
      return
      end
            
