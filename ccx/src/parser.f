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
      subroutine parser(text,textpart,free,n)
      implicit none
!
!     parser for calinput
!
      character*20 textpart(16)
      character*132 text
!
      logical free
!
      integer n,i,j,k
!
      free=.false.
      n=0
      j=0
!
      do i=1,132
        if(text(i:i).ne.',') then
          if(j.eq.0) then
            if(text(i:i).eq.' ') cycle
            n=n+1
          endif
          j=j+1
          if(j.le.20) textpart(n)(j:j)=text(i:i)
        else
          do k=j+1,20
            textpart(n)(k:k)=' '
          enddo
          j=0
          free=.true.
        endif
      enddo
!
      return
      end
