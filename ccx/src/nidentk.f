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
!
!     identifies the position id of px in an ordered array
!     x of integers; array x has two indices. The ordered array
!     is the first index 
!
!     id is such that x(1,id).le.px and x(1,id+1).gt.px
!
      SUBROUTINE nIDENTk(X,PX,N,ID,k)
      IMPLICIT none
      integer x,px,n,id,n2,m,k
      DIMENSION X(k,N)
      id=0
      if(n.eq.0) return
      N2=N+1
      do
         M=(N2+ID)/2
         IF(PX.GE.X(1,M)) then
            ID=M
         else
            N2=M
         endif
         IF((N2-ID).EQ.1) return
      enddo
      END
