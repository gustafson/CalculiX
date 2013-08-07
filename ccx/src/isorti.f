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
      subroutine isorti(nl,list,nk,key)
!
!     Sloan routine (Int.J.Num.Meth.Engng. 28,2651-2679(1989))
!
      implicit none
!
      integer nl,nk,i,j,t,value,list(nl),key(nk)
      do 20 i=2,nl
         t=list(i)
         value=key(t)
         do 10 j=i-1,1,-1
            if(value.ge.key(list(j))) then
               list(j+1)=t
               go to 20
            endif
            list(j+1)=list(j)
 10      continue
         list(1)=t
 20   continue
      end
