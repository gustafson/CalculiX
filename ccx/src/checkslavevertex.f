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
      subroutine checkslavevertex(lvertex,nvertex,pvertex,
     &  itriacornerl,xl2)
!
!     check whether triangular master vertex lies within the slave
!     surface
!
      implicit none
!
      integer nvertex,lvertex(*),nodel,i,
     &  itriacornerl(*)
!
      real*8 pvertex(3,*),xl2(3,*)
!
      
      if(nvertex.ne.0) then
         nodel=lvertex(nvertex)
      else
         nodel=0
      endif
      if(nodel.ne.0) then
!
!         S-edge lvertex(nvertex) (local number, applies to
!         the nodes as well as to the edges) was cut
!
         if(itriacornerl(nodel).eq.1) then
            nvertex=nvertex+1
            do i=1,3
               pvertex(i,nvertex)=xl2(i,nodel)
            enddo
            lvertex(nvertex)=0
!
            itriacornerl(nodel)=2
         endif
      endif
!
      return
      end
