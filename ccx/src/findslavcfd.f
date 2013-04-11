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
      subroutine findslavcfd(nmpc,labmpc,ipompc,nodempc,islav,
     &  nslav)
!
!     find the slave nodes in a cyclic symmetry constraint for
!     CFD calculations
!
      implicit none
!
      character*20 labmpc(*)
!
      integer i,j,nmpc,node,nodempc(3,*),ipompc(*),nslav,id,
     &  islav(*)
!
      nslav=0
!
!     cyclic symmetry slave nodes are the dependent nodes in a 
!     CYCLIC MPC
!
      do i=1,nmpc
         if(labmpc(i)(1:6).ne.'CYCLIC') cycle
         node=nodempc(1,ipompc(i))
         call nident(islav,node,nslav,id)
         if(id.gt.0) then
            if(islav(id).eq.node) cycle
         endif
         nslav=nslav+1
         do j=nslav,id+2,-1
            islav(j)=islav(j-1)
         enddo
         islav(id+1)=node
      enddo
!
      return
      end

