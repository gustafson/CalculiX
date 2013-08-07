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
c     Bernhardi start
      subroutine genmodes(i,kon,ipkon,lakon,ne,nk,nk_,co)
!
!     generate nodes for incompatible modes
!
      implicit none
!
      character*8 lakon(*)
c
      real*8 co(3,*)
c
      integer i,kon(*),ipkon(*),ne,nope,nopeexp,
     &  nk,nk_,j,indexe,k,nodeb(8,3)
c
      indexe=ipkon(i)
c
      if(lakon(i)(1:5).eq.'C3D8I')then
         nope=8
         nopeexp=3
      else
         write(6,*) "error wrong element type in genmodes, element=",
     &               lakon(i)
      endif
!
!     generating additional nodes for the incompatible element. 
!            
      do j=1,nopeexp
         nk=nk+1
           if(nk.gt.nk_) then
              write(*,*) '*ERROR in genmodes: increase nk_'
              stop
           endif
         kon(indexe+nope+j)=nk
         do k=1,3
            co(k,nk)=0.0d0
         enddo
      enddo
!
      return
      end
c     Bernhardi end

