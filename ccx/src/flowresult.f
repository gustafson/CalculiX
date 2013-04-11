!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine flowresult(ntg,itg,cam,vold,v,nload,sideload,
     &     nelemload,xloadact,nactdog,network,mi)
!     
      implicit none
!     
      character*20 sideload(*) 
!     
      integer i,j,nload,node,ntg,itg(*),nelemload(2,*),
     &     nactdog(0:3,*),network,mi(2)
!     
      real*8 cam(3),vold(0:mi(2),*),v(0:mi(2),*),xloadact(2,*)
!     
!     calculating the change of gas temperature: is taken
!     into account in the global convergence for purely
!     thermal networks (reason: for purely thermal networks
!     the network solution is not iterated to speed up
!     the calculation)
!
      if(network.eq.0) then
         do i=1,ntg
            node=itg(i)
            if(nactdog(0,node).eq.0) cycle
            if(dabs(vold(0,node)-v(0,node)).gt.cam(2)) then
               cam(2)=dabs(vold(0,node)-v(0,node))
               cam(5)=node+0.5d0
            endif
         enddo
      endif
!     
!     replacing vold by v (including the static temperature for
!     gases and the critical depth for liquid channels)
!
      do i=1,ntg
         node=itg(i)
c         do j=0,2
         do j=0,3
            vold(j,node)=v(j,node)
         enddo
      enddo
!     
!     updating the film boundary conditions
!     
      do i=1,nload
         if(sideload(i)(3:4).eq.'FC') then
            node=nelemload(2,i)
            xloadact(2,i)=vold(0,node)
         endif
      enddo
!     
!     updating the pressure boundary conditions
!     
      do i=1,nload
         if(sideload(i)(3:4).eq.'NP') then
            node=nelemload(2,i)
            xloadact(1,i)=vold(2,node)
         endif
      enddo
!      
      return
      end








