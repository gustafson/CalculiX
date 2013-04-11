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
      subroutine gennactdofinv(nactdof,nactdofinv,nk,mi,nodorig,
     &  ipkon,lakon,kon,ne)
!
!     inverting field nactdof, i.e. creating field nactdofinv
!     listing the node for each independent dof. For expanded
!     2-D structures this is the original 2-D node. Field 
!     nactdofinv is used for the messages listing the node in
!     which the actual deviation or residual force is maximum
!
      implicit none
!
      character*8 lakon(*),lakonl
!
      integer mi(*),nactdof(0:mi(2),*),nactdofinv(*),nk,nodorig(*),
     &  ipkon(*),i,j,l,ne,indexe,node2d,node3d,indexe2d,
     &  node3(8,3),node6(3,6),node8(3,8),kon(*),mt
!
      data node3 /1,4,8,5,12,20,16,17,9,11,15,13,
     &            0,0,0,0,2,3,7,6,10,19,14,18/
      data node6 /1,13,4,2,14,5,3,15,6,7,0,10,8,0,11,9,0,12/
      data node8 /1,17,5,2,18,6,3,19,7,4,20,8,9,0,13,10,0,14,
     &      11,0,15,12,0,16/
!
!     initialization of nodorig (node-original)
!
      do i=1,nk
         nodorig(i)=i
      enddo
!
!     replacing the 3-D nodes by their 2-D equivalents for 1d/2d
!     structures
!
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         lakonl=lakon(i)
c         if((lakonl(7:7).eq.' ').or.(lakonl(7:7).eq.'G').or.
c     &      (lakonl(1:1).ne.'C')) cycle
         if((lakonl(7:7).eq.' ').or.(lakonl(7:7).eq.'I').or.
     &      (lakonl(1:1).ne.'C')) cycle
         indexe=ipkon(i)
!
         if(lakonl(4:5).eq.'15') then
!
!        6-noded 2D element
!
            indexe2d=indexe+15
            do j=1,6
               node2d=kon(indexe2d+j)
c               do l=1,3
c                  if(node6(l,j).ne.0) then
                     node3d=kon(indexe+node6(1,j))
                     nodorig(node3d)=node2d
                     nodorig(node3d+1)=node2d
                     nodorig(node3d+2)=node2d
c                  endif
c               enddo
            enddo
         elseif(lakonl(7:7).eq.'B') then
!
!           3-noded beam element
!
            indexe2d=indexe+20
            do j=1,3
               node2d=kon(indexe2d+j)
               do l=1,8
                  if(node3(l,j).ne.0) then
                     node3d=kon(indexe+node3(l,j))
                     nodorig(node3d)=node2d
                  endif
               enddo
            enddo
         else
!
!           8-noded 2D element
!
            indexe2d=indexe+20
            do j=1,8
               node2d=kon(indexe2d+j)
c               do l=1,3
c                  if(node8(l,j).ne.0) then
                     node3d=kon(indexe+node8(1,j))
                     nodorig(node3d)=node2d
c                     write(*,*) node3d,node2d
                     nodorig(node3d+1)=node2d
c                     write(*,*) node3d+1,node2d
                     nodorig(node3d+2)=node2d
c                     write(*,*) node3d+2,node2d
c                  endif
c               enddo
            enddo
         endif
      enddo
!
!     storing the nodes (in C convention, i.e. starting with 0)
!     in field nactdofinv
!
      mt=mi(2)+1
      do i=1,nk
         do j=0,mi(2)
            if(nactdof(j,i).eq.0) cycle
            nactdofinv(nactdof(j,i))=(nodorig(i)-1)*mt+j
         enddo
      enddo
c      do i=1,25
c         write(*,*) i-1,nactdofinv(i),nactdofinv(i)/mt+1
c      enddo
!
      return
      end
