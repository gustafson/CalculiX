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
      subroutine createinum(ipkon,inum,kon,lakon,nk,ne,cflag,nelemload,
     &  nload,nodeboun,nboun,ndirboun,ithermal)
!
!     determines inum in case no extrapolation is requested in the
!     input deck (e.g. only nodal variables are requested)
!
      implicit none
!
      character*1 cflag
      character*8 lakon(*),lakonl
!
      integer ipkon(*),inum(*),kon(*),ne,indexe,nope,
     &  nk,i,j,nelemload(2,*),nload,node,nboun,
     &  nodeboun(*),ndirboun(*),ithermal(2)
!
      do i=1,nk
         inum(i)=0
      enddo
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         indexe=ipkon(i)
         lakonl=lakon(i)
!
         if(lakonl(1:1).eq.'F') then
            cycle
         elseif(lakonl(4:4).eq.'2') then
            nope=20
         elseif(lakonl(4:4).eq.'8') then
            nope=8
         elseif(lakonl(4:5).eq.'10') then
            nope=10
         elseif(lakonl(4:4).eq.'4') then
            nope=4
         elseif(lakonl(4:5).eq.'15') then
            nope=15
         elseif(lakonl(4:4).eq.'6') then
            nope=6
         elseif((lakon(i)(1:1).eq.'E').and.(lakon(i)(7:7).eq.'A'))then
            inum(kon(indexe+1))=inum(kon(indexe+1))+1
            inum(kon(indexe+2))=inum(kon(indexe+2))+1
            cycle
         else
            cycle
         endif
!
!        counting the number of elements a node belongs to
!
         do j=1,nope
            inum(kon(indexe+j))=inum(kon(indexe+j))+1
         enddo
c     Bernhardi start
c        incompatible modes elements
         if(lakonl(1:5).eq.'C3D8I') then
            do j=1,3
               inum(kon(indexe+nope+j))=inum(kon(indexe+nope+j))+1
            enddo
         endif
c     Bernhardi end
!
      enddo
c!
c!     printing values for environmental film, radiation and
c!     pressure nodes (these nodes are considered to be network
c!     nodes)
c!
c      do i=1,nload
c         node=nelemload(2,i)
c         if(node.gt.0) then
c            if(inum(node).gt.0) cycle
c            inum(node)=-1
c         endif
c      enddo
c!
c!     printing values of prescribed boundary conditions (these
c!     nodes are considered to be structural nodes)
c!
c      if(ithermal(2).gt.1) then
c         do i=1,nboun
c            node=nodeboun(i)
c            if(inum(node).ne.0) cycle
c            if((cflag.ne.' ').and.(ndirboun(i).eq.3)) cycle
c            inum(node)=1
c         enddo
c      endif
!
      return
      end
