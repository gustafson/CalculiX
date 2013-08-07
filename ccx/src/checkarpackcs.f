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
      subroutine checkarpackcs(iponoel,inoel,ne,ipkon,lakon,
     &  kon,iactnode,iactelem,iecovered,incovered,itime)
!
      implicit none
!
      character*8 lakon(*)
!
      integer iponoel(*),inoel(2,*),ne,ipkon(*),kon(*),inoelfree,
     &  nope,indexe,iactelem(*),iactnode(*),iecovered(*),nactive,
     &  itime(*),i,j,k,index,node,node1,id,iref,ielem,iact,il,ih,
     &  incovered(*),nei1,nei2,nei3,ineigh10(3,10),ineigh20(3,20)
!
      data ineigh10 /5,7,8,5,6,9,6,7,10,8,9,10,
     &               1,2,2,2,3,3,3,1,1,
     &               1,4,4,2,4,4,3,4,4/
      data ineigh20 /9,12,17,9,10,18,10,11,19,11,12,20,
     &             13,16,17,13,14,18,14,15,19,15,16,20,
     &             1,2,2,2,3,3,3,4,4,4,1,1,
     &             5,6,6,6,7,7,7,8,8,8,5,5,
     &             1,5,5,2,6,6,3,7,7,4,8,8/
!
!     determining the elements belonging to the nodes of
!     the elements
!
      inoelfree=1
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if(lakon(i)(1:1).eq.'F') cycle
         if(lakon(i)(4:4).eq.'2') then
            nope=20
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         else
            nope=15
         endif
         indexe=ipkon(i)
         do j=1,nope
            node=kon(indexe+j)
            inoel(1,inoelfree)=i
            inoel(2,inoelfree)=iponoel(node)
            iponoel(node)=inoelfree
            inoelfree=inoelfree+1
         enddo
      enddo
!
!     determining an active (element,node) set
!
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if(lakon(i)(1:1).eq.'F') cycle
         if(lakon(i)(4:4).eq.'2') then
            nope=20
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         else
            nope=15
         endif
         indexe=ipkon(i)
         node=kon(indexe+1)
         iactelem(1)=i
         iactnode(1)=node
         incovered(node)=1
         nactive=1
         exit
      enddo
!
!     covering all elements through neighboring relations  
!
      do
         if(nactive.eq.0) exit
         ielem=iactelem(1)
c         write(*,*) 'ielem ',ielem
c         do i=1,nactive
c            write(*,*) iactelem(i),iactnode(i)
c         enddo
         node=iactnode(1)
         iref=itime(node)
         indexe=ipkon(ielem)
!
!        removing the element from the active sets
!
         do i=1,nactive-1
            iactelem(i)=iactelem(i+1)
            iactnode(i)=iactnode(i+1)
         enddo
         iecovered(ielem)=1
         nactive=nactive-1
!
!        loop over all nodes belonging to the element
!
         loop:do
            do k=1,nope
               node1=kon(indexe+k)
               if(incovered(node1).eq.1) cycle
!     
!     checking for neighbors
!     
               if(nope.eq.20) then
                  nei1=kon(indexe+ineigh20(1,k))
                  nei2=kon(indexe+ineigh20(2,k))
                  nei3=kon(indexe+ineigh20(3,k))
               elseif(nope.eq.10) then
                  nei1=kon(indexe+ineigh10(1,k))
                  nei2=kon(indexe+ineigh10(2,k))
                  nei3=kon(indexe+ineigh10(3,k))
               else
                  write(*,*) '*ERROR in checkarpackcs: case not covered'
                  stop
               endif
               if(incovered(nei1).eq.1) then
                  iref=itime(nei1)
               elseif(incovered(nei2).eq.1) then
                  iref=itime(nei2)
               elseif(incovered(nei3).eq.1) then
                  iref=itime(nei3)
               else
                  cycle
               endif
               
               incovered(node1)=1
c     if(node1.eq.node) cycle
!     
!     checking for continuity of field time (to be done)
!     
               iact=itime(node1)
               il=iact
               ih=iact
               if(iact.le.iref) then
                  do
                     ih=ih+180
                     if(ih.ge.iref) exit
                     il=ih
                  enddo
               else
                  do
                     il=il-180
                     if(il.le.iref) exit
                     ih=il
                  enddo
               endif
               if((ih-iref)>(iref-il)) then
                  itime(node1)=il
               else
                  itime(node1)=ih
               endif
               write(*,*) 'check ',node1,iref,iact,il,ih,itime(node1)
!     
!     covering all elements belonging to node node1
!     
               index=iponoel(node1)
               do
                  ielem=inoel(1,index)
                  if(iecovered(ielem).eq.0) then
                     call nident(iactelem,ielem,nactive,id)
                     if(id.gt.0) then
                        if(iactelem(id).eq.ielem) then
!     
!     element already belongs to the active set
!     
                           index=inoel(2,index)
                           if(index.eq.0) exit
                           cycle
                        endif
                     endif
!     
!     new element to be added to the active set
!     
                     nactive=nactive+1
                     do j=nactive,id+2,-1
                        iactelem(j)=iactelem(j-1)
                        iactnode(j)=iactnode(j-1)
                     enddo
                     iactelem(id+1)=ielem
                     iactnode(id+1)=node1
                  endif
                  index=inoel(2,index)
                  if(index.eq.0) exit
               enddo
            enddo
            do k=1,nope
               node1=kon(indexe+k)
               if(incovered(node1).eq.0) cycle loop
            enddo
            exit
         enddo loop
      enddo
!
      return
      end
