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
      subroutine remeshsurfq(tieset,ntie,set,nset,istartset,
     &  iendset,ialset,ipkon,kon,lakon,nodface,ipoface,nface,
     &  nk,ne)
!
!     catalogueing the 8-node faces belonging to
!     quadratic elements adjacent to slave and master contact
!     surfaces 
!
      implicit none
!
      character*8 lakon(*)
      character*81 tieset(3,*),surfset,set(*)
!
      integer ipoface(*),nodface(9,*),nodes(4),nk,iaux,
     &  ne,ipkon(*),kon(*),indexe,ifaceq(9,6),index1,ialset(*),
     &  ifacew(8,5),kflag,i,j,k,l,m,n,istartset(*),iendset(*),
     &  ifree,ifreenew,ntie,ipos,ij,nface,
     &  ifour,nset,ifacet(7,4),ithree
!
!     nodes belonging to the element faces
!
      data ifaceq /4,3,2,1,11,10,9,12,21,
     &            5,6,7,8,13,14,15,16,22,
     &            1,2,6,5,9,18,13,17,23,
     &            2,3,7,6,10,19,14,18,24,
     &            3,4,8,7,11,20,15,19,25,
     &            4,1,5,8,12,17,16,20,26/
      data ifacet /1,3,2,7,6,5,11,
     &             1,2,4,5,9,8,12,
     &             2,3,4,6,10,9,13,
     &             1,4,3,8,10,7,14/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
!
      ithree=3
      ifour=4
      kflag=1
!
      ifree=1
      do i=1,6*ne-1
         nodface(9,i)=i+1
      enddo
!
!     nface is the total number of quadratic slave faces
!     (6-node or 8-node faces)
!
      nface=0
      do l=1,ntie
!     
!     check for contact conditions
!     
         if((tieset(1,l)(81:81).eq.'C').or.
     &        (tieset(1,l)(81:81).eq.'-')) then
!     
!     contact constraint (only slave surfaces are remeshed)
!     
c            do m=2,3
            do m=2,2
               surfset=tieset(m,l)
!     
!     check whether facial surface
!     
               ipos=index(surfset,' ')-1
!     
               do n=1,nset
                  if(set(n).eq.surfset) exit
               enddo
!
               if(n.le.nset) then
                  if(surfset(ipos:ipos).eq.'S') cycle
               else
                  do n=1,nset
                     if((set(n)(1:ipos-1).eq.surfset(1:ipos-1)).and.
     &                  (set(n)(ipos:ipos).eq.'T')) exit
                  enddo
               endif
!     
!     determining the quadratic contact slave faces
! 
!     the faces are catalogued by the three/four corner nodes
!     in ascending order. 
!
!     ipoface(i) points to a face for which
!     node i is the lowest end node and nodface(1,ipoface(i)),
!     nodface(2,ipoface(i)) and nodface(3,ipoface(i)) are the next 
!     lower ones. If the face is triangular nodface(3,ipoface(i))
!     is zero. 
!
!     nodface(4..7,ipoface(i)) contains the middle nodes for
!     4-node faces, nodface(3..5,ipoface(i) for 3-node faces.
!
!     nodface(8,ipoface(i)) is kept for the node number of the
!     extra node in the middle of the face
!
!     nodface(9,ipoface(i))
!     is a pointer to the next surface for which node i is the
!     lowest node; if there are no more such surfaces the pointer
!     has the value zero
!     
               do ij=istartset(n),iendset(n)
!     
                  i=int(ialset(ij)/10.d0)
                  j=ialset(ij)-10*i
                  indexe=ipkon(i)
!     
!     hexahedral element
!     
                  if((lakon(i)(4:4).eq.'2').and.
     &                 ((lakon(i)(7:7).eq.' ').or.
     &                 (lakon(i)(7:7).eq.'L').or.
     &                 (lakon(i)(7:7).eq.'B'))) then
!     
                     do k=1,4
                        nodes(k)=kon(indexe+ifaceq(k,j))
                     enddo
                     call isortii(nodes,iaux,ifour,kflag)
                     index1=ipoface(nodes(1))
                     do
!     
!     adding a surface which has not been 
!     catalogued so far
!     
                        if(index1.eq.0) then
                           ifreenew=nodface(9,ifree)
                           nodface(1,ifree)=nodes(2)
                           nodface(2,ifree)=nodes(3)
                           nodface(3,ifree)=nodes(4)
                           do k=5,8
                              nodface(k-1,ifree)=
     &                             kon(indexe+ifaceq(k,j))
                           enddo
                           nodface(9,ifree)=ipoface(nodes(1))
                           ipoface(nodes(1))=ifree
                           ifree=ifreenew
                           exit
                        endif
!     
!     updating a surface which has already
!     been catalogued
!     
                        if((nodface(1,index1).eq.nodes(2)).and.
     &                       (nodface(2,index1).eq.nodes(3)).and.
     &                       (nodface(3,index1).eq.nodes(4))) exit
                        index1=nodface(9,index1)
                     enddo
!
!                    tetrahedral element
!
                  elseif(lakon(i)(4:5).eq.'10') then
!     
                     do k=1,3
                        nodes(k)=kon(indexe+ifacet(k,j))
                     enddo
                     call isortii(nodes,iaux,ithree,kflag)
                     index1=ipoface(nodes(1))
                     do
!     
!     adding a surface which has not been 
!     catalogued so far
!     
                        if(index1.eq.0) then
                           ifreenew=nodface(9,ifree)
                           nodface(1,ifree)=nodes(2)
                           nodface(2,ifree)=nodes(3)
                           do k=4,6
                              nodface(k-1,ifree)=
     &                             kon(indexe+ifacet(k,j))
                           enddo
                           nodface(9,ifree)=ipoface(nodes(1))
                           ipoface(nodes(1))=ifree
                           ifree=ifreenew
                           exit
                        endif
!     
!     updating a surface which has already
!     been catalogued
!     
                        if((nodface(1,index1).eq.nodes(2)).and.
     &                     (nodface(2,index1).eq.nodes(3))) exit
                        index1=nodface(9,index1)
                     enddo
!
!                    wedge element, triangular face
!
                  elseif((lakon(i)(4:5).eq.'15').and.(j.le.2)) then
!     
                     do k=1,3
                        nodes(k)=kon(indexe+ifacew(k,j))
                     enddo
                     call isortii(nodes,iaux,ithree,kflag)
                     index1=ipoface(nodes(1))
                     do
!     
!     adding a surface which has not been 
!     catalogued so far
!     
                        if(index1.eq.0) then
                           ifreenew=nodface(9,ifree)
                           nodface(1,ifree)=nodes(2)
                           nodface(2,ifree)=nodes(3)
                           do k=4,6
                              nodface(k-1,ifree)=
     &                             kon(indexe+ifacew(k,j))
                           enddo
                           nodface(9,ifree)=ipoface(nodes(1))
                           ipoface(nodes(1))=ifree
                           ifree=ifreenew
                           exit
                        endif
!     
!     updating a surface which has already
!     been catalogued
!     
                        if((nodface(1,index1).eq.nodes(2)).and.
     &                     (nodface(2,index1).eq.nodes(3))) exit
                        index1=nodface(9,index1)
                     enddo
!
!                    wedge element, quadrilateral face
!
                  elseif((lakon(i)(4:5).eq.'15').and.(j.ge.3)) then
!     
                     do k=1,4
                        nodes(k)=kon(indexe+ifacew(k,j))
                     enddo
                     call isortii(nodes,iaux,ifour,kflag)
                     index1=ipoface(nodes(1))
                     do
!     
!     adding a surface which has not been 
!     catalogued so far
!     
                        if(index1.eq.0) then
                           ifreenew=nodface(9,ifree)
                           nodface(1,ifree)=nodes(2)
                           nodface(2,ifree)=nodes(3)
                           nodface(3,ifree)=nodes(4)
                           do k=5,8
                              nodface(k-1,ifree)=
     &                             kon(indexe+ifacew(k,j))
                           enddo
                           nodface(9,ifree)=ipoface(nodes(1))
                           ipoface(nodes(1))=ifree
                           ifree=ifreenew
                           exit
                        endif
!     
!     updating a surface which has already
!     been catalogued
!     
                        if((nodface(1,index1).eq.nodes(2)).and.
     &                       (nodface(2,index1).eq.nodes(3)).and.
     &                       (nodface(3,index1).eq.nodes(4))) exit
                        index1=nodface(9,index1)
                     enddo
                  endif
               enddo
            enddo
         endif
      enddo
!
      nface=ifree-1
!     
      return
      end
