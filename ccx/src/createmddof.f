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
      subroutine createmddof(imddof,nmddof,nrset,istartset,iendset,
     &            ialset,nactdof,ithermal,mi,imdnode,nmdnode,ikmpc,
     &            ilmpc,ipompc,nodempc,nmpc,nsectors,nksector,
     &            imdmpc,nmdmpc,imdboun,nmdboun,ikboun,nboun,
     &            nset,ntie,tieset,set,lakon,kon,ipkon,labmpc)
!
!     creating a set imddof containing the degrees of freedom
!     selected by the user for modal dynamic calculations. The
!     solution will be calculated for these dof's only in order
!     to speed up the calculation.
!
      implicit none
!
      character*8 lakon(*)
      character*20 labmpc(*)
      character*81 tieset(3,*),rightset,set(*),slavset
!
      integer imddof(*),nmddof,nrset,istartset(*),iendset(*),mi(2),
     &  ialset(*),nactdof(0:mi(2),*),node,ithermal,j,k,l,
     &  ikmpc(*),ilmpc(*),ipompc(*),nodempc(3,*),nmpc,nsectors,
     &  node1,nksector,imdnode(*),nmdnode,imdmpc(*),nmdmpc,
     &  imdboun(*),nmdboun,ikboun(*),nboun,index,indexe,islav,
     &  jface,nset,ntie,nnodelem,nope,nodef(8),nelem,nface,iright,
     &  ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),kon(*),
     &  ipkon(*),i
!
!     nodes per face for hex elements
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
!
!     nodes per face for tet elements
!
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
!
!     nodes per face for linear wedge elements
!
      data ifacew1 /1,3,2,0,
     &             4,5,6,0,
     &             1,2,5,4,
     &             2,3,6,5,
     &             4,6,3,1/
!
!     nodes per face for quadratic wedge elements
!
      data ifacew2 /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
!
      do j=istartset(nrset),iendset(nrset)
         if(ialset(j).gt.0) then
            node=ialset(j)
            do l=0,nsectors-1
               node1=node+l*nksector
               call addimd(imdnode,nmdnode,node1)
               if(ithermal.ne.2) then
                  do k=1,3
                     call addimdnodedof(node1,k,ikmpc,ilmpc,ipompc,
     &                    nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                    nactdof,mi,
     &                    imdmpc,nmdmpc,imdboun,nmdboun,ikboun,nboun)
                  enddo
               else
                  k=0
                  call addimdnodedof(node1,k,ikmpc,ilmpc,ipompc,
     &                 nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                    nactdof,mi,
     &                 imdmpc,nmdmpc,imdboun,nmdboun,ikboun,nboun)
               endif
            enddo
         else
            node=ialset(j-2)
            do
               node=node-ialset(j)
               if(node.ge.ialset(j-1)) exit
               do l=0,nsectors-1
                  node1=node+l*nksector
                  call addimd(imdnode,nmdnode,node1)
                  if(ithermal.ne.2) then
                     do k=1,3
                        call addimdnodedof(node1,k,ikmpc,ilmpc,
     &                       ipompc,nodempc,nmpc,imdnode,nmdnode,imddof,
     &                       nmddof,nactdof,mi,
     &                       imdmpc,nmdmpc,imdboun,nmdboun,ikboun,nboun)
                     enddo
                  else
                     k=0
                     call addimdnodedof(node1,k,ikmpc,ilmpc,ipompc,
     &                    nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                    nactdof,mi,
     &                    imdmpc,nmdmpc,imdboun,nmdboun,ikboun,nboun)
                  endif
               enddo
            enddo
         endif
      enddo
!
!     check whether all contact slave and master nodes were selected
!
      do i=1,ntie
!     
!     check for contact conditions
!     
         if(tieset(1,i)(81:81).eq.'C') then
            rightset=tieset(3,i)
!     
!     determining the master surface
!     
            do j=1,nset
               if(set(j).eq.rightset) exit
            enddo
            if(j.gt.nset) then
               write(*,*) '*ERROR in triangucont: master surface',
     &              rightset
               write(*,*) '       does not exist'
               stop
            endif
            iright=j
!     
            do j=istartset(iright),iendset(iright)
!     
               nelem=int(ialset(j)/10.d0)
               jface=ialset(j)-10*nelem
!     
               indexe=ipkon(nelem)
!     
               if(lakon(nelem)(4:4).eq.'2') then
                  nnodelem=8
                  nface=6
               elseif(lakon(nelem)(4:4).eq.'8') then
                  nnodelem=4
                  nface=6
               elseif(lakon(nelem)(4:5).eq.'10') then
                  nnodelem=6
                  nface=4
               elseif(lakon(nelem)(4:4).eq.'4') then
                  nnodelem=3
                  nface=4
               elseif(lakon(nelem)(4:5).eq.'15') then
                  if(jface.le.2) then
                     nnodelem=6
                  else
                     nnodelem=8
                  endif
                  nface=5
                  nope=15
               elseif(lakon(nelem)(4:4).eq.'6') then
                  if(jface.le.2) then
                     nnodelem=3
                  else
                     nnodelem=4
                  endif
                  nface=5
                  nope=6
               else
                  cycle
               endif
!     
!     determining the master nodes 
!     
               if(nface.eq.4) then
                  do k=1,nnodelem
                     nodef(k)=kon(indexe+ifacet(k,jface))
                  enddo
               elseif(nface.eq.5) then
                  if(nope.eq.6) then
                     do k=1,nnodelem
                        nodef(k)=kon(indexe+ifacew1(k,jface))
                     enddo
                  elseif(nope.eq.15) then
                     do k=1,nnodelem
                        nodef(k)=kon(indexe+ifacew2(k,jface))
                     enddo
                  endif
               elseif(nface.eq.6) then
                  do k=1,nnodelem
                     nodef(k)=kon(indexe+ifaceq(k,jface))
                  enddo
               endif
!
               do l=1,nnodelem
                  node=nodef(l)
                  call addimd(imdnode,nmdnode,node)
                  if(ithermal.ne.2) then
                     do k=1,3
                        call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                       nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                       nactdof,mi,
     &                       imdmpc,nmdmpc,imdboun,nmdboun,ikboun,nboun)
                     enddo
                  endif
               enddo
            enddo
!
!           determining the slave nodes 
!
            slavset=tieset(2,i)
!
!           determining the slave surface
!
            do j=1,nset
               if(set(j).eq.slavset) exit
            enddo
            if(j.gt.nset) then
               write(*,*) '*ERROR in triangucont: ',
     &           'slave nodal surface ',slavset
               write(*,*) '       does not exist'
               stop
            endif
            islav=j
!
            do j=istartset(islav),iendset(islav)
               if(ialset(j).gt.0) then
                  node=ialset(j)
                  call addimd(imdnode,nmdnode,node)
                  if(ithermal.ne.2) then
                     do k=1,3
                        call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                       nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                       nactdof,mi,
     &                       imdmpc,nmdmpc,imdboun,nmdboun,ikboun,nboun)
                     enddo
                  endif
               else
                  k=ialset(j-2)
                  do
                     k=k-ialset(j)
                     if(k.ge.ialset(j-1)) exit
                     node=k
                     call addimd(imdnode,nmdnode,node)
                     if(ithermal.ne.2) then
                        do k=1,3
                           call addimdnodedof(node,k,ikmpc,ilmpc,
     &                          ipompc,nodempc,nmpc,imdnode,nmdnode,
     &                          imddof,nmddof,nactdof,mi,imdmpc,nmdmpc,
     &                          imdboun,nmdboun,ikboun,nboun)
                        enddo
                     endif
                  enddo
               endif
            enddo
!
         endif
      enddo
!
!     adding nodes belonging to nonlinear MPC's
!      
      do i=1,nmpc
         if((labmpc(i)(1:20).ne.'                    ').and.
     &          (labmpc(i)(1:7).ne.'CONTACT').and.
     &          (labmpc(i)(1:6).ne.'CYCLIC').and.
     &          (labmpc(i)(1:9).ne.'SUBCYCLIC')) then
            index=ipompc(i)
            if(indexe.eq.0) cycle
            node=nodempc(1,index)
            call addimd(imdnode,nmdnode,node)
            if(ithermal.ne.2) then
               do k=1,3
                  call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                 nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                 nactdof,mi,
     &                 imdmpc,nmdmpc,imdboun,nmdboun,ikboun,nboun)
               enddo
            else
               k=0
               call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &              nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &              nactdof,mi,
     &              imdmpc,nmdmpc,imdboun,nmdboun,ikboun,nboun)
            endif
         endif
      enddo
!
!     subtracting 1 to comply with the C-convention
!
      do j=1,nmddof
         imddof(j)=imddof(j)-1
      enddo
!
c      write (*,*) 'nmddof, nmdnode',nmddof,nmdnode
!
      return
      end




