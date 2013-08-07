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
      subroutine createmdelem(imdnode,nmdnode,xforc,
     &              ikmpc,ilmpc,ipompc,nodempc,nmpc,imddof,nmddof,
     &              nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &              ikboun,nboun,ilboun,ithermal,imdelem,nmdelem,
     &              iponoel,inoel,prlab,prset,nprint,lakon,set,nset,
     &              ialset,ipkon,kon,istartset,iendset,nforc,
     &              ikforc,ilforc)
!
!     stores the elements for which results are requested in at
!     least one node
!
      implicit none
!
      character*6 prlab(*)
      character*8 lakon(*)
      character*81 prset(*),noset,set(*)
!
      integer iforc,node,imdnode(*),nmdnode,ikmpc(*),
     &  ilmpc(*),ipompc(*),nodempc(3,*),nmpc,imddof(*),nmddof,
     &  mi(*),nactdof(0:mi(2),*),imdmpc(*),nmdmpc,imdboun(*),nmdboun,
     &  ikboun(*),nboun,ilboun(*),ithermal,imdelem(*),nmdelem,
     &  iponoel(*),inoel(2,*),index,id,nprint,i,j,k,l,indexe,
     &  nope,nset,nrset,ialset(*),ipkon(*),kon(*),istartset(*),
     &  iendset(*),idof,m,usercload,ikforc(*),ilforc(*),nforc
!
      real*8 xforc(*)
!
!     storing all elements to which nodes in imdnode belong
!     in imdelem
!
      do m=1,nmdnode
         node=imdnode(m)
c!
c!        check whether a user concentrated load was defined
c!        in the node
c!
c         usercload=0
c         loop:do j=1,3
c            idof=8*(node-1)+j
c            call nident(ikforc,idof,nforc,id)
c            if(id.gt.0) then
c               do
c                  if(ikforc(id).eq.idof) then
c                     iforc=ilforc(id)
c                     if((xforc(iforc).lt.1.2357111318d0).and.
c     &                    (xforc(iforc).gt.1.2357111316d0)) then
c                        usercload=1
c                        exit loop
c                     endif
c                     id=id-1
c                  else
c                     exit
c                  endif
c               enddo
c            endif
c         enddo loop
!
         index=iponoel(node)
         do
            if(index.eq.0) exit
            i=inoel(1,index)
            call addimd(imdelem,nmdelem,i)
c!
c!           if a user force is defined in a node belonging
c!           to the elements the results in all nodes must
c!           be determined in order to correctly calculate the
c!           internal force in the node (input in user subroutine
c!           cload)
c!
c            if(usercload.eq.1) then
c               indexe=ipkon(i)
cc     Bernhardi start
c               if(lakon(i)(1:5).eq.'C3D8I') then
c                  nope=11
c               elseif(lakon(i)(4:4).eq.'2') then
cc     Bernhardi end
c                  nope=20
c               elseif(lakon(i)(4:4).eq.'8') then
c                  nope=8
c               elseif(lakon(i)(4:5).eq.'10') then
c                  nope=10
c               elseif(lakon(i)(4:4).eq.'4') then
c                  nope=4
c               elseif(lakon(i)(4:5).eq.'15') then
c                  nope=15
c               elseif(lakon(i)(4:4).eq.'6') then
c                  nope=6
c               elseif(lakon(i)(1:1).eq.'E') then
c                  nope=ichar(lakon(i)(8:8))-47
c               else
c                  cycle
c               endif
c!
c               do l=1,nope
c                  node=kon(indexe+l)
c                  call addimd(imdnode,nmdnode,node)
c                  if(ithermal.ne.2) then
c                     do k=1,3
c                        call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
c     &                       nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
c     &                       nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
c     &                       ikboun,nboun,ilboun)
c                     enddo
c                  else
c                     k=0
c                     call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
c     &                    nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
c     &                    nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
c     &                    ikboun,nboun,ilboun)
c                  endif
c               enddo
c            endif
!            
            index=inoel(2,index)
         enddo
      enddo
!
!     storing all nodes belonging to elements in which a user
!     load was defined in at least one node
!
      do m=1,nmdelem
         i=imdelem(m)
         indexe=ipkon(i)
c     Bernhardi start
         if(lakon(i)(1:5).eq.'C3D8I') then
            nope=11
         elseif(lakon(i)(4:4).eq.'2') then
c     Bernhardi end
            nope=20
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:5).eq.'15') then
            nope=15
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         elseif(lakon(i)(1:1).eq.'E') then
            nope=ichar(lakon(i)(8:8))-47
         else
            cycle
         endif
!
         loop: do l=1,nope
            node=kon(indexe+l)
!
!        check whether a user concentrated load was defined
!        in the node
!
            usercload=0
            do j=1,3
               idof=8*(node-1)+j
               call nident(ikforc,idof,nforc,id)
c               if(id.gt.0) then
                  do
                     if((id.eq.0).or.(ikforc(id).ne.idof)) exit
c                     if(ikforc(id).eq.idof) then
                        iforc=ilforc(id)
                        if((xforc(iforc).lt.1.2357111318d0).and.
     &                       (xforc(iforc).gt.1.2357111316d0)) then
                           usercload=1
                           exit loop
                        endif
                        id=id-1
c                     else
c                        exit
c                     endif
                  enddo
c               endif
            enddo
         enddo loop
!
!        if a user force is defined in a node belonging
!        to the elements the results in all nodes must
!        be determined in order to correctly calculate the
!        internal force in the node (input in user subroutine
!        cload)
!     
         if(usercload.eq.1) then
            do l=1,nope
               node=kon(indexe+l)
               call addimd(imdnode,nmdnode,node)
               if(ithermal.ne.2) then
                  do k=1,3
                     call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                    nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                    nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &                    ikboun,nboun,ilboun)
                  enddo
               else
                  k=0
                  call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                 nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                 nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &                 ikboun,nboun,ilboun)
               endif
            enddo
         endif
      enddo
!
!     storing the elements for which *EL PRINT was selected
!
      do m=1,nprint
         if((prlab(m)(1:4).eq.'S   ').or.
     &        (prlab(m)(1:4).eq.'E   ').or.
     &        (prlab(m)(1:4).eq.'PEEQ').or.
     &        (prlab(m)(1:4).eq.'ENER').or.
     &        (prlab(m)(1:4).eq.'SDV ').or.
     &        (prlab(m)(1:4).eq.'ELSE').or.
     &        (prlab(m)(1:4).eq.'ELKE').or.
     &        (prlab(m)(1:4).eq.'EVOL').or.
     &        (prlab(m)(1:4).eq.'HFL ')) then
            noset=prset(m)
            nrset=0
            do k=1,nset
               if(set(k).eq.noset) then
                  nrset=k
                  exit
               endif
            enddo
!
!           adding the elements belonging to nrset
!
            do j=istartset(nrset),iendset(nrset)
               if(ialset(j).gt.0) then
                  i=ialset(j)
                  call addimd(imdelem,nmdelem,i)
!
!                 in order to calculate results at the integration
!                 point of an element the results must have been
!                 determined at the nodes of this element
!
                  indexe=ipkon(i)
c     Bernhardi start
                  if(lakon(i)(1:5).eq.'C3D8I') then
                     nope=11
                  elseif(lakon(i)(4:4).eq.'2') then
c     Bernhardi end
                     nope=20
                  elseif(lakon(i)(4:4).eq.'8') then
                     nope=8
                  elseif(lakon(i)(4:5).eq.'10') then
                     nope=10
                  elseif(lakon(i)(4:4).eq.'4') then
                     nope=4
                  elseif(lakon(i)(4:5).eq.'15') then
                     nope=15
                  elseif(lakon(i)(4:4).eq.'6') then
                     nope=6
                  elseif(lakon(i)(1:1).eq.'E') then
                     nope=ichar(lakon(i)(8:8))-47
                  else
                     cycle
                  endif
!     
                  do l=1,nope
                     node=kon(indexe+l)
                     call addimd(imdnode,nmdnode,node)
                     if(ithermal.ne.2) then
                        do k=1,3
                           call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                       nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                       nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &                       ikboun,nboun,ilboun)
                        enddo
                     else
                        k=0
                        call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                       nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                       nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &                       ikboun,nboun,ilboun)
                     endif
                  enddo
               else
                  i=ialset(j-2)
                  do
                     i=i-ialset(j)
                     if(i.ge.ialset(j-1)) exit
                     call addimd(imdelem,nmdelem,i)
!
!                 in order to calculate results at the integration
!                 point of an element the results must have been
!                 determined at the nodes of this element
!
                     indexe=ipkon(i)
c     Bernhardi start
                     if(lakon(i)(1:5).eq.'C3D8I') then
                        nope=11
                     elseif(lakon(i)(4:4).eq.'2') then
c     Bernhardi end
                        nope=20
                     elseif(lakon(i)(4:4).eq.'8') then
                        nope=8
                     elseif(lakon(i)(4:5).eq.'10') then
                        nope=10
                     elseif(lakon(i)(4:4).eq.'4') then
                        nope=4
                     elseif(lakon(i)(4:5).eq.'15') then
                        nope=15
                     elseif(lakon(i)(4:4).eq.'6') then
                        nope=6
                     elseif(lakon(i)(1:1).eq.'E') then
                        nope=ichar(lakon(i)(8:8))-47
                     else
                        cycle
                     endif
!     
                     do l=1,nope
                        node=kon(indexe+l)
                        call addimd(imdnode,nmdnode,node)
                        if(ithermal.ne.2) then
                           do k=1,3
                              call addimdnodedof(node,k,ikmpc,ilmpc,
     &                             ipompc,nodempc,nmpc,imdnode,nmdnode,
     &                             imddof,nmddof,nactdof,mi,imdmpc,
     &                             nmdmpc,imdboun,nmdboun,
     &                             ikboun,nboun,ilboun)
                           enddo
                        else
                           k=0
                           call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                          nodempc,nmpc,imdnode,nmdnode,imddof,
     &                          nmddof,nactdof,mi,imdmpc,nmdmpc,imdboun,
     &                          nmdboun,ikboun,
     &                          nboun,ilboun)
                        endif
                     enddo
                  enddo
               endif
            enddo
         endif
      enddo
!
      return
      end

