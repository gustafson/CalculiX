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
      subroutine tiefaccont(lakon,ipkon,kon,ntie,tieset,nset,set,
     &  istartset,iendset,ialset,itiefac,islavsurf,islavnode,
     &  imastnode,nslavnode,nmastnode,nslavs,nmasts,ifacecount,
     &  ipe,ime,imastop,ncont,koncont,iponoels,inoels,ifreenoels)
!
!     Catalogueing the slave faces (itieface, islavsurf)
!                  the slave nodes (islavnode, nslavnode)
!                  the master nodes (imastnode, nmastnode)
!                  the opposite trangles in the triangulation
!                         (imastop)
!                  the slave faces to which the slave nodes
!                       belong
!
!     Authors: Li,Yang; Rakotonanahary, Samoela; 
!
      implicit none
!
      character*8 lakon(*)
      character*81 tieset(3,*),slavset,mastset,set(*)
!
      logical exist
!
      integer ntie,i,j,k,l,nset,istartset(*),iendset(*),ialset(*),
     &  ifaces,nelems,jfaces,ifacem,nelemm,nslavs,nmasts,
     &  jfacem,indexe,nopes,nopem,ipkon(*),kon(*),id,
     &  ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),node,
     &  itiefac(2,*),islavsurf(2,*),islavnode(*),imastnode(*),
     &  nslavnode(ntie+1),nmastnode(ntie+1),ifacecount,islav,imast,
     &  ipe(*),ime(4,*),imastop(3,*),ipos,node1,node2,index1,
     &  index1old,ifree,ncont,koncont(4,*),iponoels(*),
     &  inoels(3,*),ifreenoels,ifreenoelold
!
! nslavnode: num of slave nodes
! islavnode: all slave nodes, tie by tie, ordered within one tie constraint
! nmastnode: num of master nodes
! imastnode: all master nodes, tie by tie, ordered within one tie constraint
! islavsurf: all slave faces
! itiefac: pointer into field islavsurf
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
      ifacecount=0
      nslavs=0
      nmasts=0
      ifreenoels=0
!
!     counters for new fields islavsurf and itiefac
!       
        PRINT *, "Tiefiaccont..."
      do i=1,ntie
!
!        check for contact conditions
!
         if(tieset(1,i)(81:81).eq.'C') then
            slavset=tieset(2,i)
!
!           check whether facial slave surface; 
!
            ipos=index(slavset,' ')
            if(slavset(ipos:ipos).eq.'S') then
               nslavnode(i+1)=nslavnode(i)
               nmastnode(i+1)=nmastnode(i)
               cycle
            endif
!
            mastset=tieset(3,i)
!
!           determining the slave surface
!
            do j=1,nset
               if(set(j).eq.slavset) exit
            enddo
            if(j.gt.nset) then
               write(*,*) '*ERROR in tiefaccont: slave surface'
               write(*,*) '       does not exist'
               stop
            endif
            islav=j
            nslavnode(i)=nslavs
!
            itiefac(1,i)=ifacecount+1
            do j=istartset(islav),iendset(islav)
               if(ialset(j).gt.0) then
!
!               put all the num, made of element num and face num 
!               of slave face, into islavsurf(1,*)               
!
                  ifacecount=ifacecount+1
                  islavsurf(1,ifacecount)=ialset(j)
!               
!           Decide islavnode, and nslavnode
!
                  ifaces = ialset(j)
                  nelems = int(ifaces/10)
                  jfaces = ifaces - nelems*10
                  indexe = ipkon(nelems)
!
                  if(lakon(nelems)(4:4).eq.'2') then
                      nopes=8
                  elseif(lakon(nelems)(4:4).eq.'8') then
                      nopes=4
                  elseif(lakon(nelems)(4:5).eq.'10') then
                      nopes=6
                  elseif(lakon(nelems)(4:4).eq.'4') then
                      nopes=3
                  endif
!     
                  if(lakon(nelems)(4:4).eq.'6') then
                    if(jfaces.le.2) then
                       nopes=3
                    else
                       nopes=4
                    endif
                  endif
                  if(lakon(nelems)(4:5).eq.'15') then
                    if(jfaces.le.2) then
                       nopes=6
                    else
                       nopes=8
                    endif
                  endif   
!                  
                  do l=1,nopes
                     if((lakon(nelems)(4:4).eq.'2').or.
     &                  (lakon(nelems)(4:4).eq.'8')) then
                        node=kon(indexe+ifaceq(l,jfaces))
                     elseif((lakon(nelems)(4:4).eq.'4').or.
     &                      (lakon(nelems)(4:5).eq.'10')) then
                        node=kon(indexe+ifacet(l,jfaces))
                     elseif(lakon(nelems)(4:4).eq.'6') then
                        node=kon(indexe+ifacew1(l,jfaces))
                     elseif(lakon(nelems)(4:5).eq.'15') then
                        node=kon(indexe+ifacew2(l,jfaces))
                     endif
                     call nident(islavnode(nslavnode(i)+1),node,
     &                    nslavs-nslavnode(i),id)
                     exist=.FALSE.
                     if(id.gt.0) then
                        if(islavnode(nslavnode(i)+id).eq.node) then
                           exist=.TRUE.
                        endif
                     endif
                     if(.not.exist) then
                        nslavs=nslavs+1
                        do k=nslavs,id+2,-1
                           islavnode(k)=islavnode(k-1)
                        enddo
                        islavnode(id+1)=node
                     endif
!
!                    filling fields iponoels and inoels
!
                     ifreenoelold=iponoels(node)
                     ifreenoels=ifreenoels+1
                     iponoels(node)=ifreenoels
                     inoels(1,ifreenoels)=ifacecount
                     inoels(2,ifreenoels)=l
                     inoels(3,ifreenoels)=ifreenoelold
                  enddo
!     
               endif
            enddo
            nslavnode(ntie+1)=nslavs
            itiefac(2,i)=ifacecount
!
!           determining the master surface
!
            do j=1,nset
               if(set(j).eq.mastset) exit
            enddo
            if(j.gt.nset) then
               write(*,*) '*ERROR in tiefaccont: master surface'
               write(*,*) '       does not exist'
               stop
            endif
            imast=j
            nmastnode(i)=nmasts
!
            do j=istartset(imast),iendset(imast)
               if(ialset(j).gt.0) then
!               
!           Decide imastnode, and nmastnode
!
                  ifacem = ialset(j)
                  nelemm = int(ifacem/10)
                  jfacem = ifacem - nelemm*10
                  indexe = ipkon(nelemm)
!
                  if(lakon(nelemm)(4:4).eq.'2') then
                      nopem=8
                  elseif(lakon(nelemm)(4:4).eq.'8') then
                      nopem=4
                  elseif(lakon(nelemm)(4:5).eq.'10') then
                      nopem=6
                  elseif(lakon(nelemm)(4:4).eq.'4') then
                      nopem=3
                  endif
!     
                  if(lakon(nelemm)(4:4).eq.'6') then
                    if(jfacem.le.2) then
                       nopem=3
                    else
                       nopem=4
                    endif
                  endif
                  if(lakon(nelemm)(4:5).eq.'15') then
                    if(jfacem.le.2) then
                       nopem=6
                    else
                       nopem=8
                    endif
                  endif   
!                  
                  do l=1,nopem
                     if((lakon(nelemm)(4:4).eq.'2').or.
     &                    (lakon(nelemm)(4:4).eq.'8')) then
                        node=kon(indexe+ifaceq(l,jfacem))
                     elseif((lakon(nelemm)(4:4).eq.'4').or.
     &                       (lakon(nelemm)(4:5).eq.'10')) then
                        node=kon(indexe+ifacet(l,jfacem))
                     elseif(lakon(nelemm)(4:4).eq.'6') then
                        node=kon(indexe+ifacew1(l,jfacem))
                     elseif(lakon(nelemm)(4:5).eq.'15') then
                        node=kon(indexe+ifacew2(l,jfacem))
                     endif
                     call nident(imastnode(nmastnode(i)+1),node,
     &                    nmasts-nmastnode(i),id)
                     exist=.FALSE.
                     if(id.gt.0) then
                        if(imastnode(nmastnode(i)+id).eq.node) then
                           exist=.TRUE.
                        endif
                     endif
                     if(exist) cycle
                     nmasts=nmasts+1
                     do k=nmasts,id+2,-1
                        imastnode(k)=imastnode(k-1)
                     enddo
                     imastnode(id+1)=node
                  enddo
! 
               endif
            enddo
            nmastnode(ntie+1)=nmasts
!
         else
!
!           no contact tie
!
            nslavnode(i+1)=nslavnode(i)
            nmastnode(i+1)=nmastnode(i)
         endif 
      enddo      
!
!     catalogueing the edges in the triangulation
!     determining neighboring triangles
!
      ifree=0
      do j=1,ncont
         do k=1,3
            node1=koncont(k,j)
            if(k.eq.3) then
               node2=koncont(1,j)
            else
               node2=koncont(k+1,j)
            endif
!
            if(k.eq.1) then
               ipos=3
            else
               ipos=k-1
            endif
!
!           making sure that node1 < node2
!
            if(node1.gt.node2) then
               node=node1
               node1=node2
               node2=node
            endif    
            if(ipe(node1).eq.0) then
               ifree=ifree+1
               ipe(node1)=ifree
               ime(1,ifree)=node2
               ime(2,ifree)=j
               ime(3,ifree)=ipos
            else
               index1=ipe(node1)
               if(ime(1,index1).eq.node2) then
                  imastop(ipos,j)=ime(2,index1)
                  imastop(ime(3,index1),ime(2,index1))=j
                  ipe(node1)=ime(4,index1)
                  cycle
               endif
!
               index1old=index1
               index1=ime(4,index1)
               do
                  if(index1.eq.0) then
                     ifree=ifree+1
                     ime(4,index1old)=ifree
                     ime(1,ifree)=node2
                     ime(2,ifree)=j
                     ime(3,ifree)=ipos                       
                     exit
                  endif
                  if(ime(1,index1).eq.node2) then
                     imastop(ipos,j)=ime(2,index1)
                     imastop(ime(3,index1),ime(2,index1))=j
                     ime(4,index1old)=ime(4,index1)
                     exit
                  endif
                  index1old=index1
                  index1=ime(4,index1)
               enddo
            endif
         enddo
      enddo
!
      return
      end
