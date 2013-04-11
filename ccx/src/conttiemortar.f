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
      subroutine conttiemortar(lakon,ipkon,kon,ntie,tieset,nset,set,
     &  itiefac,islavsurf,islavnode,
     &  imastnode,nslavnode,nmastnode,nslavs,
     &  iponoels,inoels,
     &  ipoface,nodface,nk,
     &  nboun,ndirboun,nodeboun,xboun,
     &  nmpc,ipompc,nodempc,coefmpc,
     &  ikboun,ilboun,ikmpc,ilmpc,
     &  nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
     &  nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
     &  islavborder)
c>     subroutine to catalogue SPC'S and MPC'S for master and slave node
c>     and to locate slave nodes belonging to the border of the slave surface     
c>     needed for mortar contact
c>
c>     author: Sitzmann,Saskia
c>
        implicit none

      logical nodeslavsurf
!
      character*8 lakon(*)
      character*81 tieset(3,*),slavset,mastset,set(*)
!
      logical debug
!
      integer ntie,i,j,k,l,nset,
     &  ifaces,nelems,jfaces,ifacem,nelemm,nslavs,
     &  jfacem,indexe,nopes,nopem,ipkon(*),kon(*),id,
     &  ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),node,
     &  itiefac(2,*),islavsurf(2,*),islavnode(*),imastnode(*),
     &  nslavnode(ntie+1),nmastnode(ntie+1),ifacecount,islav,imast,
     &  ipos,index1,iponoels(*),inoels(3,*),ifreenoelold,
     &  numbern,numberf,iface,kflag,nk,ipoface(*),
     &  nodface(5,*),index

      integer nboun,ndirboun(*),nodeboun(*),
     &  nmpc,ipompc(*),nodempc(3,*),dof,
     &  ikboun(*),ilboun(*),ikmpc(*),ilmpc(*),
     &  nslavspc(2,*),islavspc(2,*),nsspc,nslavmpc(2,*),islavmpc(2,*),
     &  nsmpc,nmastspc(2,*),imastspc(2,*),nmspc,nmastmpc(2,*),
     &  imastmpc(2,*),nmmpc,islavborder(*),isspc,imspc,ismpc,immpc,
     &  nodem,nodes,ist
      
      real*8  xboun(*),coefmpc(*)
         
       debug=.false.

       isspc=0
       ismpc=0
       i=0
       do i=1,ntie
          do l=nslavnode(i)+1,nslavnode(i+1)
                node=islavnode(l)
!               check for SPCs
c                write(*,*)'SPC'
                nslavspc(1,l)=isspc
                do k=1,3
                 dof=8*(node-1)+k
                 call nident(ikboun, dof, 
     &                        nboun, id)
c                 write(*,*)'nodes',k,node,dof,ikboun(id)
                 if(id>0)then
                 if(dof.eq.ikboun(id))then
                  isspc=isspc+1
                  islavspc(1,isspc)=ilboun(id)

                 endif
                 endif
                enddo
                nslavspc(2,l)=isspc
!               check for MPCs
c                write(*,*)'MPC'
                nslavmpc(1,l)=ismpc
                do k=1,3
                 dof=8*(node-1)+k
                 call nident(ikmpc, dof, 
     &                        nmpc, id)
c                 write(*,*)'nodes',k,node,dof,ikmpc(id)
                 if(id>0)then
                 if(dof.eq.ikmpc(id))then
                  ismpc=ismpc+1
                  islavmpc(1,ismpc)=ipompc(ilmpc(id))
                 endif
                 endif
                enddo
                nslavmpc(2,l)=ismpc
          enddo
       enddo
       nsspc=isspc
       nsmpc=ismpc

       imspc=0
       immpc=0
       do i=1,ntie
          do l=nmastnode(i)+1,nmastnode(i+1)
                node=imastnode(l)
!               check for SPCs
c                write(*,*)'SPC'
                nmastspc(1,l)=imspc
                do k=1,3
                 dof=8*(node-1)+k
                 call nident(ikboun, dof, 
     &                        nboun, id)
c                 write(*,*)'nodem',k,node,dof,ikboun(id)
                 if(id>0)then
                 if(dof.eq.ikboun(id))then
                  imspc=imspc+1
                  imastspc(1,imspc)=ilboun(id)

                 endif
                 endif
                enddo
                nmastspc(2,l)=imspc
!               check for MPCs
c                write(*,*)'MPC'
                nmastmpc(1,l)=immpc
                do k=1,3
                 dof=8*(node-1)+k
                 call nident(ikmpc, dof, 
     &                        nmpc, id)
c                 write(*,*)'nodem',k,node,dof,ikmpc(id)
                 if(id>0)then
                 if(dof.eq.ikmpc(id))then
                  immpc=immpc+1
                  imastmpc(1,immpc)=ipompc(ilmpc(id))
                 endif
                 endif
                enddo
                nmastmpc(2,l)=immpc
          enddo
       enddo
       nmspc=imspc
       nmmpc=immpc

c       isspc=0
c       i=0
c       do 
c          i=i+1
c          if(i.gt.nboun) exit
c          node=nodeboun(ilboun(i))
c          call nident(islavnode, node, 
c     &              nslavnode(ntie+1), id)c
c
c          if(id.gt.0 .and. islavnode(id).eq.node) then
c           nodes=node
c           nslavspc(1,id)=isspc
c           isspc=isspc+1
c           islavspc(1,isspc)=ilboun(i)
c           do
c            node=nodeboun(ilboun(i+1))
c            if(node.ne.nodes)then
c              nslavspc(2,id)=isspc
c              exit
c            endif
c            i=i+1
c            isspc=isspc+1
c            islavspc(1,isspc)=ilboun(i)
c           enddo 
c          endif        
c       enddo 
c       nsspc=isspc
c
c       imspc=0
c       i=0
c       do 
c          i=i+1
c          if(i.gt.nboun) exit
c          node=nodeboun(ilboun(i))
c          call nident(imastnode, node, 
c     &              nmastnode(ntie+1), id)
c          if(id.gt.0 .and. imastnode(id).eq.node) then
c          nodem=node
c          nmastspc(1,id)=imspc
c          imspc=imspc+1
c          imastspc(imspc)=ilboun(i)
c          do
c           node=nodeboun(ilboun(i+1))
c           if(node.ne.nodem)then
c              nmastspc(2,id)=imspc
c              exit
c           endif
c           i=i+1
c           imspc=imspc+1
c           imastspc(imspc)=ilboun(i)
c          enddo
c          endif         
c       enddo 
c       nmspc=imspc
c
c       ismpc=0
c       i=0
c       do 
c          i=i+1
c          if(i.gt.nmpc) exit
c          ist=ipompc(ilmpc(i))
c          node=nodempc(1,ist)
c          call nident(islavnode, node, 
c     &              nslavnode(ntie+1), id)
c
c          write(*,*)'i',i,id,'node',node,islavnode(id)
c          if(id>0)then
c          if(islavnode(id).eq.node) then
c           nodes=node
c           nslavmpc(1,id)=ismpc
c           ismpc=ismpc+1
c           islavmpc(1,ismpc)=ipompc(ilmpc(i))
c           do
c            if(i+1.le.nmpc)then
c             ist=ipompc(ilmpc(i+1))
c             node=nodempc(1,ist)
c             if(node.ne.nodes)then
c              nslavmpc(2,id)=ismpc
c              exit
c             endif
c             i=i+1
c             ismpc=ismpc+1
c             islavmpc(1,ismpc)=ipompc(ilmpc(i))
c            else
c             exit
c            endif
c           enddo
c          endif       
c          endif  
c       enddo 
c       nsmpc=ismpc
c
c       immpc=0
c       i=0
c       do 
c          i=i+1
c          if(i.gt.nmpc) exit
c          ist=ipompc(ilmpc(i))
c          node=nodempc(1,ist)
c          call nident(imastnode, node, 
c     &              nmastnode(ntie+1), id)
c          if(id>0) then
c          if(imastnode(id).eq.node) then
c          nodem=node
c          nmastmpc(1,id)=immpc
c          immpc=immpc+1
c          imastmpc(immpc)=ipompc(ilmpc(i))
c          do
c          if(i+1.le.nmpc)then
c           ist=ipompc(ilmpc(i+1))
c           node=nodempc(1,ist)
c           if(node.ne.nodem)then
c              nmastmpc(2,id)=immpc
c              exit
c           endif
c           i=i+1
c           immpc=immpc+1
c           imastmpc(immpc)=ipompc(ilmpc(i))
c           else
c            exit
c           endif
c          enddo
c          endif
c          endif         
c       enddo 
c       nmmpc=immpc
       
       if(debug)then
       do i=1,ntie
          do l=nslavnode(i)+1,nslavnode(i+1)
           node=islavnode(l)
           write(*,*)'***node',node,'***'
           write(*,*) 'nodes-spc',node, nslavspc(1,l),nslavspc(2,l)
           do j=nslavspc(1,l)+1,nslavspc(2,l)
            ist=islavspc(1,j)
            write(*,*)' spc', ist, nodeboun(ist),ndirboun(ist)
           enddo
           write(*,*) 'nodes-mpc',node, nslavmpc(1,l),nslavmpc(2,l)
           do j=nslavmpc(1,l)+1,nslavmpc(2,l)
            ist=islavmpc(1,j)
            index=nodempc(3,ist)
            write(*,*)' mpc', ist, nodempc(2,ist), nodempc(1,index)
           enddo
          enddo
       enddo
       do i=1,ntie
          do l=nmastnode(i)+1,nmastnode(i+1)
           node=imastnode(l)
           write(*,*)'***node',node,'***'
           write(*,*) 'nodem-spc',node, nmastspc(1,l),nmastspc(2,l)
           do j=nmastspc(1,l)+1,nmastspc(2,l)
            ist=imastspc(1,j)
            write(*,*)'spc', ist, nodeboun(ist),ndirboun(ist)
           enddo
           write(*,*) 'nodem-mpc',node, nmastmpc(1,l),nmastmpc(2,l)
           do j=nmastmpc(1,l)+1,nmastmpc(2,l)
            ist=imastmpc(1,j)
            index=nodempc(3,ist)
            write(*,*)'mpc', ist, nodempc(2,ist), nodempc(1,index)
           enddo
          enddo
       enddo
       endif

c       stop
       return
       end