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
      subroutine knotmpc(ipompc,nodempc,coefmpc,irefnode,irotnode,
     &  iexpnode,
     &  labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,nk,nk_,nodeboun,ndirboun,
     &  ikboun,ilboun,nboun,nboun_,node,typeboun,co)
!
!     generates three knot MPC's for node "node" about reference
!     (translational) node irefnode and rotational node irotnode 
!
      implicit none
!
      character*1 typeboun(*)
      character*20 labmpc(*)
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,nk,nk_,ikmpc(*),
     &  ilmpc(*),node,id,mpcfreeold,j,idof,l,nodeboun(*),
     &  ndirboun(*),ikboun(*),ilboun(*),nboun,nboun_,irefnode,
     &  irotnode,iexpnode
!
      real*8 coefmpc(*),co(3,*)
!
c      data e /0.,0.,0.,0.,0.,-1.,0.,1.,0.,
c     &        0.,0.,1.,0.,0.,0.,-1.,0.,0.,
c     &        0.,-1.,0.,1.,0.,0.,0.,0.,0./
!
      nk=nk+1
      if(nk.gt.nk_) then
         write(*,*) '*ERROR in knotmpc: increase nk_'
         stop
      endif
      do j=1,3
         idof=7*(node-1)+j
         call nident(ikmpc,idof,nmpc,id)
         if(id.gt.0) then
            if(ikmpc(id).eq.idof) then
               cycle
            endif
         endif
         nmpc=nmpc+1
         if(nmpc.gt.nmpc_) then
            write(*,*) '*ERROR in knotmpc: increase nmpc_'
            stop
         endif
!
         ipompc(nmpc)=mpcfree
         labmpc(nmpc)='KNOT                '
!
         do l=nmpc,id+2,-1
            ikmpc(l)=ikmpc(l-1)
            ilmpc(l)=ilmpc(l-1)
         enddo
         ikmpc(id+1)=idof
         ilmpc(id+1)=nmpc
!
         nodempc(1,mpcfree)=node
         nodempc(2,mpcfree)=j
         coefmpc(mpcfree)=1.d0
         mpcfree=nodempc(3,mpcfree)
!
!        translation term
!
         nodempc(1,mpcfree)=irefnode
         nodempc(2,mpcfree)=j
         coefmpc(mpcfree)=-1.d0
         mpcfree=nodempc(3,mpcfree)
!
!        expansion term
!
         nodempc(1,mpcfree)=iexpnode
         nodempc(2,mpcfree)=1
         mpcfree=nodempc(3,mpcfree)
!
!        rotation terms
!
         nodempc(1,mpcfree)=irotnode
         nodempc(2,mpcfree)=1
c         coefmpc(mpcfree)=e(j,1,1)*(co(1,irefnode)-co(1,node))+
c     &        e(j,2,1)*(co(2,irefnode)-co(2,node))+
c     &        e(j,3,1)*(co(3,irefnode)-co(3,node))
         mpcfree=nodempc(3,mpcfree)
         nodempc(1,mpcfree)=irotnode
         nodempc(2,mpcfree)=2
c         coefmpc(mpcfree)=e(j,1,2)*(co(1,irefnode)-co(1,node))+
c     &        e(j,2,2)*(co(2,irefnode)-co(2,node))+
c     &        e(j,3,2)*(co(3,irefnode)-co(3,node))
         mpcfree=nodempc(3,mpcfree)
         nodempc(1,mpcfree)=irotnode
         nodempc(2,mpcfree)=3
c         coefmpc(mpcfree)=e(j,1,3)*(co(1,irefnode)-co(1,node))+
c     &        e(j,2,3)*(co(2,irefnode)-co(2,node))+
c     &        e(j,3,3)*(co(3,irefnode)-co(3,node))
         mpcfree=nodempc(3,mpcfree)
         nodempc(1,mpcfree)=nk
         nodempc(2,mpcfree)=j
         coefmpc(mpcfree)=1.d0
         mpcfreeold=mpcfree
         mpcfree=nodempc(3,mpcfree)
         nodempc(3,mpcfreeold)=0
         idof=7*(nk-1)+j
         call nident(ikboun,idof,nboun,id)
         nboun=nboun+1
         if(nboun.gt.nboun_) then
            write(*,*) '*ERROR in knotmpc: increase nboun_'
            stop
         endif
         nodeboun(nboun)=nk
         ndirboun(nboun)=j
         typeboun(nboun)='R'
         do l=nboun,id+2,-1
            ikboun(l)=ikboun(l-1)
            ilboun(l)=ilboun(l-1)
         enddo
         ikboun(id+1)=idof
         ilboun(id+1)=nboun
       enddo
!
       return
       end


