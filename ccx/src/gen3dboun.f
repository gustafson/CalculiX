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
      subroutine gen3dboun(ikboun,ilboun,nboun,nboun_,nodeboun,ndirboun,
     &  xboun,iamboun,typeboun,iponoel,inoel,iponoelmax,kon,ipkon,
     &  lakon,ne,iponor,xnor,knor,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &  mpcfree,ikmpc,ilmpc,labmpc,rig,ntrans,inotr,trab,nam,nk,nk_,co,
     &  nmethod,iperturb)
!
!     connects nodes of 1-D and 2-D elements, for which SPC's were
!     defined, to the nodes of their expanded counterparts
!
      implicit none
!
      character*1 type,typeboun(*)
      character*8 lakon(*)
      character*20 labmpc(*)
!
      integer ikboun(*),ilboun(*),nboun,nboun_,nodeboun(*),ndirboun(*),
     &  iamboun(*),iponoel(*),inoel(3,*),iponoelmax,kon(*),ipkon(*),ne,
     &  iponor(2,*),knor(*),ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,
     &  ikmpc(*),ilmpc(*),rig(*),ntrans,inotr(2,*),nbounold,i,node,
     &  index,ielem,j,indexe,indexk,idir,iamplitude,irotnode,nk,nk_,
     &  newnode,idof,id,mpcfreenew,k,nam,nmethod,iperturb
!
      real*8 xboun(*),xnor(*),coefmpc(*),trab(7,*),val,co(3,*)
!
      nbounold=nboun
      do i=1,nbounold
         node=nodeboun(i)
         if(node.gt.iponoelmax) then
            if(ndirboun(i).gt.3) then
               write(*,*) '*WARNING: in gen3dboun: node ',i,
     &              ' does not'
               write(*,*) '       belong to a beam nor shell'
               write(*,*) '       element and consequently has no'
               write(*,*) '       rotational degrees of freedom'
            endif
            cycle
         endif
         index=iponoel(node)
         if(index.eq.0) then
            if(ndirboun(i).gt.3) then
               write(*,*) '*WARNING: in gen3dboun: node ',i,
     &              ' does not'
               write(*,*) '       belong to a beam nor shell'
               write(*,*) '       element and consequently has no'
               write(*,*) '       rotational degrees of freedom'
            endif
            cycle
         endif
         ielem=inoel(1,index)
         j=inoel(2,index)
         indexe=ipkon(ielem)
         indexk=iponor(2,indexe+j)
         idir=ndirboun(i)
         val=xboun(i)
         if(nam.gt.0) iamplitude=iamboun(i)
!
         if(rig(node).ne.0) then
            if(idir.gt.3) then
               if(rig(node).lt.0) then
                  write(*,*) '*ERROR in gen3dboun: in node ',node
                  write(*,*) '       a rotational DOF is constrained'
                  write(*,*) '       by a SPC; however, the elements'
                  write(*,*) '       to which this node belongs do not'
                  write(*,*) '       have rotational DOFs'
                  stop
               endif
               j=idir-3
               irotnode=rig(node)
               type='B'
               call bounadd(irotnode,j,j,val,nodeboun,
     &              ndirboun,xboun,nboun,nboun_,iamboun,
     &              iamplitude,nam,ipompc,nodempc,coefmpc,
     &              nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &              ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &              type,typeboun,nmethod,iperturb)
            endif
         else
!
!                    2d element shell element: generate MPC's
!
            if(lakon(ielem)(7:7).eq.'L') then
               newnode=knor(indexk+1)
               idof=7*(newnode-1)+idir
               call nident(ikmpc,idof,nmpc,id)
               if((id.le.0).or.(ikmpc(id).ne.idof)) then
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                     write(*,*) 
     &                    '*ERROR in gen3dboun: increase nmpc_'
                     stop
                  endif
                  labmpc(nmpc)='                    '
                  ipompc(nmpc)=mpcfree
                  do j=nmpc,id+2,-1
                     ikmpc(j)=ikmpc(j-1)
                     ilmpc(j)=ilmpc(j-1)
                  enddo
                  ikmpc(id+1)=idof
                  ilmpc(id+1)=nmpc
                  nodempc(1,mpcfree)=newnode
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dboun: increase nmpc_'
                     stop
                  endif
                  nodempc(1,mpcfree)=knor(indexk+3)
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dboun: increase nmpc_'
                     stop
                  endif
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-2.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dboun: increase nmpc_'
                     stop
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
               endif
!
!              fixing the temperature degrees of freedom
!
               if(idir.eq.0) then
                  type='B'
                  call bounadd(knor(indexk+3),idir,idir,val,nodeboun,
     &                 ndirboun,xboun,nboun,nboun_,iamboun,
     &                 iamplitude,nam,ipompc,nodempc,coefmpc,
     &                 nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &                 ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &                 type,typeboun,nmethod,iperturb)
               endif
            elseif(lakon(ielem)(7:7).eq.'B') then
!
!                       1d beam element: generate MPC's
!
               newnode=knor(indexk+1)
               idof=7*(newnode-1)+idir
               call nident(ikmpc,idof,nmpc,id)
               if((id.le.0).or.(ikmpc(id).ne.idof)) then
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                     write(*,*) 
     &                    '*ERROR in gen3dboun: increase nmpc_'
                     stop
                  endif
                  labmpc(nmpc)='                    '
                  ipompc(nmpc)=mpcfree
                  do j=nmpc,id+2,-1
                     ikmpc(j)=ikmpc(j-1)
                     ilmpc(j)=ilmpc(j-1)
                  enddo
                  ikmpc(id+1)=idof
                  ilmpc(id+1)=nmpc
                  nodempc(1,mpcfree)=newnode
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dboun: increase nmpc_'
                     stop
                  endif
                  do k=2,4
                     nodempc(1,mpcfree)=knor(indexk+k)
                     nodempc(2,mpcfree)=idir
                     coefmpc(mpcfree)=1.d0
                     mpcfree=nodempc(3,mpcfree)
                     if(mpcfree.eq.0) then
                        write(*,*) 
     &                       '*ERROR in gen3dboun: increase nmpc_'
                        stop
                     endif
                  enddo
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-4.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dboun: increase nmpc_'
                     stop
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
               endif
!
!              fixing the temperature degrees of freedom
!
               if(idir.eq.0) then
                  type='B'
                  do k=2,4
                     call bounadd(knor(indexk+k),idir,idir,val,nodeboun,
     &                    ndirboun,xboun,nboun,nboun_,iamboun,
     &                    iamplitude,nam,ipompc,nodempc,coefmpc,
     &                    nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &                    ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &                    type,typeboun,nmethod,iperturb)
                  enddo
               endif
            else
!
!                       2d plane stress, plane strain or axisymmetric
!                       element: SPC
!
               newnode=knor(indexk+2)
               idof=7*(newnode-1)+idir
               call nident(ikmpc,idof,nmpc,id)
               if(((id.le.0).or.(ikmpc(id).ne.idof)).and.
     &              (idir.ne.3)) then
                  nmpc=nmpc+1
                  if(nmpc.gt.nmpc_) then
                     write(*,*) 
     &                    '*ERROR in gen3dmpc: increase nmpc_'
                     stop
                  endif
                  labmpc(nmpc)='                    '
                  ipompc(nmpc)=mpcfree
                  do j=nmpc,id+2,-1
                     ikmpc(j)=ikmpc(j-1)
                     ilmpc(j)=ilmpc(j-1)
                  enddo
                  ikmpc(id+1)=idof
                  ilmpc(id+1)=nmpc
                  nodempc(1,mpcfree)=newnode
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dmpc: increase nmpc_'
                     stop
                  endif
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-1.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dmpc: increase nmpc_'
                     stop
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
               endif
c     node=knor(indexk+2)
c     type='B'
c     call bounadd(node,idir,idir,val,nodeboun,
c     &              ndirboun,xboun,nboun,nboun_,iamboun,
c     &              iamplitude,nam,ipompc,nodempc,coefmpc,
c     &              nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
c     &              ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
c     &              type,typeboun,nmethod,iperturb)
            endif
         endif
      enddo
!
      return
      end


