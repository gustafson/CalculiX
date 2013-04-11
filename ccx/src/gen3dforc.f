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
      subroutine gen3dforc(ikforc,ilforc,nforc,nforc_,nodeforc,
     &  ndirforc,xforc,iamforc,ntrans,inotr,trab,rig,ipompc,nodempc,
     &  coefmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,labmpc,iponoel,inoel,
     &  iponoelmax,kon,ipkon,lakon,ne,iponor,xnor,knor,nam,nk,nk_,
     &  co,thicke)
!
!     connects nodes of 1-D and 2-D elements, for which 
!     concentrated forces were
!     defined, to the nodes of their expanded counterparts
!
      implicit none
!
      logical add
!
      character*8 lakon(*)
      character*20 labmpc(*)
!
      integer ikforc(*),ilforc(*),nodeforc(2,*),ndirforc(*),iamforc(*),
     &  nforc,nforc_,ntrans,inotr(2,*),rig(*),ipompc(*),nodempc(3,*),
     &  nmpc,nmpc_,mpcfree,ikmpc(*),ilmpc(*),iponoel(*),inoel(3,*),
     &  iponoelmax,kon(*),ipkon(*),ne,iponor(2,*),knor(*),nforcold,
     &  i,node,index,ielem,j,indexe,indexk,nam,iamplitude,idir,
     &  irotnode,nk,nk_,newnode,idof,id,mpcfreenew,k,isector
!
      real*8 xforc(*),trab(7,*),coefmpc(*),xnor(*),val,co(3,*),
     &  thicke(2,*),pi
!
      add=.false.
      pi=4.d0*datan(1.d0)
      isector=0
!
      nforcold=nforc
      do i=1,nforcold
         node=nodeforc(1,i)
         if(node.gt.iponoelmax) then
            if(ndirforc(i).gt.3) then
               write(*,*) '*WARNING: in gen3dforc: node ',i,
     &              ' does not'
               write(*,*) '       belong to a beam nor shell'
               write(*,*) '       element and consequently has no'
               write(*,*) '       rotational degrees of freedom'
            endif
            cycle
         endif
         index=iponoel(node)
         if(index.eq.0) then
            if(ndirforc(i).gt.3) then
               write(*,*) '*WARNING: in gen3dforc: node ',i,
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
         if(nam.gt.0) iamplitude=iamforc(i)
         idir=ndirforc(i)
!
         if(rig(node).ne.0) then
            if(idir.gt.3) then
               if(rig(node).lt.0) then
                  write(*,*) '*ERROR in gen3dforc: in node ',node
                  write(*,*) '       a rotational DOF is loaded;'
                  write(*,*) '       however, the elements to which'
                  write(*,*) '       this node belongs do not have'
                  write(*,*) '       rotational DOFs'
                  stop
               endif
               val=xforc(i)
               j=idir-3
               irotnode=rig(node)
               call forcadd(irotnode,j,val,nodeforc,
     &              ndirforc,xforc,nforc,nforc_,iamforc,
     &              iamplitude,nam,ntrans,trab,inotr,co,
     &              ikforc,ilforc,isector,add)
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
     &                    '*ERROR in gen3dforc: increase nmpc_'
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
     &                    '*ERROR in gen3dforc: increase nmpc_'
                     stop
                  endif
                  nodempc(1,mpcfree)=knor(indexk+3)
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=1.d0
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dforc: increase nmpc_'
                     stop
                  endif
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-2.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dforc: increase nmpc_'
                     stop
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
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
     &                    '*ERROR in gen3dforc: increase nmpc_'
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
     &                    '*ERROR in gen3dforc: increase nmpc_'
                     stop
                  endif
                  do k=2,4
                     nodempc(1,mpcfree)=knor(indexk+k)
                     nodempc(2,mpcfree)=idir
                     coefmpc(mpcfree)=1.d0
                     mpcfree=nodempc(3,mpcfree)
                     if(mpcfree.eq.0) then
                        write(*,*) 
     &                       '*ERROR in gen3dforc: increase nmpc_'
                        stop
                     endif
                  enddo
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=idir
                  coefmpc(mpcfree)=-4.d0
                  mpcfreenew=nodempc(3,mpcfree)
                  if(mpcfreenew.eq.0) then
                     write(*,*) 
     &                    '*ERROR in gen3dforc: increase nmpc_'
                     stop
                  endif
                  nodempc(3,mpcfree)=0
                  mpcfree=mpcfreenew
               endif
            else
!
!              2d plane strain, plane stress or axisymmetric 
!              element
!
               node=knor(indexk+2)
               val=xforc(i)
c               if(lakon(ielem)(7:7).eq.'A') then
c                  val=val*thicke(1,indexe+j)/(2.d0*pi)
c               endif
               call forcadd(node,idir,val,nodeforc,
     &              ndirforc,xforc,nforc,nforc_,iamforc,
     &              iamplitude,nam,ntrans,trab,inotr,co,
     &              ikforc,ilforc,isector,add)
            endif
         endif
      enddo
!
      return
      end


