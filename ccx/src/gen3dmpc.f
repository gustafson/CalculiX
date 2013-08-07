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
      subroutine gen3dmpc(ipompc,nodempc,coefmpc,nmpc,nmpc_,mpcfree,
     &  ikmpc,ilmpc,labmpc,iponoel,inoel,iponoelmax,kon,ipkon,lakon,
     &  ne,iponor,xnor,knor,rig)
!
!     connects nodes of 1-D and 2-D elements, for which MPC's were
!     defined, to the nodes of their expanded counterparts
!
      implicit none
!
      character*8 lakon(*)
      character*20 labmpc(*)
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,ikmpc(*),
     &  ilmpc(*),iponoel(*),inoel(3,*),iponoelmax,kon(*),ipkon(*),
     &  ne,iponor(2,*),knor(*),rig(*),i,index1,node,index2,ielem,
     &  indexe,j,indexk,newnode,idir,idof,id,mpcfreenew,k
!
      real*8 coefmpc(*),xnor(*)
!
      do i=1,nmpc
         index1=ipompc(i)
         do
            node=nodempc(1,index1)
            if(node.le.iponoelmax) then
               if(rig(node).ne.0) then
c                  if(nodempc(2,index1).gt.3) then
                  if(nodempc(2,index1).gt.4) then
                     if(rig(node).lt.0) then
                        write(*,*) '*ERROR in gen3dmpc: in node ',node
                        write(*,*) '  a rotational DOF is constrained'
                        write(*,*) '  by a SPC; however, the elements'
                        write(*,*) '  to which this node belongs do not'
                        write(*,*) '  have rotational DOFs'
                        stop
                     endif
                     nodempc(1,index1)=rig(node)
c                     nodempc(2,index1)=nodempc(2,index1)-3
                     nodempc(2,index1)=nodempc(2,index1)-4
                  endif
               else
                  index2=iponoel(node)
c
c                 check for nodes not belonging to 1d or 2d elements
c
                  if(index2.eq.0) then
                     index1=nodempc(3,index1)
                     if(index1.eq.0) exit
                     cycle
                  endif
c
                  ielem=inoel(1,index2)
                  indexe=ipkon(ielem)
                  j=inoel(2,index2)
                  indexk=iponor(2,indexe+j)
!
!                    2d element shell element
!
                  if(lakon(ielem)(7:7).eq.'L') then
                     newnode=knor(indexk+1)
                     idir=nodempc(2,index1)
                     idof=8*(newnode-1)+idir
                     call nident(ikmpc,idof,nmpc,id)
                     if((id.le.0).or.(ikmpc(id).ne.idof)) then
                        nmpc=nmpc+1
                        if(nmpc.gt.nmpc_) then
                           write(*,*) 
     &                          '*ERROR in gen3dmpc: increase nmpc_'
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
     &                          '*ERROR in gen3dmpc: increase nmpc_'
                           stop
                        endif
                        nodempc(1,mpcfree)=knor(indexk+3)
                        nodempc(2,mpcfree)=idir
                        coefmpc(mpcfree)=1.d0
                        mpcfree=nodempc(3,mpcfree)
                        if(mpcfree.eq.0) then
                           write(*,*) 
     &                          '*ERROR in gen3dmpc: increase nmpc_'
                           stop
                        endif
                        nodempc(1,mpcfree)=node
                        nodempc(2,mpcfree)=idir
                        coefmpc(mpcfree)=-2.d0
                        mpcfreenew=nodempc(3,mpcfree)
                        if(mpcfreenew.eq.0) then
                           write(*,*) 
     &                          '*ERROR in gen3dmpc: increase nmpc_'
                           stop
                        endif
                        nodempc(3,mpcfree)=0
                        mpcfree=mpcfreenew
                     endif
                  elseif(lakon(ielem)(7:7).eq.'B') then
!
!                       1d beam element
!
                     newnode=knor(indexk+1)
                     idir=nodempc(2,index1)
                     idof=8*(newnode-1)+idir
                     call nident(ikmpc,idof,nmpc,id)
                     if((id.le.0).or.(ikmpc(id).ne.idof)) then
                        nmpc=nmpc+1
                        if(nmpc.gt.nmpc_) then
                           write(*,*) 
     &                          '*ERROR in gen3dmpc: increase nmpc_'
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
     &                          '*ERROR in gen3dmpc: increase nmpc_'
                           stop
                        endif
                        do k=2,4
                           nodempc(1,mpcfree)=knor(indexk+k)
                           nodempc(2,mpcfree)=idir
                           coefmpc(mpcfree)=1.d0
                           mpcfree=nodempc(3,mpcfree)
                           if(mpcfree.eq.0) then
                              write(*,*) 
     &                             '*ERROR in gen3dmpc: increase nmpc_'
                              stop
                           endif
                        enddo
                        nodempc(1,mpcfree)=node
                        nodempc(2,mpcfree)=idir
                        coefmpc(mpcfree)=-4.d0
                        mpcfreenew=nodempc(3,mpcfree)
                        if(mpcfreenew.eq.0) then
                           write(*,*) 
     &                          '*ERROR in gen3dmpc: increase nmpc_'
                           stop
                        endif
                        nodempc(3,mpcfree)=0
                        mpcfree=mpcfreenew
                     endif
                  else
!
!                       2d plane stress, plane strain or axisymmetric
!                       element
!
                     newnode=knor(indexk+2)
                     idir=nodempc(2,index1)
                     idof=8*(newnode-1)+idir
                     call nident(ikmpc,idof,nmpc,id)
                     if(((id.le.0).or.(ikmpc(id).ne.idof)).and.
     &                    (idir.ne.3)) then
                        nmpc=nmpc+1
                        if(nmpc.gt.nmpc_) then
                           write(*,*) 
     &                          '*ERROR in gen3dmpc: increase nmpc_'
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
     &                          '*ERROR in gen3dmpc: increase nmpc_'
                           stop
                        endif
                        nodempc(1,mpcfree)=node
                        nodempc(2,mpcfree)=idir
                        coefmpc(mpcfree)=-1.d0
                        mpcfreenew=nodempc(3,mpcfree)
                        if(mpcfreenew.eq.0) then
                           write(*,*) 
     &                          '*ERROR in gen3dmpc: increase nmpc_'
                           stop
                        endif
                        nodempc(3,mpcfree)=0
                        mpcfree=mpcfreenew
                     endif
                  endif
               endif
            endif
            index1=nodempc(3,index1)
            if(index1.eq.0) exit
         enddo
      enddo
!
      return
      end


