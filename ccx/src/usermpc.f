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
      subroutine usermpc(ipompc,nodempc,coefmpc,
     &  labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,nk,nk_,nodeboun,ndirboun,
     &  ikboun,ilboun,nboun,nboun_,xboun,inode,node,co,label,typeboun)
!
!     initializes mpc fields for a user MPC
!
      implicit none
!
      character*1 typeboun(*)
      character*20 labmpc(*),label
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,nk,nk_,ikmpc(*),
     &  ilmpc(*),node,id,mpcfreeold,idof,l,nodeboun(*),
     &  ndirboun(*),ikboun(*),ilboun(*),nboun,nboun_,inode
!
      real*8 coefmpc(3,*),co(3,*),xboun(*)
!
      if(node.ne.0) then
         if(inode.eq.1) then
!
!           define a new MPC
!           default for the dependent DOF direction is 1
!
            idof=7*(node-1)+1
!
            call nident(ikmpc,idof,nmpc,id)
            if(id.gt.0) then
               if(ikmpc(id).eq.idof) then
                  write(*,*) '*WARNING in usermpc: DOF for node ',node
                  write(*,*) '         in direction 1 has been used'
                  write(*,*) '         on the dependent side of another'
                  write(*,*) '         MPC. ',label
                  write(*,*) '         constraint cannot be applied'
                  return
               endif
            endif
            nmpc=nmpc+1
            if(nmpc.gt.nmpc_) then
               write(*,*) '*ERROR in usermpc: increase nmpc_'
               stop
            endif
!
            ipompc(nmpc)=mpcfree
            labmpc(nmpc)=label
!
            do l=nmpc,id+2,-1
               ikmpc(l)=ikmpc(l-1)
               ilmpc(l)=ilmpc(l-1)
            enddo
            ikmpc(id+1)=idof
            ilmpc(id+1)=nmpc
         endif
!
!        general case: add a term to the MPC
!
         nodempc(1,mpcfree)=node
c
         if(inode.eq.1) then
            nodempc(2,mpcfree)=1
         else
            nodempc(2,mpcfree)=0
         endif
c         if(label(1:7).eq.'MEANROT') 
c     &         nodempc(2,mpcfree)=inode-3*int((inode-0.5)/3.)
         mpcfree=nodempc(3,mpcfree)
      else
!
!        MPC definition finished: add a nonhomogeneous term
!
         nk=nk+1
         if(nk.gt.nk_) then
            write(*,*) '*ERROR in usermpc: increase nk_'
            stop
         endif
!
         nodempc(1,mpcfree)=nk
         nodempc(2,mpcfree)=1
         mpcfreeold=mpcfree
         mpcfree=nodempc(3,mpcfree)
         nodempc(3,mpcfreeold)=0
         idof=7*(nk-1)+1
         call nident(ikboun,idof,nboun,id)
         nboun=nboun+1
         if(nboun.gt.nboun_) then
            write(*,*) '*ERROR in usermpc: increase nboun_'
            stop
         endif
         nodeboun(nboun)=nk
         ndirboun(nboun)=1
         typeboun(nboun)='U'
         do l=nboun,id+2,-1
            ikboun(l)=ikboun(l-1)
            ilboun(l)=ilboun(l-1)
         enddo
         ikboun(id+1)=idof
         ilboun(id+1)=nboun
      endif
!
      return
      end


