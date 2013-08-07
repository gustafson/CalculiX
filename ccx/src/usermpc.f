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
      subroutine usermpc(ipompc,nodempc,coefmpc,
     &  labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,nk,nk_,nodeboun,ndirboun,
     &  ikboun,ilboun,nboun,nboun_,inode,node,co,label,typeboun,
     &  iperturb)
!
!     initializes mpc fields for a user MPC
!
      implicit none
!
      character*1 typeboun(*)
      character*20 labmpc(*),label
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,nk,nk_,
     &  ikmpc(*),
     &  ilmpc(*),node,id,mpcfreeold,idof,l,nodeboun(*),iperturb(2),
     &  ndirboun(*),ikboun(*),ilboun(*),nboun,nboun_,inode,nodevector,
     &  index,index1,node1,i,j,imax,nkn,idofrem,idofins
!
      real*8 coefmpc(*),co(3,*),aa(3),dd,cgx(3),pi(3),c1,c4,c9,
     &  c10,a(3),a1,amax
!
      save nodevector
!
      if(node.ne.0) then
         if(inode.eq.1) then
!
!           define a new MPC
!           default for the dependent DOF direction is 1
!
            idof=8*(node-1)+1
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
!
!        nodevector: additional node such that:
!         - the coordinates of this node are the axis direction
!         - the 1st DOF is reserved for the mean rotation value
!
         if((labmpc(nmpc)(1:7).eq.'MEANROT').or.
     &        (labmpc(nmpc)(1:1).eq.'1')) then
            nodevector=node
            labmpc(nmpc)(1:7)='MEANROT'
         endif
!
         if(inode.eq.1) then
            nodempc(2,mpcfree)=1
         else
            nodempc(2,mpcfree)=0
         endif
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
c
         coefmpc(mpcfree)=1.d0
c
         mpcfreeold=mpcfree
         mpcfree=nodempc(3,mpcfree)
         nodempc(3,mpcfreeold)=0
         idof=8*(nk-1)+1
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
!
!        calculating the MPC coefficients for linear applications
!
         if((labmpc(nmpc)(1:7).eq.'MEANROT').or.
     &        (labmpc(nmpc)(1:1).eq.'1')) then
            nkn=(inode-1)/3
            if(3*nkn.ne.inode-1) then
               write(*,*)
     &              '*ERROR in usermpc: MPC has wrong number of terms'
               stop
            endif
!
!           normal along the rotation axis
!
            dd=0.d0
            do i=1,3
               aa(i)=co(i,nodevector)
               dd=dd+aa(i)**2
            enddo
            dd=dsqrt(dd)
            if(dd.lt.1.d-10) then
               write(*,*) 
     &           '*ERROR in usermpc: rotation vector has zero length'
               stop
            endif
            do i=1,3
               aa(i)=aa(i)/dd
            enddo
!     
!     finding the center of gravity of the position and the
!     displacements of the nodes involved in the MPC
!
            do i=1,3
               cgx(i)=0.d0
            enddo
!
            index=ipompc(nmpc)
            do
               node=nodempc(1,index)
               if(node.eq.nodevector) exit
               do j=1,3
                  cgx(j)=cgx(j)+co(j,node)
               enddo
               index=nodempc(3,nodempc(3,nodempc(3,index)))
            enddo
!
            do i=1,3
               cgx(i)=cgx(i)/nkn
            enddo
!
!           calculating the derivatives
!
            index=ipompc(nmpc)
            do
               node=nodempc(1,index)
               if(node.eq.nodevector) exit
!
!              relative positions
!               
               do j=1,3
                  pi(j)=co(j,node)-cgx(j)
               enddo
               c1=pi(1)*pi(1)+pi(2)*pi(2)+pi(3)*pi(3)
               if(c1.lt.1.d-20) then
                  write(*,*)'*WARNING in usermpc: node on rotation axis'
                  index=nodempc(3,nodempc(3,nodempc(3,index)))
                  cycle
               endif
!
               do j=1,3
                  if(j.eq.1) then
                     c4=aa(2)*pi(3)-aa(3)*pi(2)
                  elseif(j.eq.2) then
                     c4=aa(3)*pi(1)-aa(1)*pi(3)
                  else
                     c4=aa(1)*pi(2)-aa(2)*pi(1)
                  endif
                  c9=c4/c1
!
                  index1=ipompc(nmpc)
                  do
                     node1=nodempc(1,index1)
                     if(node1.eq.nodevector) exit
                     if(node1.eq.node) then
                        c10=c9*(1.d0-1.d0/real(nkn))
                     else
                        c10=-c9/real(nkn)
                     endif
                     if(j.eq.1) then
                        coefmpc(index1)=coefmpc(index1)+c10
                     elseif(j.eq.2) then
                        coefmpc(nodempc(3,index1))=
     &                  coefmpc(nodempc(3,index1))+c10
                     else
                        coefmpc(nodempc(3,nodempc(3,index1)))=
     &                  coefmpc(nodempc(3,nodempc(3,index1)))+c10
                     endif
                     index1=nodempc(3,nodempc(3,nodempc(3,index1)))
                  enddo
               enddo
               index=nodempc(3,nodempc(3,nodempc(3,index)))
            enddo
            coefmpc(index)=-nkn
!
!     assigning the degrees of freedom
!
            j=0
            index=ipompc(nmpc)
            do
               j=j+1
               if(j.gt.3) j=1
               nodempc(2,index)=j
               index=nodempc(3,index)
               if(nodempc(1,index).eq.nk) exit
            enddo
!
!     looking for the maximum tangent to decide which DOF should be
!     taken to be the dependent one
!
            index=ipompc(nmpc)
            if(dabs(coefmpc(index)).lt.1.d-5) then
!
!              changing the DOF of the dependent degree of freedom
!
               amax=dabs(coefmpc(index))
               imax=1
               a(1)=coefmpc(index)
               do i=2,3
                  index=nodempc(3,index)
                  a(i)=coefmpc(index)
                  if(dabs(a(i)).gt.amax) then
                     amax=dabs(a(i))
                     imax=i
                  endif
               enddo
!
               index=ipompc(nmpc)
               nodempc(2,index)=imax
               a1=a(1)
               coefmpc(index)=a(imax)
               do i=2,3
                  index=nodempc(3,index)
                  if(i.eq.imax) then
                     nodempc(2,index)=1
                     coefmpc(index)=a1
                  else
                     nodempc(2,index)=i
                  endif
               enddo
!
!              updating ikmpc and ilmpc
!               
               index=ipompc(nmpc)
               idofrem=8*(nodempc(1,index)-1)+1
               idofins=8*(nodempc(1,index)-1)+imax
               call changedepterm(ikmpc,ilmpc,nmpc,nmpc,idofrem,idofins)
            endif
         elseif(labmpc(nmpc)(1:4).eq.'DIST') then
            iperturb(2)=1
            write(*,*) '*INFO in usermpc: nonlinear geometric'
            write(*,*) '      effects are turned on'
            write(*,*)
            if(iperturb(1).eq.0) iperturb(1)=2
         elseif(labmpc(nmpc)(1:3).eq.'GAP') then
            iperturb(2)=1
            write(*,*) '*INFO in usermpc: nonlinear geometric'
            write(*,*) '      effects are turned on'
            write(*,*)
            if(iperturb(1).eq.0) iperturb(1)=2
         elseif(labmpc(nmpc)(1:4).eq.'USER') then
            iperturb(2)=1
            write(*,*) '*INFO in usermpc: nonlinear geometric'
            write(*,*) '      effects are turned on'
            write(*,*)
            if(iperturb(1).eq.0) iperturb(1)=2
         else
            write(*,*) '*ERROR in usermpc: mpc of type',labmpc(nmpc)
            write(*,*) '       is unknown'
            stop
         endif
      endif
!
      return
      end


