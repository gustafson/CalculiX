!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine gencycsymelemcfd(ics,cs,icscp,xcs,ycs,zcs,islav,
     &  nslav,islavcp,xslav,yslav,zslav,nface,nelemface,sideface,
     &  nk,co,ne,ipkon,lakon,kon,nkon,ielmat,mi)
!
!     creates an additional layer of elements on the slave and the
!     master surface of a cyclic symmetric structure
!
      implicit none
!
      logical master,slave
!
      character*1 sideface(*)
      character*8 lakon(*)
!
      integer ics(*),i,j,ncs,node,icscp(*),islav(*),nslav,islavcp(*),
     &  nface,nelemface(*),nk,ne,ipkon(*),kon(*),nkon,nelem,iface,id,
     &  ifaceq(8,6),ifacet(6,4),ifacew(8,5),itopology(8),jmax,mi(*),
     &  ielmat(mi(3),*),indexe
!
      real*8 cs(17,*),phibasis,phi,xa(3),xn(3),dd,cphi,sphi,dphi,
     &  xcs(*),ycs(*),zcs(*),c(3,3),d(3,3),xp(3),xq(3),xslav(*),
     &  yslav(*),zslav(*),al,co(3,*)
!
!     nodes belonging to the element faces
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
!
!     mcs=1 is assumed
!     phi is the angle of the cyclic basis sector
!
      phibasis=360.d0/cs(1,1)
!
!     number of nodes on the independent side
!
      ncs=int(cs(4,1))
!
!     xn is an normed vector along the axis
!     xa is a point on the axis
!
      do i=1,3
         xa(i)=cs(5+i,1)
         xn(i)=cs(8+i,1)-xa(i)
      enddo
      dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
      do i=1,3
         xn(i)=xn(i)/dd
      enddo
!
!     new nodes on the master side
!
!     phi is the thickness of the layer in radians
!
c      phi=phibasis/10.d0
      phi=phibasis
      cphi=dcos(phi)
      sphi=dsin(phi)
      dphi=1.d0-cphi
!
!     rotation matrix
!
      c(1,1)=cphi+dphi*xn(1)*xn(1)
      c(1,2)=     dphi*xn(1)*xn(2)-sphi*xn(3)
      c(1,3)=     dphi*xn(1)*xn(3)+sphi*xn(2)
      c(2,1)=     dphi*xn(2)*xn(1)+sphi*xn(3)
      c(2,2)=cphi+dphi*xn(2)*xn(2)
      c(2,3)=     dphi*xn(2)*xn(3)-sphi*xn(1)
      c(3,1)=     dphi*xn(3)*xn(1)-sphi*xn(2)
      c(3,2)=     dphi*xn(3)*xn(2)+sphi*xn(1)
      c(3,3)=cphi+dphi*xn(3)*xn(3)
!
!     cyclic position of these nodes
!
c      phi=-9.d0*phibasis/10.d0
      phi=0.d0
      cphi=dcos(phi)
      sphi=dsin(phi)
      dphi=1.d0-cphi
!
!     rotation matrix
!
      d(1,1)=cphi+dphi*xn(1)*xn(1)
      d(1,2)=     dphi*xn(1)*xn(2)-sphi*xn(3)
      d(1,3)=     dphi*xn(1)*xn(3)+sphi*xn(2)
      d(2,1)=     dphi*xn(2)*xn(1)+sphi*xn(3)
      d(2,2)=cphi+dphi*xn(2)*xn(2)
      d(2,3)=     dphi*xn(2)*xn(3)-sphi*xn(1)
      d(3,1)=     dphi*xn(3)*xn(1)-sphi*xn(2)
      d(3,2)=     dphi*xn(3)*xn(2)+sphi*xn(1)
      d(3,3)=cphi+dphi*xn(3)*xn(3)
!
      do i=1,ncs
         node=ics(i)
         do j=1,3
            xp(j)=co(j,node)
         enddo
         al=(xp(1)-xa(1))*xn(1)+
     &      (xp(2)-xa(2))*xn(2)+
     &      (xp(3)-xa(3))*xn(3)
!
!        xq is the orthogonal projection of xp on the axis
!
         do j=1,3
            xq(j)=xa(j)+al*xn(j)
         enddo
!
         nk=nk+1
         icscp(i)=nk
         do j=1,3
            co(j,nk)=xq(j)+c(j,1)*(xp(1)-xq(1))+
     &                     c(j,2)*(xp(2)-xq(2))+
     &                     c(j,3)*(xp(3)-xq(3))
         enddo
         xcs(i)=xq(1)+d(1,1)*(xp(1)-xq(1))+
     &        d(1,2)*(xp(2)-xq(2))+
     &        d(1,3)*(xp(3)-xq(3))
         ycs(i)=xq(2)+d(2,1)*(xp(1)-xq(1))+
     &        d(2,2)*(xp(2)-xq(2))+
     &        d(2,3)*(xp(3)-xq(3))
         zcs(i)=xq(3)+d(3,1)*(xp(1)-xq(1))+
     &        d(3,2)*(xp(2)-xq(2))+
     &        d(3,3)*(xp(3)-xq(3))
      enddo
!
!     new nodes on the slave side
!
!     phi is the thickness of the layer in radians
!
c      phi=-phibasis/10.d0
      phi=-phibasis
      cphi=dcos(phi)
      sphi=dsin(phi)
      dphi=1.d0-cphi
!
!     rotation matrix
!
      c(1,1)=cphi+dphi*xn(1)*xn(1)
      c(1,2)=     dphi*xn(1)*xn(2)-sphi*xn(3)
      c(1,3)=     dphi*xn(1)*xn(3)+sphi*xn(2)
      c(2,1)=     dphi*xn(2)*xn(1)+sphi*xn(3)
      c(2,2)=cphi+dphi*xn(2)*xn(2)
      c(2,3)=     dphi*xn(2)*xn(3)-sphi*xn(1)
      c(3,1)=     dphi*xn(3)*xn(1)-sphi*xn(2)
      c(3,2)=     dphi*xn(3)*xn(2)+sphi*xn(1)
      c(3,3)=cphi+dphi*xn(3)*xn(3)
!
!     cyclic position of these nodes
!
c      phi=9.d0*phibasis/10.d0
      phi=0.d0
      cphi=dcos(phi)
      sphi=dsin(phi)
      dphi=1.d0-cphi
!
!     rotation matrix
!
      d(1,1)=cphi+dphi*xn(1)*xn(1)
      d(1,2)=     dphi*xn(1)*xn(2)-sphi*xn(3)
      d(1,3)=     dphi*xn(1)*xn(3)+sphi*xn(2)
      d(2,1)=     dphi*xn(2)*xn(1)+sphi*xn(3)
      d(2,2)=cphi+dphi*xn(2)*xn(2)
      d(2,3)=     dphi*xn(2)*xn(3)-sphi*xn(1)
      d(3,1)=     dphi*xn(3)*xn(1)-sphi*xn(2)
      d(3,2)=     dphi*xn(3)*xn(2)+sphi*xn(1)
      d(3,3)=cphi+dphi*xn(3)*xn(3)
!
      do i=1,nslav
         node=islav(i)
         do j=1,3
            xp(j)=co(j,node)
         enddo
         al=(xp(1)-xa(1))*xn(1)+
     &      (xp(2)-xa(2))*xn(2)+
     &      (xp(3)-xa(3))*xn(3)
!
!        xq is the orthogonal projection of xp on the axis
!
         do j=1,3
            xq(j)=xa(j)+al*xn(j)
         enddo
!
         nk=nk+1
         islavcp(i)=nk
         do j=1,3
            co(j,nk)=xq(j)+c(j,1)*(xp(1)-xq(1))+
     &                     c(j,2)*(xp(2)-xq(2))+
     &                     c(j,3)*(xp(3)-xq(3))
         enddo
         xslav(i)=xq(1)+d(1,1)*(xp(1)-xq(1))+
     &        d(1,2)*(xp(2)-xq(2))+
     &        d(1,3)*(xp(3)-xq(3))
         yslav(i)=xq(2)+d(2,1)*(xp(1)-xq(1))+
     &        d(2,2)*(xp(2)-xq(2))+
     &        d(2,3)*(xp(3)-xq(3))
         zslav(i)=xq(3)+d(3,1)*(xp(1)-xq(1))+
     &        d(3,2)*(xp(2)-xq(2))+
     &        d(3,3)*(xp(3)-xq(3))
      enddo
!
!     determining the slave and master faces
!     generating new elements
!
      do i=1,nface
         nelem=nelemface(i)
         indexe=ipkon(nelem)
         read(sideface(i)(1:1),'(i1)') iface
!
!        tetrahedra
!
         if(lakon(nelem)(4:4).eq.'4') then
            master=.true.
            do j=1,3
               node=kon(indexe+ifacet(j,iface))
               call nident(ics,node,ncs,id)
               if(id.le.0) then
                  master=.false.
                  exit
               elseif(ics(id).ne.node) then
                  master=.false.
                  exit
               endif
               itopology(j)=node
               itopology(j+3)=icscp(id)
            enddo
!
!           nodes on the axis are not considered yet
!
            if(.not.master) then
               slave=.true.
               do j=1,3
                  node=kon(indexe+ifacet(j,iface))
                  call nident(islav,node,nslav,id)
                  if(id.le.0) then
                     slave=.false.
                     exit
                  elseif(islav(id).ne.node) then
                     slave=.false.
                     exit
                  endif
                  itopology(j)=node
                  itopology(j+3)=icscp(id)
               enddo
            endif
!
            if((master).or.(slave)) then
               ne=ne+1
               ipkon(ne)=nkon
               lakon(ne)='F3D6    '
               ielmat(1,ne)=ielmat(1,nelem)
               do j=1,6
                  kon(nkon+j)=itopology(j)
               enddo
               nkon=nkon+6
            endif
!
!           wedges
!
         elseif(lakon(nelem)(4:4).eq.'6') then
            master=.true.
            if(iface.le.2) then
               jmax=3
            else
               jmax=4
            endif
!
            do j=1,jmax
               node=kon(indexe+ifacew(j,iface))
               call nident(ics,node,ncs,id)
               if(id.le.0) then
                  master=.false.
                  exit
               elseif(ics(id).ne.node) then
                  master=.false.
                  exit
               endif
               itopology(j)=node
               itopology(j+jmax)=icscp(id)
            enddo
!
!           nodes on the axis are not considered yet
!
            if(.not.master) then
               slave=.true.
               do j=1,jmax
                  node=kon(indexe+ifacew(j,iface))
                  call nident(islav,node,nslav,id)
                  if(id.le.0) then
                     slave=.false.
                     exit
                  elseif(islav(id).ne.node) then
                     slave=.false.
                     exit
                  endif
                  itopology(j)=node
                  itopology(j+jmax)=icscp(id)
               enddo
            endif
!
            if((master).or.(slave)) then
               ne=ne+1
               ipkon(ne)=nkon
               if(jmax.eq.3) then
                  lakon(ne)='F3D6    '
               else
                  lakon(ne)='F3D8    '
               endif
               ielmat(1,ne)=ielmat(1,nelem)
               do j=1,2*jmax
                  kon(nkon+j)=itopology(j)
               enddo
               nkon=nkon+2*jmax
            endif
!
!           hexahedra
!
         elseif(lakon(nelem)(4:4).eq.'8') then
            master=.true.
            do j=1,4
               node=kon(indexe+ifaceq(j,iface))
               call nident(ics,node,ncs,id)
               if(id.le.0) then
                  master=.false.
                  exit
               elseif(ics(id).ne.node) then
                  master=.false.
                  exit
               endif
               itopology(j)=node
               itopology(j+4)=icscp(id)
            enddo
!
!           nodes on the axis are not considered yet
!
            if(.not.master) then
               slave=.true.
               do j=1,4
                  node=kon(indexe+ifaceq(j,iface))
                  call nident(islav,node,nslav,id)
                  if(id.le.0) then
                     slave=.false.
                     exit
                  elseif(islav(id).ne.node) then
                     slave=.false.
                     exit
                  endif
                  itopology(j)=node
                  itopology(j+4)=icscp(id)
               enddo
            endif
!
            if((master).or.(slave)) then
               ne=ne+1
               ipkon(ne)=nkon
               lakon(ne)='F3D8    '
               ielmat(1,ne)=ielmat(1,nelem)
               do j=1,8
                  kon(nkon+j)=itopology(j)
               enddo
               nkon=nkon+8
            endif
         endif
      enddo
!
      return
      end

