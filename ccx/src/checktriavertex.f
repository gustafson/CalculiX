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
      subroutine checktriavertex(inodesin,nnodesin,node,nvertex,pvertex,
     &  lvertex,pnodesin,inodesout,nnodesout,nopes,slavstraight,
     &  xn,co,xl2,vold,mi)
!
!     check whether triangular master vertex lies within the slave
!     surface
!
      implicit none
!
      logical in
!
      integer  inodesin(*),nnodesin,node,idi,nvertex,lvertex(*),i,
     &  inodesout(*),nnodesout,ido,nopes,j,mi(*)
!
      real*8 pvertex(3,*),pnodesin(3,*),slavstraight(20),xn(3),co(3,*),
     &  al,xl2(3,*),ratio(8),dist,xil,etl,vold(0:mi(2),*)
!
       in=.false.
!
      do
!
!        check whether nodes was already calatogued as being
!        inside the slave surface
!
         call nident(inodesin,node,nnodesin,idi)
         if(idi.gt.0) then
            if(inodesin(idi).eq.node) then
               in=.true.
               nvertex=nvertex+1
               do i=1,3
                  pvertex(i,nvertex)=pnodesin(i,idi)
               enddo
               lvertex(nvertex)=0
               exit
            endif
         endif
!
!        check whether nodes was already calatogued as being
!        outside the slave surface
!
         call nident(inodesout,node,nnodesout,ido)
         if(ido.gt.0) then
            if(inodesout(ido).eq.node) exit
         endif
!
!        node is not catalogued: check whether node is inside
!        or outside the slave surface
!
         do i=1,nopes
c            if((slavstraight(i*4-3)*co(1,node)+
c     &          slavstraight(i*4-2)*co(2,node)+
c     &          slavstraight(i*4-1)*co(3,node)+
            if((slavstraight(i*4-3)*(co(1,node)+vold(1,node))+
     &          slavstraight(i*4-2)*(co(2,node)+vold(2,node))+
     &          slavstraight(i*4-1)*(co(3,node)+vold(3,node))+
     &          slavstraight(i*4)).gt.0.d0) exit
            if(i.eq.nopes) in=.true.
         enddo
         if(in) then
            nvertex=nvertex+1
            lvertex(nvertex)=0
!
!           projecting the node on the mean slave plane
!
c            al=-xn(1)*co(1,node)-xn(2)*co(2,node)-xn(3)*co(3,node)-
c     &          slavstraight(nopes*4+4)
            al=-xn(1)*(co(1,node)+vold(1,node))-xn(2)*
     &   (co(2,node)+vold(2,node))-xn(3)*(co(3,node)+vold(3,node))-
     &          slavstraight(nopes*4+4)
            do i=1,3
               pvertex(i,nvertex)=co(i,node)+vold(i,node)+al*xn(i)
            enddo
!
!           projecting the node on the slave surface
!
            call attachline(xl2,pvertex,nopes,ratio,dist,xil,etl,xn)
!
!           cataloguein the node in inodesin
!
            nnodesin=nnodesin+1
            do j=nnodesin,idi+2,-1
               inodesin(j)=inodesin(j-1)
               do i=1,3
                  pnodesin(i,j)=pnodesin(i,j-1)
               enddo
            enddo
            inodesin(idi+1)=node
            do i=1,3
c               pnodesin(i,idi+1)=co(i,node)+vold(i,node)
               pnodesin(i,idi+1)=pvertex(i,nvertex)
            enddo
            exit
         else
!
!           cataloguein the node in inodesout
!
            nnodesout=nnodesout+1
            do j=nnodesout,ido+2,-1
               inodesout(j)=inodesout(j-1)
            enddo
            inodesout(ido+1)=node
            exit
         endif
      enddo
!
      return
      end
