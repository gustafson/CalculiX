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
      subroutine map3dto1d2d(yn,ipkon,inum,kon,lakon,nfield,nk,
     &  ne,cflag,co)
!
!     interpolates 3d field nodal values to 1d/2d nodal locations
!
!     the number of internal state variables is limited to 999
!     (cfr. array field)
!
      implicit none
!
      character*1 cflag
      character*8 lakon(*),lakonl
!
      integer ipkon(*),inum(*),kon(*),ne,indexe,nfield,nk,i,j,k,l,m,
     &  node3(4,3),node6(2,6),node8(2,8),node2d,node3d,indexe2d,ne1d2d,
     &  node3m(8,3),node(8),m1,m2,nodea,nodeb,nodec,iflag
!
      real*8 yn(nfield,*),cg(3),p(3),pcg(3),t(3),xl(0:3,8),shp(4,8),
     &  xsj(3),e1(3),e2(3),e3(3),s(6),dd,xi,et,co(3,*),xs(3,2)
!
      real*8 gauss2d1(2,1),gauss2d2(2,4),gauss2d3(2,9),gauss2d4(2,1),
     &  gauss2d5(2,3),gauss3d1(3,1),gauss3d2(3,8),gauss3d3(3,27),
     &  gauss3d4(3,1),gauss3d5(3,4),gauss3d6(3,15),gauss3d7(3,2),
     &  gauss3d8(3,9),gauss3d9(3,18),weight2d1(1),weight2d2(4),
     &  weight2d3(9),weight2d4(1),weight2d5(3),weight3d1(1),
     &  weight3d2(8),weight3d3(27),weight3d4(1),weight3d5(4),
     &  weight3d6(15),weight3d7(2),weight3d8(9),weight3d9(18)
!
      include "gauss.f"
!
      data node3 /1,4,8,5,9,11,15,13,2,3,7,6/
      data node3m /1,5,8,4,17,16,20,12,
     &             0,0,0,0,0,0,0,0,
     &             3,7,6,2,19,14,18,10/
      data node6 /1,4,2,5,3,6,7,10,8,11,9,12/
      data node8 /1,5,2,6,3,7,4,8,9,13,10,14,11,15,12,16/
      data iflag /2/
!
!     removing any results in 1d/2d nodes
!
      ne1d2d=0
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         lakonl=lakon(i)
         if((lakonl(7:7).eq.' ').or.(lakonl(7:7).eq.'G').or.
     &      (lakonl(1:1).ne.'C')) cycle
         ne1d2d=1
         indexe=ipkon(i)
!
!        inactivating the 3d expansion nodes of 1d/2d elements
!
         do j=1,20
            inum(kon(indexe+j))=0
         enddo
!
         if(lakonl(4:5).eq.'15') then
            indexe2d=indexe+15
            do j=1,6
               node2d=kon(indexe2d+j)
               inum(node2d)=0
               do k=1,nfield
                  yn(k,node2d)=0.d0
               enddo
            enddo
         elseif(lakonl(7:7).eq.'B') then
            indexe2d=indexe+20
            do j=1,3
               node2d=kon(indexe2d+j)
               inum(node2d)=0
               do k=1,nfield
                  yn(k,node2d)=0.d0
               enddo
            enddo
         else
            indexe2d=indexe+20
            do j=1,8
               node2d=kon(indexe2d+j)
               inum(node2d)=0
               do k=1,nfield
                  yn(k,node2d)=0.d0
               enddo
            enddo
         endif
!
      enddo
!
!     if no 1d/2d elements return
!
      if(ne1d2d.eq.0) return
!
!     interpolation of 3d results on 1d/2d nodes
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         lakonl=lakon(i)
         if((lakonl(7:7).eq.' ').or.(lakonl(7:7).eq.'G').or.
     &      (lakonl(1:1).ne.'C')) cycle
         indexe=ipkon(i)
!
         if(lakonl(4:5).eq.'15') then
            indexe2d=indexe+15
            do j=1,6
               node2d=kon(indexe2d+j)
               do l=1,2
                  inum(node2d)=inum(node2d)-1
                  node3d=kon(indexe+node6(l,j))
                  do k=1,nfield
                     yn(k,node2d)=yn(k,node2d)+yn(k,node3d)
                  enddo
               enddo
            enddo
         elseif(lakonl(7:7).eq.'B') then
            indexe2d=indexe+20
            if(cflag.ne.'M') then
!
!              mean values for beam elements
!
               do j=1,3
                  node2d=kon(indexe2d+j)
                  do l=1,4
                     inum(node2d)=inum(node2d)-1
                     node3d=kon(indexe+node3(l,j))
                     do k=1,nfield
                        yn(k,node2d)=yn(k,node2d)+yn(k,node3d)
                     enddo
                  enddo
               enddo
            else
!
!              section forces for beam elements
!
               do j=1,3,2
                  node2d=kon(indexe2d+j)
                  inum(node2d)=inum(node2d)-1
!
!                 coordinates of the nodes belonging to the section
!
                  do l=1,8
                     node(l)=kon(indexe+node3m(l,j))
                     do m=1,3
                        xl(m,l)=co(m,node(l))
                     enddo
                  enddo
!
!                 center of gravity and unit vectors 1 and 2
!
                  do m=1,3
                     cg(m)=(xl(m,6)+xl(m,8))/2.d0
                     if(j.eq.1) then
                        e1(m)=(xl(m,8)-xl(m,6))
                     else
                        e1(m)=(xl(m,6)-xl(m,8))
                     endif
                     e2(m)=(xl(m,7)-xl(m,5))
                  enddo
                  dd=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
                  do m=1,3
                     e1(m)=e1(m)/dd
                  enddo
                  dd=dsqrt(e2(1)*e2(1)+e2(2)*e2(2)+e2(3)*e2(3))
                  do m=1,3
                     e2(m)=e2(m)/dd
                  enddo
!
!                 loop over the integration points (2x2)
!                  
                  do l=1,4
c                  do l=1,9
                     xi=gauss2d2(1,l)
                     et=gauss2d2(2,l)
c                     xi=gauss2d3(1,l)
c                     et=gauss2d3(2,l)
                     call shape8q(xi,et,xl,xsj,xs,shp,iflag)
!
!                    local unit normal (only once per section)
!
                     if(l.eq.1) then
                        dd=dsqrt(xsj(1)*xsj(1)+xsj(2)*xsj(2)+
     &                        xsj(3)*xsj(3))
                        do m=1,3
                           e3(m)=xsj(m)/dd
                        enddo
                     endif
!
!                    local stress tensor
!
                     do m1=1,6
                        s(m1)=0.d0
                        do m2=1,8
                           s(m1)=s(m1)+shp(4,m2)*yn(m1,node(m2))
                        enddo
                     enddo
!
!                    local coordinates
!
                     do m1=1,3
                        p(m1)=0.d0
                        do m2=1,8
                           p(m1)=p(m1)+shp(4,m2)*xl(m1,m2)
                        enddo
                        pcg(m1)=p(m1)-cg(m1)
                     enddo
!
!                    local stress vector on section
!
                     t(1)=s(1)*xsj(1)+s(4)*xsj(2)+s(5)*xsj(3)
                     t(2)=s(4)*xsj(1)+s(2)*xsj(2)+s(6)*xsj(3)
                     t(3)=s(5)*xsj(1)+s(6)*xsj(2)+s(3)*xsj(3)
!
!                    section forces
!
                     yn(1,node2d)=yn(1,node2d)+
     &                   (e1(1)*t(1)+e1(2)*t(2)+e1(3)*t(3))
c     &                   *weight2d3(l)
                     yn(2,node2d)=yn(2,node2d)+
     &                   (e2(1)*t(1)+e2(2)*t(2)+e2(3)*t(3))
c     &                   *weight2d3(l)
                     yn(3,node2d)=yn(3,node2d)+
     &                   (e3(1)*t(1)+e3(2)*t(2)+e3(3)*t(3))
c     &                   *weight2d3(l)
!
!                    section moments
!
                     yn(4,node2d)=yn(4,node2d)+
     &                     (e1(1)*pcg(2)*t(3)+e1(2)*pcg(3)*t(1)+
     &                     e1(3)*pcg(1)*t(2)-e1(3)*pcg(2)*t(1)-
     &                     e1(1)*pcg(3)*t(2)-e1(2)*pcg(1)*t(3))
c     &                     *weight2d3(l)
                     yn(5,node2d)=yn(5,node2d)+
     &                     (e2(1)*pcg(2)*t(3)+e2(2)*pcg(3)*t(1)+
     &                     e2(3)*pcg(1)*t(2)-e2(3)*pcg(2)*t(1)-
     &                     e2(1)*pcg(3)*t(2)-e2(2)*pcg(1)*t(3))
c     &                     *weight2d3(l)
                     yn(6,node2d)=yn(6,node2d)+
     &                     (e3(1)*pcg(2)*t(3)+e3(2)*pcg(3)*t(1)+
     &                     e3(3)*pcg(1)*t(2)-e3(3)*pcg(2)*t(1)-
     &                     e3(1)*pcg(3)*t(2)-e3(2)*pcg(1)*t(3))
c     &                     *weight2d3(l)
!
                  enddo
               enddo
!
            endif
         else
            indexe2d=indexe+20
            do j=1,8
               node2d=kon(indexe2d+j)
               do l=1,2
                  inum(node2d)=inum(node2d)-1
                  node3d=kon(indexe+node8(l,j))
                  do k=1,nfield
                     yn(k,node2d)=yn(k,node2d)+yn(k,node3d)
                  enddo
c                  write(*,*) node2d,node3d,yn(4,node2d),yn(4,node3d)
               enddo
            enddo
         endif
!
      enddo
!
!     taking the mean of nodal contributions coming from different
!     elements having the node in common
!
      do i=1,nk
         if(inum(i).lt.0) then
            inum(i)=-inum(i)
            do j=1,nfield
               yn(j,i)=yn(j,i)/inum(i)
            enddo
         endif
      enddo
!
!     beam section forces in the middle nodes
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         lakonl=lakon(i)
         if((lakonl(7:7).eq.' ').or.(lakonl(7:7).eq.'G').or.
     &      (lakonl(1:1).ne.'C')) cycle
         indexe=ipkon(i)
!
         if(lakonl(7:7).eq.'B') then
            indexe2d=indexe+20
            if(cflag.eq.'M') then
!
!              section forces in the middle node are the mean
!              of those in the end nodes
!
               nodea=kon(indexe2d+1)
               nodeb=kon(indexe2d+2)
               nodec=kon(indexe2d+3)
               inum(nodeb)=1
               do j=1,6
                  yn(j,nodeb)=yn(j,nodeb)+(yn(j,nodea)+yn(j,nodec))/2.d0
               enddo
!
            endif
         endif
!
      enddo
!
      return
      end
