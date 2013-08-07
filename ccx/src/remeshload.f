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
      subroutine remeshload(ipkon,kon,lakon,nelemload,sideload,iamload,
     &  xload,nload,ne,t1,iamt1,nam,ithermal,vold,mi,xloadold)
!
!     remeshing quadratic elements adjacent to contact surfaces
!     using appropriate linear elements (C3D8I,C3D4 and C3D6)
!
      implicit none
!
      character*8 lakon(*)
      character*20 sideload(*)
!
      integer ne,ipkon(*),kon(*),indexe,ifaceq(9,6),k,node,
     &  ifacew(8,5),i,j,ielface10(4,4),konl(27),indexer,nodes(8),
     &  ielface15(4,5),ielface20(4,6),id,ig,nelemload(2,*),
     &  iamload(2,*),nload,kon10(4,8),kon15(6,8),kon20(8,8),
     &  nelem,iamt1(*),nam,ithermal(2),iamplitude,mi(*),iload,
     &  nload_,nelement
!
      real*8 xload(2,*),t1(*),vold(0:mi(2),*),xloadold(2,*)
!
!     nodes belonging to the element faces
!
      data ifaceq /4,3,2,1,11,10,9,12,21,
     &            5,6,7,8,13,14,15,16,22,
     &            1,2,6,5,9,18,13,17,23,
     &            2,3,7,6,10,19,14,18,24,
     &            3,4,8,7,11,20,15,19,25,
     &            4,1,5,8,12,17,16,20,26/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
      data ielface10 /1,2,3,4,
     &                1,2,7,8,
     &                2,3,5,8,
     &                1,3,6,8/
      data ielface15 /1,2,3,4,
     &                5,6,7,8,
     &                1,2,5,6,
     &                2,3,6,7,
     &                1,3,5,7/
      data ielface20 /1,2,3,4,
     &                5,6,7,8,
     &                1,2,5,6,
     &                2,4,6,8,
     &                3,4,7,8,
     &                1,3,5,7/
      data kon10 /1,5,7,8,
     &            5,2,6,9,
     &            7,6,3,10,
     &            5,6,7,10,
     &            5,6,10,9,
     &            7,5,10,8,
     &            5,9,10,8,
     &            8,9,10,4/
      data kon15 /1,7,9,13,16,18,
     &            7,2,8,16,14,17,
     &            9,8,3,18,17,15,
     &            7,8,9,16,17,18,
     &            13,16,18,4,10,12,
     &            16,14,17,10,5,11,
     &            18,17,15,12,11,6,
     &            16,17,18,10,11,12/
      data kon20 /1,9,21,12,17,23,27,26,
     &            9,2,10,21,23,18,24,27,
     &            12,21,11,4,26,27,25,20,
     &            21,10,3,11,27,24,19,25,
     &            17,23,27,26,5,13,22,16,
     &            23,18,24,27,13,6,14,22,
     &            26,27,25,20,16,22,15,8,
     &            27,24,19,25,22,14,7,15/
!
      nload_=5*nload
!
      do i=1,ne
         indexe=ipkon(i)
         if(indexe.lt.-1) then
!
!           element has been remeshed
!
            indexe=-indexe-2
!
!           looking for any distributed load
!
            call nident2(nelemload,i,nload,id)
            do
               if((id.eq.0).or.(nelemload(1,id).ne.i)) exit
               read(sideload(id)(2:2),'(i1)') ig
!
!              kon(indexe+1) is the first subelement of element i
!
               nelem=kon(indexe+1)-1
               do j=1,4
c                  nload=nload+1
                  if(lakon(i)(4:4).eq.'2') then
                     nelement=nelem+ielface20(j,ig)
                  elseif(lakon(i)(4:5).eq.'10') then
                     nelement=nelem+ielface10(j,ig)
                  elseif(lakon(i)(4:5).eq.'15') then
                     nelement=nelem+ielface15(j,ig)
                  endif
                  call loadaddt(nelement,sideload(id),xload(1,id),
     &                 xload(2,id),nelemload,sideload,xload,nload,
     &                 nload_,iamload,iamload(2,id),iamload(1,id),
     &                 nam,nelemload(2,id),iload)
c                  nelemload(2,nload)=nelemload(2,id)
c                  sideload(nload)=sideload(id)
c                  if(nam.gt.0) then
c                     iamload(1,nload)=iamload(1,id)
c                     iamload(2,nload)=iamload(2,id)
c                  endif
c                  xload(1,nload)=xload(1,id)
c                  xload(2,nload)=xload(2,id)
                  xloadold(1,iload)=xloadold(1,id)
                  xloadold(2,iload)=xloadold(2,id)
               enddo
               id=id-1
            enddo
!
!     temperature interpolation in the extra nodes of the
!     remeshed elements
!
            if(ithermal(1).gt.0) then
               if(lakon(i)(4:4).eq.'2') then
                  nelem=kon(indexe+1)
                  indexer=ipkon(nelem)
!
!                 determining the nodes belonging to the element
!                 cave: C3D8I are 11-node elements!
!
                  konl(1)=kon(indexer+1)
                  do j=2,20
                     konl(j)=kon(indexe+j)
                  enddo
                  konl(21)=kon(indexer+3)
                  konl(22)=kon(indexer+51)
                  konl(23)=kon(indexer+6)
                  konl(24)=kon(indexer+18)
                  konl(25)=kon(indexer+29)
                  konl(26)=kon(indexer+8)
                  konl(27)=kon(indexer+7)
!
!                 check whether a temperature amplitude applies
!
                  if(nam.gt.0) then
                     iamplitude=iamt1(konl(1))
                     do j=2,20
                        if(iamt1(konl(j)).ne.iamplitude) then
                           write(*,*) '*ERROR in remeshload:'
                           write(*,*) 
     &'different temperature amplitudes are applied to different nodes'
                           write(*,*) 
     &'of one and the same quadratic element which needs to be remeshed'
                           write(*,*) 
     &'because of contact. This is not allowed. Element nr:',i
                           stop
                        endif
                     enddo
                     do j=21,27
                        iamt1(konl(j))=iamt1(konl(1))
                     enddo
                  endif
!
                  if(ithermal(1).eq.1) then
!
!                    facial nodes
!
                     do j=1,6
                        node=konl(20+j)
                        do k=1,8
                           nodes(k)=konl(ifaceq(k,j))
                        enddo
                        t1(node)=0.d0
                        do k=1,4
                           t1(node)=t1(node)-t1(nodes(k))
                        enddo
                        do k=5,8
                           t1(node)=t1(node)+2.d0*t1(nodes(k))
                        enddo
                        t1(node)=t1(node)/4.d0
                     enddo
!
!                    volumetric node
!                     
                     node=konl(27)
                     t1(node)=0.d0
                     do k=1,8
                        t1(node)=t1(node)-t1(konl(k))
                     enddo
                     do k=9,20
                        t1(node)=t1(node)+t1(konl(k))
                     enddo
                     t1(node)=t1(node)/4.d0
!
!                    dummy C3D8I-nodes: dummy temperature
!                    (else unitialized value in tempload.f)
!
                     do k=0,7
                        do j=9,11
                           node=kon(indexer+11*k+j)
                           t1(node)=0.d0
                        enddo
                     enddo
!
                  endif
               elseif(lakon(i)(4:5).eq.'15') then
                  nelem=kon(indexe+1)
                  indexer=ipkon(nelem)
!
!                 determining the nodes belonging to the element
!
                  konl(1)=kon(indexer+1)
                  do j=2,15
                     konl(j)=kon(indexe+j)
                  enddo
                  konl(16)=kon(indexer+5)
                  konl(17)=kon(indexer+12)
                  konl(18)=kon(indexer+6)
!
!                 check whether a temperature amplitude applies
!
                  if(nam.gt.0) then
                     iamplitude=iamt1(konl(1))
                     do j=2,15
                        if(iamt1(konl(j)).ne.iamplitude) then
                           write(*,*) '*ERROR in remeshload:'
                           write(*,*) 
     &'different temperature amplitudes are applied to different nodes'
                           write(*,*) 
     &'of one and the same quadratic element which needs to be remeshed'
                           write(*,*) 
     &'because of contact. This is not allowed. Element nr:',i
                           stop
                        endif
                     enddo
                     do j=16,18
                        iamt1(konl(j))=iamt1(konl(1))
                     enddo
                  endif
!
                  if(ithermal(1).eq.1) then
!
!                    facial nodes
!
                     do j=3,5
                        node=konl(13+j)
                        do k=1,8
                           nodes(k)=konl(ifacew(k,j))
                        enddo
                        t1(node)=0.d0
                        do k=1,4
                           t1(node)=t1(node)-t1(nodes(k))
                        enddo
                        do k=5,8
                           t1(node)=t1(node)+2.d0*t1(nodes(k))
                        enddo
                     enddo
                     t1(node)=t1(node)/4.d0
                  endif
               endif
            endif
!
         endif
      enddo
!     
      return
      end
