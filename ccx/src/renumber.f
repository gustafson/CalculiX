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
      subroutine renumber(nk,kon,ipkon,lakon,ne,ipompc,nodempc,nmpc,nnn,
     &  npn,adj,xadj,iw,mmm,xnpn,inum1,inum2)
!
!     renumbers the nodes to reduce the profile length
!
      implicit none
!
      character*8 lakon(*)
!
      integer kon(*),ipompc(*),nodempc(3,*),npn(*),inum1(*),inum2(*),
     &  nnn(*),iw(*),mmm(*),xnpn(*),adj(*),xadj(*),ipkon(*),node
!
      integer nne,inpn,iadj,nk,ne,nmpc,i,j,nterm,e2,oldpro,newpro,
     &  index,kflag,nope,indexe,oldpro_exp,newpro_exp,nknew
!
      kflag=2
      nne=0
      inpn=0
      iadj=0
!
!     taking the elements into account
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         indexe=ipkon(i)
         if(lakon(i)(4:4).eq.'2') then
            nope=20
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:5).eq.'15') then
            nope=15
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         elseif(lakon(i)(1:1).eq.'E') then
            read(lakon(i)(8:8),'(i1)') nope
            nope=nope+1
         elseif(lakon(i)(1:1).eq.'D') then
            cycle
         endif
!
         nne=nne+1
         xnpn(nne)=inpn+1
         do j=1,nope
           node=kon(indexe+j)
           npn(inpn+j)=node
           inum1(node)=1
         enddo
         inpn=inpn+nope
         iadj=iadj+nope*(nope-1)
      enddo
!
!     taking the equations into account
!
      do i=1,nmpc
         nne=nne+1
         xnpn(nne)=inpn+1
         index=ipompc(i)
         nterm=0
         do
            nterm=nterm+1
            node=nodempc(1,index)
            npn(inpn+nterm)=node
            inum1(node)=1
            index=nodempc(3,index)
            if(index.eq.0) exit
         enddo
         inpn=inpn+nterm
         iadj=iadj+nterm*(nterm-1)
      enddo
!
      xnpn(nne+1)=inpn+1
!
!     numbering the node which are really used and changing the
!     numbers in npn
!
      nknew=0
      do i=1,nk
         if(inum1(i).gt.0) then
            nknew=nknew+1
            inum1(i)=nknew
         endif
      enddo
      do i=1,inpn
         npn(i)=inum1(npn(i))
      enddo
!
      call graph(nknew,nne,inpn,npn,xnpn,iadj,adj,xadj)
!
      e2=xadj(nknew+1)-1
!
      call label(nknew,e2,adj,xadj,mmm,iw,oldpro,newpro,oldpro_exp,
     &  newpro_exp)
!
      write(*,*) 'old profile = ',oldpro_exp,'*2147483647+',oldpro
      write(*,*) 'new profile = ',newpro_exp,'*2147483647+',newpro
      write(*,*)
!
!     restoring the original numbering
!
      do i=1,nk
         if(inum1(i).ne.0) then
            inum2(inum1(i))=i
         endif
      enddo
      index=0
      do i=1,nk
         if(inum1(i).eq.0) then
            inum1(i)=i
         else
            index=index+1
            inum1(i)=inum2(mmm(index))
         endif
      enddo
!
      call isortii(inum1,nnn,nk,kflag)
!
      return
      end
