!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2014 Guido Dhondt
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
      subroutine precfd(ne,ipkon,kon,lakon,ipnei,neifa,neiel,ipoface,
     &  nodface,ielfa,nkonnei,nface,ifaext,nfaext,
     &  isolidsurf,nsolidsurf,set,nset,istartset,iendset,ialset,
     &  vel,vfa,vold,mi)
!
!
      implicit none
!
      character*8 lakon(*)
      character*81 set(*),noset
!
      integer ne,ipkon(*),ipnei(*),ipoface(*),nodface(5,*),neifa(*),
     &  ielfa(4,*),nkonnei,nface,i,j,k,index,indexe,iface,neiel(*),
     &  nfaext,ifaext(*),isolidsurf(*),nsolidsurf,indexf,jopposite(6),
     &  nset,istartset(*),iendset(*),ialset(*),iaux,kflag,ifour,
     &  ifaceq(8,6),ifacet(7,4),ifacew(8,5),kon(*),nodes(4),iel1,
     &  indexold,ifree,ifreenew,ifreenei,iloc1,mi(*)
!
      real*8 vel(0:4,*),vfa(0:4,*),vold(0:mi(2),*)
!
!     nodes belonging to the element faces
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,11,
     &             1,2,4,5,9,8,12,
     &             2,3,4,6,10,9,13,
     &             1,4,3,8,10,7,14/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
      data jopposite /2,1,5,6,3,4/
!
      kflag=1
      ifour=4
!
!     determining the external element faces of the fluid mesh 
!     the faces are catalogued by the three lowes nodes numbers
!     in ascending order. ipoface(i) points to a face for which
!     node i is the lowest node and nodface(1,ipoface(i)) and
!     nodface(2,ipoface(i)) are the next lower ones. 
!     nodface(3,ipoface(i)) contains the element number,
!     nodface(4,ipoface(i)) the face number and nodface(5,ipoface(i))
!     is a pointer to the next surface for which node i is the
!     lowest node; if there are no more such surfaces the pointer
!     has the value zero
!     An external element face is one which belongs to one element
!     only
!
      ifree=1
      ifreenei=0
!
      do i=1,6*ne-1
         nodface(5,i)=i+1
      enddo
      do i=1,ne
         if(lakon(i)(1:1).ne.'F') cycle
         indexe=ipkon(i)
         if(lakon(i)(4:4).eq.'8') then
            ipnei(i)=ifreenei
            do j=1,6
               do k=1,4
                  nodes(k)=kon(indexe+ifaceq(k,j))
               enddo
               call isortii(nodes,iaux,ifour,kflag)
               indexold=0
               index=ipoface(nodes(1))
               do
!
!                 adding a surface which has not been 
!                 catalogued so far
!
                  if(index.eq.0) then
                     ifreenew=nodface(5,ifree)
                     nodface(1,ifree)=nodes(2)
                     nodface(2,ifree)=nodes(3)
                     nodface(3,ifree)=i
                     nodface(4,ifree)=j
                     nodface(5,ifree)=ipoface(nodes(1))
                     ipoface(nodes(1))=ifree
                     ifreenei=ifreenei+1
                     neifa(ifreenei)=ifree
                     ielfa(1,ifree)=i
                     ielfa(4,ifree)=j
                     ifree=ifreenew
                     exit
                  endif
!
!                 removing a surface which has already
!                 been catalogued
!
                  if((nodface(1,index).eq.nodes(2)).and.
     &               (nodface(2,index).eq.nodes(3))) then
                     ifreenei=ifreenei+1
                     neifa(ifreenei)=index
                     ielfa(2,index)=i
                  endif
                  indexold=index
                  index=nodface(5,index)
               enddo
            enddo
         endif
      enddo
!
      nkonnei=ifreenei
      nface=ifree
!
!     determining neighboring elements of element i
!
      do i=1,ne
         if(lakon(i)(1:1).ne.'F') cycle
         index=ipnei(i)
         if(lakon(i)(4:4).eq.'8') then
            do j=1,6
               iface=neifa(index+j)
               if(ielfa(1,iface).eq.i) then
                  neiel(index+j)=ielfa(2,iface)
               else
                  neiel(index+j)=ielfa(1,iface)
               endif
            enddo
         endif
      enddo
!
!     catalogueing external faces
!
      nfaext=0
!
      do i=1,nface
         if(ielfa(2,i).ne.0) cycle
         nfaext=nfaext+1
         ifaext(nfaext)=i
         j=ielfa(4,i)
         iel1=ielfa(1,i)
         indexf=ipnei(iel1)
         j=jopposite(j)
         ielfa(3,i)=neiel(indexf+j)
      enddo
!
!     faces belonging to solid surfaces
!
      noset(1:13)='SOLIDSURFACET'
      do i=1,nset
         if(set(i)(1:13).eq.noset(1:13)) exit
      enddo
      if(i.gt.nset) then
         write(*,*) '*WARNING in precfd: node set SOLID SURFACE '
         write(*,*) '         has not been defined.'
      else
         do j=istartset(i),iendset(i)
            nsolidsurf=nsolidsurf+1
            isolidsurf(nsolidsurf)=ialset(j)
         enddo
         call isortii(isolidsurf,iaux,nsolidsurf,kflag)
      endif
!
!     initial conditions: element values
!     vel was initialized at allocation
!
      do i=1,ne
         if(lakon(i)(1:1).ne.'F') cycle
         indexe=ipkon(i)
         if(lakon(i)(4:4).eq.'8') then
            do j=0,4
               do k=1,8
                  vel(j,i)=vel(j,i)+vold(j,kon(indexe+k))
               enddo
               vel(j,i)=vel(j,i)/8.d0
            enddo
         endif
      enddo
!
!     initial conditions: facial values
!     vfa was initialized at allocation
!
      do i=1,nface
         indexe=ipkon(ielfa(1,i))
         iloc1=ielfa(4,i)
         if(lakon(i)(4:4).eq.'8') then
            do j=1,4
               nodes(j)=kon(indexe+ifaceq(j,iloc1))
               do k=0,4
                  vfa(k,j)=vfa(k,j)+vold(k,nodes(j))
               enddo
            enddo
            vfa(k,j)=vfa(k,j)/4.d0
         endif
      enddo
!
      return
      end
