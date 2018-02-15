!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2017 Guido Dhondt
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
      subroutine createialdesi(ndesi,nodedesi,iponoel,inoel,istartdesi,
     &                         ialdesi,lakon,ipkon,kon,nodedesiinv,
     &                         icoordinate,noregion)
!
      implicit none
!
      character*8 lakon(*)
!
      integer ndesi,node,nodedesi(*),iponoel(*),inoel(2,*),
     &   istartdesi(*),ialdesi(*),ifree,index,i,ipkon(*),kon(*),
     &   nodedesiinv(*),icoordinate,indexe,nopedesi,nnodes,nelem,
     &   m,nope,noregion
!
!     determining the elements belonging to a given design
!     variable i. They are stored in ialdesi(istartdesi(i))..
!     ...up to..... ialdesi(istartdesi(i+1)-1)
!
      ifree=1
!
      if(icoordinate.eq.1) then
!
!        coordinates as design variables
!
!        an element is taken into account if more than nopedesign
!        nodes in the element are design variables (important for
!        design nodes on the border of the design domain)
!
         do i=1,ndesi
            istartdesi(i)=ifree
            node=nodedesi(i)
            index=iponoel(node)
            do
               if(index.eq.0) exit
               nelem=inoel(1,index)
!
               if(lakon(nelem)(4:4).eq.'8') then
                  nopedesi=3
                  nope=8
               elseif(lakon(nelem)(4:5).eq.'20') then
                  nopedesi=5
                  nope=20
               elseif(lakon(nelem)(4:5).eq.'10') then
                  nopedesi=3
                  nope=10
               elseif(lakon(nelem)(4:4).eq.'4') then
                  nopedesi=3
                  nope=4
               elseif(lakon(nelem)(4:4).eq.'6') then
                  nopedesi=3
                  nope=6
               elseif(lakon(nelem)(4:5).eq.'15') then
                  nopedesi=3
                  nope=15
               endif
               if(noregion.eq.1) nopedesi=0
!
               indexe=ipkon(nelem)
!
!              summing the design variables in the element
!
               nnodes=0
               do m=1,nope
                  if(nodedesiinv(kon(indexe+m)).eq.1) then
                     nnodes=nnodes+1
                  endif
               enddo
!
               if(nnodes.ge.nopedesi) then
                  ialdesi(ifree)=nelem
                  ifree=ifree+1
               endif
               index=inoel(2,index)
            enddo
         enddo
         istartdesi(ndesi+1)=ifree
      else
!         
!        orientation as design variables
!
         do i=1,ndesi
            istartdesi(i)=ifree
            node=nodedesi(i)
            index=iponoel(node)
            do
               if(index.eq.0) exit
               ialdesi(ifree)=inoel(1,index)
               ifree=ifree+1
               index=inoel(2,index)
            enddo
         enddo
         istartdesi(ndesi+1)=ifree
      endif
!
      return
      end
