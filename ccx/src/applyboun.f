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
      subroutine applyboun(ifaext,nfaext,ielfa,ikboun,ilboun,
     &  nboun,typeboun,nelemload,nload,sideload,isolidsurf,nsolidsurf,
     &  ifabou,nfabou)
!
!     stores pointers to ifabou in ielfa(2,*) at those locations
!     which are zero (external faces)
!     stores pointers to the boundary conditions in field ifabou
!
      implicit none
!
      character*1 typeboun(*)
      character*20 sideload(*)
!
      integer nfabou,ifaext(*),nfaext,ielem,ielfa(4,*),iface,jface,
     &  j,ikbounfa,ikboun(*),nboun,id,ilboun(*),iboun,nelemload(2,*),
     &  nload,isolidsurf(*),nsolidsurf,ifabou(*),i
!
      nfabou=1
!
      do i=1,nfaext
!
!        number of the face in field ielfa
!
         iface=ifaext(i)
!
!        adjacent element number
!
         ielem=ielfa(1,iface)
!
!        face label used to apply the SPC
!
         jface=10*ielem+ielfa(4,iface)
!
!        loop over the degrees of freedom
!
         do j=0,4
            ikbounfa=10*jface+j
            call nident(ikboun,ikbounfa,nboun,id)
            if(id.gt.0) then
               if(ikboun(id).eq.ikbounfa) then
                  iboun=ilboun(id)
                  if(typeboun(iboun).ne.'F') cycle
                  if(ielfa(2,iface).eq.0) then
                     ielfa(2,iface)=-nfabou
                     nfabou=nfabou+7
                  endif
                  ifabou(-ielfa(2,iface)+j)=iboun
               endif
            endif
         enddo
!
!        heat flux
!
         call nident2(nelemload,ielem,nload,id)
!
         do
            if(id.gt.0) then
               if(nelemload(1,id).eq.ielem) then
                  if(sideload(id)(1:1).eq.'S') then
                     if(ielfa(2,iface).eq.0) then
                        ielfa(2,iface)=-nfabou
                        nfabou=nfabou+7
                     endif
                     ifabou(-ielfa(2,iface)+6)=id
                  endif
                  id=id-1
                  cycle
               else
                  exit
               endif
            endif
         enddo
!
!        wall
!
         call nident(isolidsurf,jface,nsolidsurf,id)
         if(id.gt.0) then
            if(isolidsurf(id).eq.jface) then
               if(ielfa(2,iface).eq.0) then
                  ielfa(2,iface)=-nfabou
                  nfabou=nfabou+7
               endif
               ifabou(-ielfa(2,iface)+5)=1
            endif
         endif
      enddo
!
!     dimension of field ifabou containing the pointers to the
!     boundary conditions
!
      nfabou=nfabou-1
!     
      return
      end
