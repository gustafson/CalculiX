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
      subroutine e_c3d_tlhs(co,nk,konl,lakonl,sm,nelem,ipvar,var)
!
!     computation of the velocity element matrix for the element with
!     the topology in konl
!
      implicit none
!
      character*8 lakonl
!
      integer konl(20),nk,nelem,i,j,ii,jj,kk,
     &  nope,mint3d,iflag,index,ipvar(*)
!
      real*8 co(3,*),xl(3,20),shp(4,20),xi,et,ze,xsj,
     &   sm(78,78),weight,var(*)
!
      include "gauss.f"
!
      data iflag /2/
!
      if(lakonl(4:4).eq.'2') then
         nope=20
      elseif(lakonl(4:4).eq.'8') then
         nope=8
      elseif(lakonl(4:5).eq.'10') then
         nope=10
      elseif(lakonl(4:4).eq.'4') then
         nope=4
      elseif(lakonl(4:5).eq.'15') then
         nope=15
      elseif(lakonl(4:4).eq.'6') then
         nope=6
      endif
!
      if(lakonl(4:5).eq.'8R') then
         mint3d=1
      elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) then
         if((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'S').or.
     &        (lakonl(7:7).eq.'E')) then
            mint3d=4
         else
            mint3d=8
         endif
      elseif(lakonl(4:4).eq.'2') then
         mint3d=27
      elseif(lakonl(4:5).eq.'10') then
         mint3d=4
      elseif(lakonl(4:4).eq.'4') then
         mint3d=1
      elseif(lakonl(4:5).eq.'15') then
         mint3d=9
      elseif(lakonl(4:5).eq.'6 ') then
         mint3d=2
      elseif(lakonl(4:5).eq.'6R') then
         mint3d=1
      endif
!
!     initialisation of sm
!
      do i=1,nope
         do j=1,nope
            sm(i,j)=0.d0
         enddo
      enddo
!
!     computation of the matrix: loop over the Gauss points
!
      index=ipvar(nelem)
      do kk=1,mint3d
         if(lakonl(4:5).eq.'8R') then
            weight=weight3d1(kk)
         elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) 
     &           then
            weight=weight3d2(kk)
         elseif(lakonl(4:4).eq.'2') then
            weight=weight3d3(kk)
         elseif(lakonl(4:5).eq.'10') then
            weight=weight3d5(kk)
         elseif(lakonl(4:4).eq.'4') then
            weight=weight3d4(kk)
         elseif(lakonl(4:5).eq.'15') then
            weight=weight3d8(kk)
         elseif(lakonl(4:5).eq.'6 ') then
            weight=weight3d7(kk)
         elseif(lakonl(4:5).eq.'6R') then
            weight=weight3d11(kk)
         endif
!
!        copying the shape functions, their derivatives and the
!        Jacobian determinant from field var
!
         do jj=1,nope
            do ii=1,4
               index=index+1
               shp(ii,jj)=var(index)
            enddo
         enddo
         index=index+1
         xsj=var(index)
         
         index=index+nope+14
!
         weight=weight*xsj
!     
         do jj=1,nope
!     
            do ii=1,jj
!     
!              lhs temperature and turbulence matrix
!     
               sm(ii,jj)=sm(ii,jj)
     &              +shp(4,ii)*shp(4,jj)*weight
            enddo
         enddo
      enddo
!
      return
      end

