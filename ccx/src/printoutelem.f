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
      subroutine printoutelem(prlab,ipkon,lakon,kon,co,
     &     ener,mint_,ii,nelem,energytot,volumetot)
!
!     stores whole element results for element "nelem" in the .dat file
!
      implicit none
!
      character*6 prlab(*)
      character*8 lakon(*)
!
      integer ipkon(*),nelem,ii,kon(*),mint_,nope,indexe,j,k,konl(20),
     &  mint3d,jj,nener,iflag
!
      real*8 ener(mint_,*),energytot,volumetot,energy,volume,co(3,*),
     &  xl(3,20),xi,et,ze,xsj,shp(4,20),weight
!
      include "gauss.f"
!
      data iflag /2/
!
      if(ipkon(nelem).lt.0) return
!
      if(prlab(ii)(1:4).eq.'ELSE') then
         nener=1
      else
         nener=0
      endif
!
      indexe=ipkon(nelem)
!
      if(lakon(nelem)(4:4).eq.'2') then
         nope=20
      elseif(lakon(nelem)(4:4).eq.'8') then
         nope=8
      elseif(lakon(nelem)(4:5).eq.'10') then
         nope=10
      elseif(lakon(nelem)(4:4).eq.'4') then
         nope=4
      elseif(lakon(nelem)(4:5).eq.'15') then
         nope=15
      else
         nope=6
      endif
!
      do j=1,nope
         konl(j)=kon(indexe+j)
         do k=1,3
            xl(k,j)=co(k,konl(j))
         enddo
      enddo
!
      energy=0.d0
      volume=0.d0
!
      if(lakon(nelem)(4:5).eq.'8R') then
         mint3d=1
      elseif((lakon(nelem)(4:4).eq.'8').or.
     &        (lakon(nelem)(4:6).eq.'20R')) then
         mint3d=8
      elseif(lakon(nelem)(4:4).eq.'2') then
         mint3d=27
      elseif(lakon(nelem)(4:5).eq.'10') then
         mint3d=4
      elseif(lakon(nelem)(4:4).eq.'4') then
         mint3d=1
      elseif(lakon(nelem)(4:5).eq.'15') then
         mint3d=9
      else
         mint3d=2
      endif
!
      do jj=1,mint3d
         if(lakon(nelem)(4:5).eq.'8R') then
            xi=gauss3d1(1,jj)
            et=gauss3d1(2,jj)
            ze=gauss3d1(3,jj)
            weight=weight3d1(jj)
         elseif((lakon(nelem)(4:4).eq.'8').or.
     &           (lakon(nelem)(4:6).eq.'20R'))
     &           then
            xi=gauss3d2(1,jj)
            et=gauss3d2(2,jj)
            ze=gauss3d2(3,jj)
            weight=weight3d2(jj)
         elseif(lakon(nelem)(4:4).eq.'2') then
            xi=gauss3d3(1,jj)
            et=gauss3d3(2,jj)
            ze=gauss3d3(3,jj)
            weight=weight3d3(jj)
         elseif(lakon(nelem)(4:5).eq.'10') then
            xi=gauss3d5(1,jj)
            et=gauss3d5(2,jj)
            ze=gauss3d5(3,jj)
            weight=weight3d5(jj)
         elseif(lakon(nelem)(4:4).eq.'4') then
            xi=gauss3d4(1,jj)
            et=gauss3d4(2,jj)
            ze=gauss3d4(3,jj)
            weight=weight3d4(jj)
         elseif(lakon(nelem)(4:5).eq.'15') then
            xi=gauss3d8(1,jj)
            et=gauss3d8(2,jj)
            ze=gauss3d8(3,jj)
            weight=weight3d8(jj)
         else
            xi=gauss3d7(1,jj)
            et=gauss3d7(2,jj)
            ze=gauss3d7(3,jj)
            weight=weight3d7(jj)
         endif
!
         if(nope.eq.20) then
            call shape20h(xi,et,ze,xl,xsj,shp,iflag)
         elseif(nope.eq.8) then
            call shape8h(xi,et,ze,xl,xsj,shp,iflag)
         elseif(nope.eq.10) then
            call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
         elseif(nope.eq.4) then
            call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
         elseif(nope.eq.15) then
            call shape15w(xi,et,ze,xl,xsj,shp,iflag)
         else
            call shape6w(xi,et,ze,xl,xsj,shp,iflag)
         endif
!
         if(nener.eq.1) energy=energy+weight*xsj*ener(jj,nelem)
         volume=volume+weight*xsj
      enddo
!
      volumetot=volumetot+volume
      energytot=energytot+energy
!
!     writing to file
!
      if((prlab(ii)(1:5).eq.'ELSE ').or.
     &   (prlab(ii)(1:5).eq.'ELSET')) then
         write(5,'(i6,1p,1x,e11.4)') nelem,energy
      elseif((prlab(ii)(1:5).eq.'EVOL ').or.
     &       (prlab(ii)(1:5).eq.'EVOLT')) then
         write(5,'(i6,1p,1x,e11.4)') nelem,volume
      endif
!
      return
      end
