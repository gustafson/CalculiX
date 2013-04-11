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
      subroutine printoutint(prlab,ipkon,lakon,stx,eme,xstate,ener,
     &  mint_,nstate_,l1,lb,ii,nelem,qfx,orab,ielorien,norien,co,kon)
!
!     stores integration point results for element "nelem" in the .dat file
!
      implicit none
!
      character*6 prlab(*)
      character*8 lakon(*)
!
      integer ipkon(*),mint_,nstate_,nelem,l,lb,ii,mint3d,j,k,nope,
     &  ielorien(*),norien,kon(*),konl,indexe,m,iorien,iflag,l1
!
      real*8 stx(6,mint_,*),eme(6,mint_,*),xstate(nstate_,mint_,*),
     &  ener(mint_,*),qfx(3,mint_,*),xi,et,ze,xl(3,20),xsj,shp(4,20),
     &  coords(3,27),weight,orab(7,*),co(3,*),a(3,3),b(3,3),c(3,3),
     &  qfxl(3)
!
      include "gauss.f"
!
      data iflag /1/
!
      if(ipkon(nelem).lt.0) return
!
!     check whether transformation is necessary (if orientation
!     is applied and output in local system is requested)
!
      if((norien.eq.0).or.(prlab(ii)(6:6).eq.'G')) then
         iorien=0
      else
         iorien=ielorien(nelem)
      endif
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
      elseif(lakon(nelem)(4:4).eq.'6') then
         mint3d=2
      else
         return
      endif
!
!     calculation of the integration point coordinates for
!     output in the local system
!
      if(iorien.ne.0) then
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
         elseif(lakon(nelem)(4:4).eq.'6') then
            nope=6
         endif
!
         indexe=ipkon(nelem)
         do j=1,nope
            konl=kon(indexe+j)
            do k=1,3
               xl(k,j)=co(k,konl)
            enddo
         enddo
!
         do j=1,mint3d
            if(lakon(nelem)(4:5).eq.'8R') then
               xi=gauss3d1(1,j)
               et=gauss3d1(2,j)
               ze=gauss3d1(3,j)
               weight=weight3d1(j)
            elseif((lakon(nelem)(4:4).eq.'8').or.
     &              (lakon(nelem)(4:6).eq.'20R'))
     &              then
               xi=gauss3d2(1,j)
               et=gauss3d2(2,j)
               ze=gauss3d2(3,j)
               weight=weight3d2(j)
            elseif(lakon(nelem)(4:4).eq.'2') then
               xi=gauss3d3(1,j)
               et=gauss3d3(2,j)
               ze=gauss3d3(3,j)
               weight=weight3d3(j)
            elseif(lakon(nelem)(4:5).eq.'10') then
               xi=gauss3d5(1,j)
               et=gauss3d5(2,j)
               ze=gauss3d5(3,j)
               weight=weight3d5(j)
            elseif(lakon(nelem)(4:4).eq.'4') then
               xi=gauss3d4(1,j)
               et=gauss3d4(2,j)
               ze=gauss3d4(3,j)
               weight=weight3d4(j)
            elseif(lakon(nelem)(4:5).eq.'15') then
               xi=gauss3d8(1,j)
               et=gauss3d8(2,j)
               ze=gauss3d8(3,j)
               weight=weight3d8(j)
            elseif(lakon(nelem)(4:4).eq.'6') then
               xi=gauss3d7(1,j)
               et=gauss3d7(2,j)
               ze=gauss3d7(3,j)
               weight=weight3d7(j)
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
            do k=1,3
               coords(k,j)=0.d0
               do l=1,nope
                  coords(k,j)=coords(k,j)+xl(k,l)*shp(4,l)
               enddo
            enddo
         enddo
      endif
!
      if(prlab(ii)(1:4).eq.'S   ') then
         if(iorien.eq.0) then
            do j=1,mint3d
               write(5,'(i6,1x,i3,1p,6(1x,e11.4))') nelem,j,
     &              (stx(k,j,nelem),k=1,6)
            enddo
         else
            do j=1,mint3d
               call transformatrix(orab(1,iorien),coords(1,j),a)
               b(1,1)=stx(1,j,nelem)
               b(2,2)=stx(2,j,nelem)
               b(3,3)=stx(3,j,nelem)
               b(1,2)=stx(4,j,nelem)
               b(1,3)=stx(5,j,nelem)
               b(2,3)=stx(6,j,nelem)
               b(2,1)=b(1,2)
               b(3,1)=b(1,3)
               b(3,2)=b(2,3)
               do k=1,3
                  do l=1,3
                     c(k,l)=0.d0
                     do m=1,3
                        c(k,l)=c(k,l)+b(k,m)*a(m,l)
                     enddo
                  enddo
               enddo
               do k=1,3
                  do l=k,3
                     b(k,l)=0.d0
                     do m=1,3
                        b(k,l)=b(k,l)+a(m,k)*c(m,l)
                     enddo
                  enddo
               enddo
               write(5,'(i6,1x,i3,1p,6(1x,e11.4))') nelem,j,
     &              b(1,1),b(2,2),b(3,3),b(1,2),b(1,3),b(2,3)
            enddo
         endif
      elseif(prlab(ii)(1:4).eq.'E   ') then
         if(iorien.eq.0) then
            do j=1,mint3d
               write(5,'(i6,1x,i3,1p,6(1x,e11.4))') nelem,j,
     &              (eme(k,j,nelem),k=1,6)
            enddo
         else
            do j=1,mint3d
               call transformatrix(orab(1,iorien),coords(1,j),a)
               b(1,1)=eme(1,j,nelem)
               b(2,2)=eme(2,j,nelem)
               b(3,3)=eme(3,j,nelem)
               b(1,2)=eme(4,j,nelem)
               b(1,3)=eme(5,j,nelem)
               b(2,3)=eme(6,j,nelem)
               b(2,1)=b(1,2)
               b(3,1)=b(1,3)
               b(3,2)=b(2,3)
               do k=1,3
                  do l=1,3
                     do m=1,3
                        c(k,l)=b(k,m)*a(m,l)
                     enddo
                  enddo
               enddo
               do k=1,3
                  do l=k,3
                     do m=1,3
                        b(k,l)=a(m,k)*c(m,l)
                     enddo
                  enddo
               enddo
               write(5,'(i6,1x,i3,1p,6(1x,e11.4))') nelem,j,
     &              b(1,1),b(2,2),b(3,3),b(1,2),b(1,3),b(2,3)
            enddo
         endif
      elseif(prlab(ii)(1:4).eq.'PEEQ') then
         do j=1,mint3d
            write(5,'(i6,1x,i3,1p,6(1x,e11.4))') nelem,j,
     &           xstate(1,j,nelem)
         enddo
      elseif(prlab(ii)(1:4).eq.'ENER') then
         do j=1,mint3d
            write(5,'(i6,1x,i3,1p,6(1x,e11.4))') nelem,j,
     &           ener(j,nelem)
         enddo
      elseif(prlab(ii)(1:4).eq.'SDV ') then
         if(iorien.ne.0) then
            write(*,*) '*WARNING in printoutint: SDV cannot be'
            write(*,*) '         printed in the local system'
            write(*,*) '         results are in the global system'
         endif
         do j=1,mint3d
            write(5,'(i6,1x,i3,1p,99(1x,e11.4))') nelem,j,
     &           (xstate(k,j,nelem),k=1,nstate_)
         enddo
      elseif(prlab(ii)(1:4).eq.'HFL ') then
         if(iorien.eq.0) then
            do j=1,mint3d
               write(5,'(i6,1x,i3,1p,3(1x,e11.4))') nelem,j,
     &              (qfx(k,j,nelem),k=1,3)
            enddo
         else
            do j=1,mint3d
               do k=1,3
                  qfxl(k)=qfx(k,j,nelem)
               enddo
               call transformatrix(orab(1,iorien),coords(1,j),a)
               write(5,'(i6,1x,i3,1p,3(1x,e11.4))') nelem,j,
     &              qfxl(1)*a(1,1)+qfxl(2)*a(2,1)+qfxl(3)*a(3,1),
     &              qfxl(1)*a(1,2)+qfxl(2)*a(2,2)+qfxl(3)*a(3,2),
     &              qfxl(1)*a(1,3)+qfxl(2)*a(2,3)+qfxl(3)*a(3,3)
            enddo
         endif
      endif
!
      return
      end
