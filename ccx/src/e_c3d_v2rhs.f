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
      subroutine e_c3d_v2rhs(co,nk,konl,lakonl,
     &  ff,nelem,nmethod,vold,v,dtime,theta2,iexplicit,mi)
!
!     computation of the velocity element matrix and rhs for the element with
!     element with the topology in konl: step 3 (correction **)
!
!     ff: rhs 
!
      implicit none
!
      character*8 lakonl
!
      integer konl(20),nk,nelem,i,j,i1,j1,nmethod,jj,jj1,kk,
     &  nope,mint3d,iflag,iexplicit,mi(2)
!
      real*8 co(3,*),xl(3,20),shp(4,20),ff(60),xsjmod,vl(0:mi(2),20),
     &  vel(3),div,voldl(0:mi(2),20),v(0:mi(2),*),vold(0:mi(2),*),dtime,
     &  term,
     &  xi,et,ze,xsj,weight,shpv(20),theta2,dpress(3),ddpress(3)
!
      include "gauss.f"
!
      data iflag /3/
c      data iperm /13,14,-15,16,17,-18,19,20,-21,22,23,-24,
c     &            1,2,-3,4,5,-6,7,8,-9,10,11,-12,
c     &            37,38,-39,40,41,-42,43,44,-45,46,47,-48,
c     &            25,26,-27,28,29,-30,31,32,-33,34,35,-36,
c     &            49,50,-51,52,53,-54,55,56,-57,58,59,-60/
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
      elseif(lakonl(4:4).eq.'6') then
         mint3d=2
      else
         mint3d=0
      endif
!     
!     computation of the coordinates of the local nodes
!     
      do i=1,nope
         do j=1,3
            xl(j,i)=co(j,konl(i))
         enddo
      enddo
!     
!     initialisation for distributed forces
!     
      do i=1,3*nope
         ff(i)=0.d0
      enddo
!     
!     temperature, velocity and auxiliary variables
!     (rho*energy density, rho*velocity and rho)
!     
      do i1=1,nope
         do j1=1,4
            voldl(j1,i1)=vold(j1,konl(i1))
         enddo
c      write(*,*) 'pressure ',nelem,i1,voldl(4,i1)
      enddo
      if(iexplicit.eq.0) then
         do i1=1,nope
            vl(4,i1)=v(4,konl(i1))
         enddo
      endif
!     
!     computation of the matrix: loop over the Gauss points
!     
      do kk=1,mint3d
         if(lakonl(4:5).eq.'8R') then
            xi=gauss3d1(1,kk)
            et=gauss3d1(2,kk)
            ze=gauss3d1(3,kk)
            weight=weight3d1(kk)
         elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) 
     &           then
            xi=gauss3d2(1,kk)
            et=gauss3d2(2,kk)
            ze=gauss3d2(3,kk)
            weight=weight3d2(kk)
         elseif(lakonl(4:4).eq.'2') then
            xi=gauss3d3(1,kk)
            et=gauss3d3(2,kk)
            ze=gauss3d3(3,kk)
            weight=weight3d3(kk)
         elseif(lakonl(4:5).eq.'10') then
            xi=gauss3d5(1,kk)
            et=gauss3d5(2,kk)
            ze=gauss3d5(3,kk)
            weight=weight3d5(kk)
         elseif(lakonl(4:4).eq.'4') then
            xi=gauss3d4(1,kk)
            et=gauss3d4(2,kk)
            ze=gauss3d4(3,kk)
            weight=weight3d4(kk)
         elseif(lakonl(4:5).eq.'15') then
            xi=gauss3d8(1,kk)
            et=gauss3d8(2,kk)
            ze=gauss3d8(3,kk)
            weight=weight3d8(kk)
         elseif(lakonl(4:4).eq.'6') then
            xi=gauss3d7(1,kk)
            et=gauss3d7(2,kk)
            ze=gauss3d7(3,kk)
            weight=weight3d7(kk)
         endif
!     
!     calculation of the shape functions and their derivatives
!     in the gauss point
!     
         if(nope.eq.20) then
            if(lakonl(7:7).eq.'A') then
               call shape20h_ax(xi,et,ze,xl,xsj,shp,iflag)
            elseif((lakonl(7:7).eq.'E').or.(lakonl(7:7).eq.'S')) then
               call shape20h_pl(xi,et,ze,xl,xsj,shp,iflag)
            else
               call shape20h(xi,et,ze,xl,xsj,shp,iflag)
            endif
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
!     check the jacobian determinant
!     
         if(xsj.lt.1.d-20) then
            write(*,*) '*WARNING in e_c3d: nonpositive jacobian'
            write(*,*) '         determinant in element',nelem
            write(*,*)
            xsj=dabs(xsj)
            nmethod=0
         endif
!     
         xsjmod=dtime*xsj*weight
c         write(*,*) 'xsjmodin ec3dv2rhs ',xsjmod
!     
!     calculating the temperature temp, the velocity vel, the
!     divergence of the velocity div and the divergence
!     of the shape function times the velocity shpv(*)
!     in the integration point
!     
         do i1=1,3
            vel(i1)=0.d0
            dpress(i1)=0.d0
         enddo
c         div=0.d0
         do i1=1,nope
            do j1=1,3
               vel(j1)=vel(j1)+shp(4,i1)*voldl(j1,i1)
               dpress(j1)=dpress(j1)+shp(j1,i1)*voldl(4,i1)
c               div=div+shp(j1,i1)*voldl(j1,i1)
            enddo
         enddo
         do i1=1,nope
            shpv(i1)=shp(1,i1)*vel(1)+shp(2,i1)*vel(2)+
     &           shp(3,i1)*vel(3)
c            shpv(i1)=shp(1,i1)*vel(1)+shp(2,i1)*vel(2)+
c     &           shp(3,i1)*vel(3)+shp(4,i1)*div
         enddo
!
!        only for the semi-implicit procedure: calculate ddpress
!
         if(iexplicit.eq.0) then
            do i1=1,3
               ddpress(i1)=0.d0
            enddo
            div=0.d0
            do i1=1,nope
               do j1=1,3
                  ddpress(j1)=ddpress(j1)+shp(j1,i1)*vl(4,i1)
               enddo
            enddo
         endif
c         write(*,*) 'v2rhs.. ',nelem,kk,(dpress(i1),i1=1,3)
!     
!     determination of rhs
!     
         if(iexplicit.eq.1) then
            jj1=1
            do jj=1,nope
               term=xsjmod*(shp(4,jj)+dtime*shpv(jj)/2.d0)
               ff(jj1)=ff(jj1)-dpress(1)*term
               ff(jj1+1)=ff(jj1+1)-dpress(2)*term
               ff(jj1+2)=ff(jj1+2)-dpress(3)*term
               jj1=jj1+3
            enddo
         else
            jj1=1
            do jj=1,nope
!
!               without stability term 
!
c               ff(jj1)=ff(jj1)-xsjmod*(dpress(1)*(shp(4,jj)
c     &              )+theta2*shp(4,jj)*ddpress(1))
c               ff(jj1+1)=ff(jj1+1)-xsjmod*(dpress(2)*(shp(4,jj)
c     &              )+theta2*shp(4,jj)*ddpress(2))
c               ff(jj1+2)=ff(jj1+2)-xsjmod*(dpress(3)*(shp(4,jj)
c     &              )+theta2*shp(4,jj)*ddpress(3))
!
!               with stability term 
!
               ff(jj1)=ff(jj1)-xsjmod*(dpress(1)*(shp(4,jj)+
     &              (1.d0-theta2)*
     &              dtime*shpv(jj)/2.d0)+theta2*shp(4,jj)*ddpress(1))
               ff(jj1+1)=ff(jj1+1)-xsjmod*(dpress(2)*(shp(4,jj)+
     &              (1.d0-theta2)*
     &              dtime*shpv(jj)/2.d0)+theta2*shp(4,jj)*ddpress(2))
               ff(jj1+2)=ff(jj1+2)-xsjmod*(dpress(3)*(shp(4,jj)+
     &              (1.d0-theta2)*
     &              dtime*shpv(jj)/2.d0)+theta2*shp(4,jj)*ddpress(3))
               jj1=jj1+3
            enddo
         endif
!
      enddo
c      jj1=1
c      do jj=1,nope
c         write(*,*) 'force in ec3dv2rhs ',nelem,jj,
c     &          ff(jj1),ff(jj1+1),ff(jj1+2)
c         jj1=jj1+3
c      enddo
!     
!     for axially symmetric and plane stress/strain elements: 
!     complete s and sm
!
c      if((lakonl(6:7).eq.'RA').or.(lakonl(6:7).eq.'RS').or.
c     &   (lakonl(6:7).eq.'RE')) then
c!
c         if((nload.ne.0).or.(nbody.ne.0)) then
c            do i=1,60
c               k=abs(iperm(i))
c               ffax(i)=ff(k)*iperm(i)/k
c            enddo
c            do i=1,60
c               ff(i)=ff(i)+ffax(i)
c            enddo
c         endif
c!
c      endif
!
      return
      end

