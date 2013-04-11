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
      subroutine e_c3d_v1rhs(co,nk,konl,lakonl,p1,p2,omx,bodyfx,
     &  nbody,ff,nelem,nmethod,rhcon,nrhcon,ielmat,ntmat_,vold,voldaux,
     &  idist,dtime,matname,mi,
     &  ttime,time,istep,iinc,shcon,nshcon,
     &  turbulent,voldtu,yy,nelemface,sideface,nface,compressible)
!
!     computation of the velocity element matrix and rhs for the element with
!     element with the topology in konl: step 1 (correction *)
!
!     ff: rhs 
!
      implicit none
!
      integer turbulent,compressible
!
      character*1 sideface(*)
      character*8 lakonl
      character*80 matname(*),amat
!
      integer konl(20),ifaceq(8,6),nk,nbody,nelem,
     &  idist,i,j,k,i1,i2,j1,nmethod,ii,jj,jj1,id,ipointer,
     &  ig,kk,nrhcon(*),ielmat(*),nshcon(*),ntmat_,nope,nopes,imat,
     &  mint2d,mint3d,mi(2),ifacet(6,4),nopev,ifacew(8,5),istep,iinc,
     &  iflag,k1,nelemface(*),nface
!
      real*8 co(3,*),xl(3,20),shp(4,20),xs2(3,7),dvi,p1(3),p2(3),
     &  bodyf(3),bodyfx(3),ff(60),bf(3),q(3),c1,c2,xsjmod,
     &  rhcon(0:1,ntmat_,*),vel(3),div,shcon(0:3,ntmat_,*),
     &  voldl(0:mi(2),20),xl2(3,8),xsj2(3),shp2(7,8),omcor,
     &  vold(0:mi(2),*),om,omx,xi,et,ze,const,xsj,temp,theta2,
     &  voldaux(0:4,*),voldauxl(0:4,20),rho,weight,shpv(20),t(3,3),
     &  voldtu(2,*),voldtul(2,20),cvel(3),vkl(3,3),corio(3),xkin,
     &  xtuf,vort,un,yy(*),yyl(20),y,f2,unt,umt,a1,arg2,xlocal20(3,9,6),
     &  xlocal4(3,1,4),xlocal10(3,3,4),xlocal6(3,1,5),dpress(3),
     &  xlocal15(3,4,5),xlocal8(3,4,6),xlocal8r(3,1,6),xi3d,et3d,ze3d,
     &  shpvc(20)
!
      real*8 dtime,ttime,time,tvar(2)
!
      include "gauss.f"
      include "xlocal.f"
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
      data iflag /3/
      data a1 /0.31d0/
c      data iperm /13,14,-15,16,17,-18,19,20,-21,22,23,-24,
c     &            1,2,-3,4,5,-6,7,8,-9,10,11,-12,
c     &            37,38,-39,40,41,-42,43,44,-45,46,47,-48,
c     &            25,26,-27,28,29,-30,31,32,-33,34,35,-36,
c     &            49,50,-51,52,53,-54,55,56,-57,58,59,-60/
!
!     for pressure stability term (nithiarasu)
!
      theta2=.5d0
!
      tvar(1)=time
      tvar(2)=ttime+dtime
!
      imat=ielmat(nelem)
      amat=matname(imat)
!
      if(lakonl(4:4).eq.'2') then
         nope=20
         nopev=8
         nopes=8
      elseif(lakonl(4:4).eq.'8') then
         nope=8
         nopev=8
         nopes=4
      elseif(lakonl(4:5).eq.'10') then
         nope=10
         nopev=4
         nopes=6
      elseif(lakonl(4:4).eq.'4') then
         nope=4
         nopev=4
         nopes=3
      elseif(lakonl(4:5).eq.'15') then
         nope=15
         nopev=6
      elseif(lakonl(4:4).eq.'6') then
         nope=6
         nopev=6
      endif
!
      if(lakonl(4:5).eq.'8R') then
         mint2d=1
         mint3d=1
      elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) then
         if((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'S').or.
     &        (lakonl(7:7).eq.'E')) then
            mint2d=2
            mint3d=4
         else
            mint2d=4
            mint3d=8
         endif
      elseif(lakonl(4:4).eq.'2') then
         mint2d=9
         mint3d=27
      elseif(lakonl(4:5).eq.'10') then
         mint2d=3
         mint3d=4
      elseif(lakonl(4:4).eq.'4') then
         mint2d=1
         mint3d=1
      elseif(lakonl(4:5).eq.'15') then
         mint3d=9
      elseif(lakonl(4:4).eq.'6') then
         mint3d=2
      endif
!     
!     computation of the coordinates of the local nodes
!     
      do i=1,nope
         do j=1,3
            xl(j,i)=co(j,konl(i))
         enddo
c     write(*,*) i,konl(i),(xl(j,i),j=1,3)
      enddo
!     
!     initialisation for distributed forces
!     
      do i=1,3*nope
         ff(i)=0.d0
      enddo
!     
!     temperature, velocity, auxiliary variables
!     (rho*velocity and rho) and if turbulent 
!     rho*turbulence variables
!     
      do i1=1,nope
c         do i2=0,3
         do i2=0,4
            voldl(i2,i1)=vold(i2,konl(i1))
         enddo
         do i2=1,4
            voldauxl(i2,i1)=voldaux(i2,konl(i1))
         enddo
         if(turbulent.ne.0) then
            voldtul(1,i1)=voldtu(1,konl(i1))
            voldtul(2,i1)=voldtu(2,konl(i1))
            yyl(i1)=yy(konl(i1))
         endif
c         write(*,*) nelem,i1,konl(i1),(voldl(j,i1),j=0,3)
c         write(*,*) nelem,i1,konl(i1),(voldauxl(j,i1),j=1,4)
      enddo
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
!     
!        calculating of
!        the temperature temp
!        the velocity vel
!        the auxiliary variable cvel
!        the velocity gradient vkl
!        the divergence of the velocity div 
!        the divergence of the shape function times the velocity shpv(*)
!              in the integration point
!     
         temp=0.d0
         do i1=1,3
            vel(i1)=0.d0
            cvel(i1)=0.d0
            dpress(i1)=0.d0
            do j1=1,3
               vkl(i1,j1)=0.d0
            enddo
         enddo
         div=0.d0
         do i1=1,nope
            temp=temp+shp(4,i1)*voldl(0,i1)
            do j1=1,3
               vel(j1)=vel(j1)+shp(4,i1)*voldl(j1,i1)
               dpress(j1)=dpress(j1)+shp(j1,i1)*voldl(4,i1)
               do k1=1,3
                  vkl(j1,k1)=vkl(j1,k1)+shp(k1,i1)*voldl(j1,i1)
               enddo
            enddo
         enddo
         if(compressible.eq.1) then
            div=vkl(1,1)+vkl(2,2)+vkl(3,3)
         else
            div=0.d0
         endif
!
         do i1=1,nope
            shpv(i1)=shp(1,i1)*vel(1)+shp(2,i1)*vel(2)+
     &           shp(3,i1)*vel(3)+shp(4,i1)*div
cv            shpv(i1)=shp(1,i1)*voldl(1,i1)+shp(2,i1)*voldl(2,i1)+
cv     &           shp(3,i1)*voldl(3,i1)
            do j1=1,3
               cvel(j1)=cvel(j1)+shpv(i1)*voldauxl(j1,i1)
            enddo
         enddo
!     
!     material data (density and dynamic viscosity)
!     
c         call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,r,dvi,
c     &        rhcon,nrhcon,rho)
         call materialdata_dvi(imat,ntmat_,temp,shcon,nshcon,dvi)
!
!     determining the dissipative stress 
!
         do i1=1,3
            do j1=i1,3
               t(i1,j1)=vkl(i1,j1)+vkl(j1,i1)
            enddo
cccc
            if(compressible.eq.1) t(i1,i1)=t(i1,i1)-2.d0*div/3.d0
cccc
         enddo
!     
!     calculation of the density for gases
!     
!     calculation of the turbulent kinetic energy, turbulence
!     frequency and their spatial derivatives for gases and liquids
!
         if(compressible.eq.1) then
!
!           gas
!
            rho=0.d0
            do i1=1,nope
               rho=rho+shp(4,i1)*voldauxl(4,i1)
            enddo
            if(turbulent.ne.0) then
               xkin=0.d0
               xtuf=0.d0
               y=0.d0
               do i1=1,nope
                  xkin=xkin+shp(4,i1)*voldtul(1,i1)
                  xtuf=xtuf+shp(4,i1)*voldtul(2,i1)
                  y=y+shp(4,i1)*yyl(i1)
               enddo
               xkin=xkin/rho
               xtuf=xtuf/rho
            endif
         else
!
!           liquid
!
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_)
            if(turbulent.ne.0) then
               xkin=0.d0
               xtuf=0.d0
               y=0.d0
               do i1=1,nope
                  xkin=xkin+shp(4,i1)*voldtul(1,i1)
                  xtuf=xtuf+shp(4,i1)*voldtul(2,i1)
                  y=y+shp(4,i1)*yyl(i1)
               enddo
               xkin=xkin/rho
               xtuf=xtuf/rho
            endif
         endif
!
!        calculation of turbulent auxiliary variables
!
         if(turbulent.ne.0) then
!     
!           vorticity
!     
            vort=dsqrt((vkl(3,2)-vkl(2,3))**2+
     &           (vkl(1,3)-vkl(3,1))**2+
     &           (vkl(2,1)-vkl(1,2))**2)
!     
!           kinematic viscosity
!     
            un=dvi/rho
!     
!           factor F2
!     
            c1=dsqrt(xkin)/(0.09d0*xtuf*y)
            c2=500.d0*un/(y*y*xtuf)
            arg2=max(2.d0*c1,c2)
            f2=dtanh(arg2*arg2)
!     
!     kinematic and dynamic turbulent viscosity
!     
            unt=a1*xkin/max(a1*xtuf,vort*f2)
            umt=unt*rho
            write(*,*) 'e_c3d_v1rhs ',dvi,umt
!
            do i1=1,3
               do j1=i1,3
                  t(i1,j1)=(dvi+umt)*t(i1,j1)
               enddo
               t(i1,i1)=t(i1,i1)-2.d0*rho*xkin/3.d0
            enddo
         else
!
            do i1=1,3
               do j1=i1,3
                  t(i1,j1)=dvi*t(i1,j1)
               enddo
            enddo
         endif
         t(2,1)=t(1,2)
         t(3,1)=t(1,3)
         t(3,2)=t(2,3)
!     
!     initialisation for the body forces
!     
         om=omx*rho
         omcor=2.d0*rho*dsqrt(omx)
         if(nbody.ne.0) then
            do ii=1,3
               bodyf(ii)=bodyfx(ii)*rho
            enddo
         endif
!     
!     determination of lhs and rhs
!     
         jj1=1
         do jj=1,nope
c
c           convective + diffusive
c
            ff(jj1)=ff(jj1)-xsjmod*
     &            (cvel(1)*(shp(4,jj)+dtime*shpv(jj)/2.d0)
     &             +shp(1,jj)*t(1,1)+shp(2,jj)*t(1,2)+shp(3,jj)*t(1,3))
            ff(jj1+1)=ff(jj1+1)-xsjmod*
     &            (cvel(2)*(shp(4,jj)+dtime*shpv(jj)/2.d0)
     &             +shp(1,jj)*t(2,1)+shp(2,jj)*t(2,2)+shp(3,jj)*t(2,3))
            ff(jj1+2)=ff(jj1+2)-xsjmod*
     &            (cvel(3)*(shp(4,jj)+dtime*shpv(jj)/2.d0)
     &             +shp(1,jj)*t(3,1)+shp(2,jj)*t(3,2)+shp(3,jj)*t(3,3))
c
c           convective1
c
c            ff(jj1)=ff(jj1)-xsjmod*
c     &            (cvel(1)*(shp(4,jj)))
c            ff(jj1+1)=ff(jj1+1)-xsjmod*
c     &            (cvel(2)*(shp(4,jj)))
c            ff(jj1+2)=ff(jj1+2)-xsjmod*
c     &            (cvel(3)*(shp(4,jj)))
c
c           convective2
c
c            ff(jj1)=ff(jj1)-xsjmod*
c     &            (cvel(1)*(dtime*shpv(jj)/2.d0))
c            ff(jj1+1)=ff(jj1+1)-xsjmod*
c     &            (cvel(2)*(dtime*shpv(jj)/2.d0))
c            ff(jj1+2)=ff(jj1+2)-xsjmod*
c     &            (cvel(3)*(dtime*shpv(jj)/2.d0))
c
c           diffusive
c
c            ff(jj1)=ff(jj1)-xsjmod*
c     &             (shp(1,jj)*t(1,1)+shp(2,jj)*t(1,2)+shp(3,jj)*t(1,3))
c            ff(jj1+1)=ff(jj1+1)-xsjmod*
c     &             (shp(1,jj)*t(2,1)+shp(2,jj)*t(2,2)+shp(3,jj)*t(2,3))
c            ff(jj1+2)=ff(jj1+2)-xsjmod*
c     &             (shp(1,jj)*t(3,1)+shp(2,jj)*t(3,2)+shp(3,jj)*t(3,3))
c
c           pressure stability term (nithiarasu)
c
c            ff(jj1)=ff(jj1)-xsjmod*dpress(1)*(1.d0-theta2)*
c     &              dtime*shpv(jj)/2.d0
c            ff(jj1+1)=ff(jj1+1)-xsjmod*dpress(2)*(1.d0-theta2)*
c     &              dtime*shpv(jj)/2.d0
c            ff(jj1+2)=ff(jj1+2)-xsjmod*dpress(3)*(1.d0-theta2)*
c     &              dtime*shpv(jj)/2.d0
            jj1=jj1+3
         enddo
!     
!     computation of contribution due to body forces
!     
         if(nbody.ne.0) then
            if(om.gt.0.d0) then
               do i1=1,3
!     
!     computation of the global coordinates of the gauss
!     point
!     
                  q(i1)=0.d0
                  do j1=1,nope
                     q(i1)=q(i1)+shp(4,j1)*xl(i1,j1)
                  enddo
!     
                  q(i1)=q(i1)-p1(i1)
               enddo
               const=q(1)*p2(1)+q(2)*p2(2)+q(3)*p2(3)
!
!              Coriolis forces
!
               corio(1)=vel(2)*p2(3)-vel(3)*p2(2)
               corio(2)=vel(3)*p2(1)-vel(1)*p2(3)
               corio(3)=vel(1)*p2(2)-vel(2)*p2(1)
!     
!     inclusion of the centrifugal force into the body force
!     
               do i1=1,3
c
                  bf(i1)=bodyf(i1)+(q(i1)-const*p2(i1))*om+
     &                   corio(i1)*omcor
c                  bf(i1)=bodyf(i1)+(q(i1)-const*p2(i1))*om
               enddo
            else
               do i1=1,3
                  bf(i1)=bodyf(i1)
               enddo
            endif
            jj1=1
            do jj=1,nope
               ff(jj1)=ff(jj1)+xsjmod*bf(1)*(shp(4,jj)+
     &              dtime*shpv(jj)/2.d0)
               ff(jj1+1)=ff(jj1+1)+xsjmod*bf(2)*(shp(4,jj)+
     &              dtime*shpv(jj)/2.d0)
               ff(jj1+2)=ff(jj1+2)+xsjmod*bf(3)*(shp(4,jj)+
     &              dtime*shpv(jj)/2.d0)
               jj1=jj1+3
            enddo
         endif
!     
      enddo
!     
      if(nface.ne.0) then
!     
!        free stream or solid surface boundaries
!     
         call nident(nelemface,nelem,nface,id)
         do
            if((id.eq.0).or.(nelemface(id).ne.nelem)) exit
c            read(sideface(id)(1:1),'(i1)') ig
            ig=ichar(sideface(id)(1:1))-48
!     
!     treatment of wedge faces
!     
            if(lakonl(4:4).eq.'6') then
               mint2d=1
               if(ig.le.2) then
                  nopes=3
               else
                  nopes=4
               endif
            endif
            if(lakonl(4:5).eq.'15') then
               if(ig.le.2) then
                  mint2d=3
                  nopes=6
               else
                  mint2d=4
                  nopes=8
               endif
            endif
!     
            if((nope.eq.20).or.(nope.eq.8)) then
               do i=1,nopes
                  do j=1,3
                     xl2(j,i)=co(j,konl(ifaceq(i,ig)))
                  enddo
               enddo
            elseif((nope.eq.10).or.(nope.eq.4)) then
               do i=1,nopes
                  do j=1,3
                     xl2(j,i)=co(j,konl(ifacet(i,ig)))
                  enddo
               enddo
            else
               do i=1,nopes
                  do j=1,3
                     xl2(j,i)=co(j,konl(ifacew(i,ig)))
                  enddo
               enddo
            endif
!     
            do i=1,mint2d
!     
!              local coordinates of the surface integration
!              point within the surface local coordinate system
!     
               if((lakonl(4:5).eq.'8R').or.
     &              ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
                  xi=gauss2d1(1,i)
                  et=gauss2d1(2,i)
                  weight=weight2d1(i)
               elseif((lakonl(4:4).eq.'8').or.
     &                 (lakonl(4:6).eq.'20R').or.
     &                 ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
                  xi=gauss2d2(1,i)
                  et=gauss2d2(2,i)
                  weight=weight2d2(i)
               elseif(lakonl(4:4).eq.'2') then
                  xi=gauss2d3(1,i)
                  et=gauss2d3(2,i)
                  weight=weight2d3(i)
               elseif((lakonl(4:5).eq.'10').or.
     &                 ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
                  xi=gauss2d5(1,i)
                  et=gauss2d5(2,i)
                  weight=weight2d5(i)
               elseif((lakonl(4:4).eq.'4').or.
     &                 ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
                  xi=gauss2d4(1,i)
                  et=gauss2d4(2,i)
                  weight=weight2d4(i)
               endif
!     
!              local surface normal
!
               if(nopes.eq.8) then
                  call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
               elseif(nopes.eq.4) then
                  call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
               elseif(nopes.eq.6) then
                  call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
               else
                  call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
               endif
!
!              local coordinates of the surface integration
!              point within the element local coordinate system
!     
               if(lakonl(4:5).eq.'8R') then
                  xi3d=xlocal8r(1,i,ig)
                  et3d=xlocal8r(2,i,ig)
                  ze3d=xlocal8r(3,i,ig)
                  call shape8h(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
               elseif(lakonl(4:4).eq.'8') then
                  xi3d=xlocal8(1,i,ig)
                  et3d=xlocal8(2,i,ig)
                  ze3d=xlocal8(3,i,ig)
                  call shape8h(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
               elseif(lakonl(4:6).eq.'20R') then
                  xi3d=xlocal8(1,i,ig)
                  et3d=xlocal8(2,i,ig)
                  ze3d=xlocal8(3,i,ig)
                  call shape20h(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
               elseif(lakonl(4:4).eq.'2') then
                  xi3d=xlocal20(1,i,ig)
                  et3d=xlocal20(2,i,ig)
                  ze3d=xlocal20(3,i,ig)
                  call shape20h(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
               elseif(lakonl(4:5).eq.'10') then
                  xi3d=xlocal10(1,i,ig)
                  et3d=xlocal10(2,i,ig)
                  ze3d=xlocal10(3,i,ig)
                  call shape10tet(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
               elseif(lakonl(4:4).eq.'4') then
                  xi3d=xlocal4(1,i,ig)
                  et3d=xlocal4(2,i,ig)
                  ze3d=xlocal4(3,i,ig)
                  call shape4tet(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
               elseif(lakonl(4:5).eq.'15') then
                  xi3d=xlocal15(1,i,ig)
                  et3d=xlocal15(2,i,ig)
                  ze3d=xlocal15(3,i,ig)
                  call shape15w(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
               elseif(lakonl(4:4).eq.'6') then
                  xi3d=xlocal6(1,i,ig)
                  et3d=xlocal6(2,i,ig)
                  ze3d=xlocal6(3,i,ig)
                  call shape6w(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
               endif
!
ce               xsjmod=dtime*weight
!     
!     calculating of
!     the temperature temp
!     the velocity gradient vkl
!     in the integration point
!     
               temp=0.d0
               do i1=1,3
ce                  vel(i1)=0.d0
ce                  cvel(i1)=0.d0
                  do j1=1,3
                     vkl(i1,j1)=0.d0
                  enddo
               enddo
               do i1=1,nope
                  temp=temp+shp(4,i1)*voldl(0,i1)
                  do j1=1,3
ce                     vel(j1)=vel(j1)+shp(4,i1)*voldl(j1,i1)
                     do k1=1,3
                        vkl(j1,k1)=vkl(j1,k1)+shp(k1,i1)*voldl(j1,i1)
                     enddo
                  enddo
               enddo
               if(compressible.eq.1) then
                  div=vkl(1,1)+vkl(2,2)+vkl(3,3)
               else
                  div=0.d0
               endif
c               div=vkl(1,1)+vkl(2,2)+vkl(3,3)
ce               do i1=1,nope
ce                  shpv(i1)=shp(1,i1)*vel(1)+shp(2,i1)*vel(2)+
ce     &                 shp(3,i1)*vel(3)
cec     shpv(i1)=shp(1,i1)*vel(1)+shp(2,i1)*vel(2)+
cec     &           shp(3,i1)*vel(3)+shp(4,i1)*div
ce                  do j1=1,3
ce                    cvel(j1)=cvel(j1)+shpv(i1)*voldauxl(j1,i1)
ce                  enddo
ce               enddo
!     
!     material data (density and dynamic viscosity)
!     
c               call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,r,
c     &              dvi,rhcon,nrhcon,rho)
               call materialdata_dvi(imat,ntmat_,temp,shcon,nshcon,dvi)
!     
!     determining the dissipative stress 
!     
               do i1=1,3
                  do j1=i1,3
                     t(i1,j1)=vkl(i1,j1)+vkl(j1,i1)
                  enddo
                  if(compressible.eq.1) t(i1,i1)=t(i1,i1)-2.d0*div/3.d0
c                  t(i1,i1)=t(i1,i1)-2.d0*div/3.d0
c                  t(i1,i1)=t(i1,i1)
               enddo
!     
!              calculation of the density for gases
!     
!     calculation of the turbulent kinetic energy, turbulence
!     frequency and their spatial derivatives for gases and liquids
!     
               if(compressible.eq.1) then
!     
!                 gas
!     
                  rho=0.d0
                  do i1=1,nope
                     rho=rho+shp(4,i1)*voldauxl(4,i1)
                  enddo
                  if(turbulent.ne.0) then
                     xkin=0.d0
                     xtuf=0.d0
                     y=0.d0
                     do i1=1,nope
                        xkin=xkin+shp(4,i1)*voldtul(1,i1)/voldauxl(4,i1)
                        xtuf=xtuf+shp(4,i1)*voldtul(2,i1)/voldauxl(4,i1)
                        y=y+shp(4,i1)*yyl(i1)
                     enddo
                  endif
               else
!     
!                 liquid
!     
                  call materialdata_rho(rhcon,nrhcon,imat,rho,
     &                 temp,ntmat_)
                  if(turbulent.ne.0) then
                     xkin=0.d0
                     xtuf=0.d0
                     y=0.d0
                     do i1=1,nope
                        xkin=xkin+shp(4,i1)*voldtul(1,i1)/rho
                        xtuf=xtuf+shp(4,i1)*voldtul(2,i1)/rho
                        y=y+shp(4,i1)*yyl(i1)
                     enddo
                  endif
               endif
!     
!              calculation of turbulent auxiliary variables
!     
               if(turbulent.ne.0) then
!     
!                 vorticity
!     
                  vort=dsqrt((vkl(3,2)-vkl(2,3))**2+
     &                 (vkl(1,3)-vkl(3,1))**2+
     &                 (vkl(2,1)-vkl(1,2))**2)
!     
!                 kinematic viscosity
!     
                  un=dvi/rho
!     
!                 factor F2
!     
                  c1=dsqrt(xkin)/(0.09d0*xtuf*y)
                  c2=500.d0*un/(y*y*xtuf)
                  arg2=max(2.d0*c1,c2)
                  f2=dtanh(arg2*arg2)
!     
!                 kinematic and dynamic turbulent viscosity
!     
                  unt=a1*xkin/max(a1*xtuf,vort*f2)
                  umt=unt*rho
!     
                  do i1=1,3
                     do j1=i1,3
                        t(i1,j1)=(dvi+umt)*t(i1,j1)
                     enddo
                     t(i1,i1)=t(i1,i1)-2.d0*rho*xkin/3.d0
                  enddo
               else
!     
                  do i1=1,3
                     do j1=i1,3
                        t(i1,j1)=dvi*t(i1,j1)
                     enddo
                  enddo
               endif
               t(2,1)=t(1,2)
               t(3,1)=t(1,3)
               t(3,2)=t(2,3)
!     
!     computation of contribution due to body forces
!     
ce               if(nbody.ne.0) then
ce                  if(om.gt.0.d0) then
ce                     do i1=1,3
ce!     
ce!     computation of the global coordinates of the gauss
ce!     point
ce!     
ce                        q(i1)=0.d0
ce                        do j1=1,nope
ce                           q(i1)=q(i1)+shp(4,j1)*xl(i1,j1)
ce                        enddo
ce!     
ce                        q(i1)=q(i1)-p1(i1)
ce                     enddo
ce                    const=q(1)*p2(1)+q(2)*p2(2)+q(3)*p2(3)
ce!     
ce!     Coriolis forces
ce!     
ce                     corio(1)=vel(2)*p2(3)-vel(3)*p2(2)
ce                     corio(2)=vel(3)*p2(1)-vel(1)*p2(3)
ce                     corio(3)=vel(1)*p2(2)-vel(2)*p2(1)
ce!     
ce!     inclusion of the centrifugal force into the body force
ce!     
ce                     do i1=1,3
ce                        bf(i1)=bodyf(i1)+(q(i1)-const*p2(i1)+corio(i1))*om
ce                     enddo
ce                  else
ce                     do i1=1,3
ce                        bf(i1)=bodyf(i1)
ce                     enddo
ce                  endif
ce               endif
ce
ce               do i1=1,3
ce                  tm(i1)=0.d0
ce                  do j1=1,3
ce                     tm(i1)=tm(i1)+xsjmod*(t(i1,j1)+dtime*vel(j1)*
ce     &                      (cvel(i1)-bf(i1))/2.d0)*xsj2(j1)
ce                  enddo
ce               enddo
!               
               do k=1,nopes
                  if((nope.eq.20).or.(nope.eq.8)) then
                     ipointer=(ifaceq(k,ig)-1)*3
                  elseif((nope.eq.10).or.(nope.eq.4)) then
                     ipointer=(ifacet(k,ig)-1)*3
                  else
                     ipointer=(ifacew(k,ig)-1)*3
                  endif
                  ff(ipointer+1)=ff(ipointer+1)+shp2(4,k)*
     &              (t(1,1)*xsj2(1)+t(1,2)*xsj2(2)+t(1,3)*xsj2(3))*
     &              weight*dtime
                  ff(ipointer+2)=ff(ipointer+2)+shp2(4,k)*
     &              (t(2,1)*xsj2(1)+t(2,2)*xsj2(2)+t(2,3)*xsj2(3))*
     &              weight*dtime
                  ff(ipointer+3)=ff(ipointer+3)+shp2(4,k)*
     &              (t(3,1)*xsj2(1)+t(3,2)*xsj2(2)+t(3,3)*xsj2(3))*
     &              weight*dtime
ce
ce                  the 9 lines above have to be replaced by the
ce                  next three
ce
ce                  ff(ipointer+1)=ff(ipointer+1)+shp2(4,k)*tm(1)
ce                  ff(ipointer+2)=ff(ipointer+2)+shp2(4,k)*tm(2)
ce                  ff(ipointer+3)=ff(ipointer+3)+shp2(4,k)*tm(3)
               enddo
            enddo
            id=id-1
         enddo
      endif
!     
!     
!     for axially symmetric and plane stress/strain elements: 
!     complete s 
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
c!
      return
      end

