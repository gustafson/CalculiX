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
      subroutine e_c3d_krhs(co,nk,konl,lakonl,ffk,fft,nelem,nmethod,
     &  rhcon,nrhcon,ielmat,ntmat_,vold,voldaux,dtime,matname,mint_,
     &  shcon,nshcon,voldtu,compressible,yy,nelemface,sideface,nface,
     &  turbulent)
!
!     computation of the turbulence element matrix and rhs for the
!     element with the topology in konl: step 4
!
!     ffk and fft: rhs (x 2: kinetic term and turbulence frequency term):
!
      implicit none
!
      character*1 sideface(*)
      character*8 lakonl
      character*80 matname(*),amat
!
      integer konl(20),ifaceq(8,6),nk,nelem,nload,i,j,k,i1,i2,j1,k1,
     &  nmethod,ii,jj,id,ipointer,ig,kk,nrhcon(*),ielmat(*),nshcon(*),
     &  ntmat_,nope,nopes,imat,mint2d,mint3d,mint_,ifacet(6,4),nopev,
     &  ifacew(8,5),istep,iinc,layer,kspt,jltyp,iflag,iscale,
     &  compressible,idf,igl,nelemface(*),nface,turbulent
!
      real*8 co(3,*),xl(3,20),shp(4,20),xs2(3,7),dvi,
     &  ffk(60),xsjmod,vkl(3,3),rhcon(0:1,ntmat_,*),reltime,
     &  t(3,3),bfv,press,vel(3),div,shcon(0:3,ntmat_,*),pgauss(3),
     &  xkin,xtuf,voldl(0:4,20),yyl(20),tvk(3),tvt(3),
     &  xl2(0:3,8),xsj2(3),shp2(7,8),vold(0:4,*),tvnk,tvnt,
     &  om,omx,xi,et,ze,const,xsj,fft(60),dxkin(3),
     &  temp,voldaux(0:4,*),voldauxl(0:4,20),rho,dxtuf(3),
     &  weight,shpv(20),rhokin,rhotuf,y,vort,c1,c2,arg2,f2,
     &  a1,unt,umt,cdktuf,arg1,f1,skin,skin1,skin2,stuf,stuf1,
     &  stuf2,beta,beta1,beta2,betas,gamm,gamm1,xkappa,un,
     &  gamm2,umsk,umst,tu,tuk,tut,voldtu(2,*),voldtul(2,20),
     &  f1m,yy(*),xsjmodk,xsjmodt,xi3d,et3d,ze3d,xlocal20(3,9,6),
     &  xlocal4(3,1,4),xlocal10(3,3,4),xlocal6(3,1,5),
     &  xlocal15(3,4,5),xlocal8(3,4,6),xlocal8r(3,1,6)
!
      real*8 dtime,ttime,time,tvar(2),
     &  coords(3)
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
c      data iperm /13,14,-15,16,17,-18,19,20,-21,22,23,-24,
c     &            1,2,-3,4,5,-6,7,8,-9,10,11,-12,
c     &            37,38,-39,40,41,-42,43,44,-45,46,47,-48,
c     &            25,26,-27,28,29,-30,31,32,-33,34,35,-36,
c     &            49,50,-51,52,53,-54,55,56,-57,58,59,-60/
!
!     turbulence constants
!
      data skin1 /0.85d0/
      data skin2 /1.d0/
      data stuf1 /0.5d0/
      data stuf2 /0.856d0/
      data beta1 /0.075d0/
      data beta2 /0.0828d0/
      data a1 /0.31d0/
      data betas /0.09d0/
      data xkappa / 0.41d0/
!
      gamm1=beta1/betas-stuf1*xkappa*xkappa/dsqrt(betas)
      gamm2=beta2/betas-stuf2*xkappa*xkappa/dsqrt(betas)
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
!       initialisation for distributed forces
!
      do i=1,nope
         ffk(i)=0.d0
         fft(i)=0.d0
      enddo
!
!     temperature, velocity and auxiliary variables
!     (rho*energy density, rho*velocity and rho)
!
         do i1=1,nope
            do i2=0,3
               voldl(i2,i1)=vold(i2,konl(i1))
            enddo
            voldtul(1,i1)=voldtu(1,konl(i1))
            voldtul(2,i1)=voldtu(2,konl(i1))
            yyl(i1)=yy(konl(i1))
         enddo
         if(compressible.eq.1) then
            do i1=1,nope
               voldauxl(4,i1)=voldaux(4,konl(i1))
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
!     
!        calculating of
!        rho times turbulent kinetic energy times shpv(*): rhokin
!        rho times turbulence frequency times shpv(*): rhotuf
!        distance from solid surface: y
!        the velocity vel
!        the velocity gradient vkl
!        the divergence of the shape function times the velocity shpv(*)
!              in the integration point
!     
         temp=0.d0
         rhokin=0.d0
         rhotuf=0.d0
         y=0.d0
         do i1=1,3
            vel(i1)=0.d0
            do j1=1,3
               vkl(i1,j1)=0.d0
            enddo
         enddo
         do i1=1,nope
            temp=temp+shp(4,i1)*voldl(0,i1)
            y=y+shp(4,i1)*yyl(i1)
            do j1=1,3
               vel(j1)=vel(j1)+shp(4,i1)*voldl(j1,i1)
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
         do i1=1,nope
            shpv(i1)=shp(1,i1)*vel(1)+shp(2,i1)*vel(2)+
     &           shp(3,i1)*vel(3)+shp(4,i1)*div
            rhokin=rhokin+shpv(i1)*voldtul(1,i1)
            rhotuf=rhotuf+shpv(i1)*voldtul(2,i1)
         enddo
!     
!     material data (density and dynamic viscosity)
!     
c         call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,r,dvi,
c     &        rhcon,nrhcon,rho)
         call materialdata_dvi(imat,ntmat_,temp,shcon,nshcon,dvi)
!
!     determining the stress
!
         do i1=1,3
            do j1=1,3
               t(i1,j1)=vkl(i1,j1)+vkl(j1,i1)
            enddo
            if(compressible.eq.1) t(i1,i1)=t(i1,i1)-2.d0*div/3.d0
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
            xkin=0.d0
            xtuf=0.d0
            do j1=1,3
               dxkin(j1)=0.d0
               dxtuf(j1)=0.d0
            enddo
            do i1=1,nope
               rho=rho+shp(4,i1)*voldauxl(4,i1)
               xkin=xkin+shp(4,i1)*voldtul(1,i1)
               xtuf=xtuf+shp(4,i1)*voldtul(2,i1)
               do j1=1,3
                  dxkin(j1)=dxkin(j1)+shp(j1,i1)*voldtul(1,i1)
                  dxtuf(j1)=dxtuf(j1)+shp(j1,i1)*voldtul(2,i1)
               enddo
            enddo
            xkin=xkin/rho
            xtuf=xtuf/rho
            do j1=1,3
               dxkin(j1)=dxkin(j1)/rho
               dxtuf(j1)=dxtuf(j1)/rho
            enddo
         else
!
!           liquid
!
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_)
            xkin=0.d0
            xtuf=0.d0
            do j1=1,3
               dxkin(j1)=0.d0
               dxtuf(j1)=0.d0
            enddo
            do i1=1,nope
               xkin=xkin+shp(4,i1)*voldtul(1,i1)
               xtuf=xtuf+shp(4,i1)*voldtul(2,i1)
               do j1=1,3
                  dxkin(j1)=dxkin(j1)+shp(j1,i1)*voldtul(1,i1)
                  dxtuf(j1)=dxtuf(j1)+shp(j1,i1)*voldtul(2,i1)
               enddo
            enddo
            xkin=xkin/rho
            xtuf=xtuf/rho
            do j1=1,3
               dxkin(j1)=dxkin(j1)/rho
               dxtuf(j1)=dxtuf(j1)/rho
            enddo
         endif
!
!        calculation of turbulent auxiliary variables
!
!        vorticity
!
         vort=dsqrt((vkl(3,2)-vkl(2,3))**2+
     &              (vkl(1,3)-vkl(3,1))**2+
     &              (vkl(2,1)-vkl(1,2))**2)
!
!        kinematic viscosity
!
         un=dvi/rho
!
!        factor F2
!
         c1=dsqrt(xkin)/(0.09d0*xtuf*y)
         c2=500.d0*un/(y*y*xtuf)
         arg2=max(2.d0*c1,c2)
         f2=dtanh(arg2*arg2)
!
!        kinematic and dynamic turbulent viscosity
!
         unt=a1*xkin/max(a1*xtuf,vort*f2)
         umt=unt*rho
!
!        factor F1
!     
         if(turbulent.eq.1) then
!
!           k-epsilon model
!
            f1=1.d0
         elseif(turbulent.eq.2) then
!
!           q-omega model
!
            f1=0.d0
         else
!
!           SST model
!
            cdktuf=max(2.d0*rho*stuf2*
     &         (dxkin(1)*dxtuf(1)+dxkin(2)*dxtuf(2)+dxkin(3)*dxtuf(3))/
     &         xtuf,1.d-20)
            arg1=min(max(c1,c2),4.d0*rho*stuf2*xkin/(cdktuf*y*y))
            f1=dtanh(arg1**4.d0)
         endif
         f1m=1.d0-f1
!
!        interpolation of the constants
!
         skin=f1*skin1+f1m*skin2
         stuf=f1*stuf1+f1m*stuf2
         beta=f1*beta1+f1m*beta2
         gamm=f1*gamm1+f1m*gamm2
!
!        auxiliary quantities
!
         umsk=dvi+skin*umt
         umst=dvi+stuf*umt
         tu=umt*(t(1,1)*vkl(1,1)+t(1,2)*vkl(1,2)+t(1,3)*vkl(1,3)+
     &           t(2,1)*vkl(2,1)+t(2,2)*vkl(2,2)+t(2,3)*vkl(2,3)+
     &           t(3,1)*vkl(3,1)+t(3,2)*vkl(3,2)+t(3,3)*vkl(3,3))
!
         if(compressible.eq.1) then
            tu=tu-2.d0*rho*xkin*div/3.d0
         endif
!
         tuk=tu-betas*rho*xtuf*xkin
         tut=gamm*tu/unt-beta*rho*xtuf*xtuf+2.d0*f1m*rho*stuf2*
     &       (dxkin(1)*dxtuf(1)+dxkin(2)*dxtuf(2)+dxkin(3)*dxtuf(3))/
     &       xtuf
         do i1=1,3
            dxkin(i1)=dxkin(i1)*umsk
            dxtuf(i1)=dxtuf(i1)*umst
         enddo
!     
!     determination of lhs and rhs
!     
         do jj=1,nope
!     
            ffk(jj)=ffk(jj)-xsjmod*((shp(4,jj)+dtime*shpv(jj)/2.d0)*
     &           (rhokin-tuk)+(shp(1,jj)*dxkin(1)+shp(2,jj)*dxkin(2)
     &              +shp(3,jj)*dxkin(3)))
            fft(jj)=fft(jj)-xsjmod*((shp(4,jj)+dtime*shpv(jj)/2.d0)*
     &           (rhotuf-tut)+(shp(1,jj)*dxtuf(1)+shp(2,jj)*dxtuf(2)
     &              +shp(3,jj)*dxtuf(3)))
         enddo
         
!     
      enddo
!     
      if(nface.ne.0) then
!     
!        free stream or solid surface boundaries
!     
         call nident(nelemface,nelem,nface,idf)
         do
            if((idf.eq.0).or.(nelemface(idf).ne.nelem)) exit
            read(sideface(idf)(1:1),'(i1)') ig
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
!     material data (density and dynamic viscosity)
!     
c               call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,r,
c     &              dvi,rhcon,nrhcon,rho)
               call materialdata_dvi(imat,ntmat_,temp,shcon,nshcon,dvi)
!     
!     calculation of the density for gases
!     
!     calculation of the turbulent kinetic energy, turbulence
!     frequency and their spatial derivatives for gases and liquids
!
               if(compressible.eq.1) then
!     
!     gas
!     
                  rho=0.d0
                  xkin=0.d0
                  xtuf=0.d0
                  do j1=1,3
                     dxkin(j1)=0.d0
                     dxtuf(j1)=0.d0
                  enddo
                  do i1=1,nope
                     rho=rho+shp(4,i1)*voldauxl(4,i1)
                     xkin=xkin+shp(4,i1)*voldtul(1,i1)
                     xtuf=xtuf+shp(4,i1)*voldtul(2,i1)
                     do j1=1,3
                        dxkin(j1)=dxkin(j1)+shp(j1,i1)*voldtul(1,i1)
                        dxtuf(j1)=dxtuf(j1)+shp(j1,i1)*voldtul(2,i1)
                     enddo
                  enddo
                  xkin=xkin/rho
                  xtuf=xtuf/rho
                  do j1=1,3
                     dxkin(j1)=dxkin(j1)/rho
                     dxtuf(j1)=dxtuf(j1)/rho
                  enddo
               else
!     
!     liquid
!     
                  call materialdata_rho(rhcon,nrhcon,imat,rho,
     &                 temp,ntmat_)
                  xkin=0.d0
                  xtuf=0.d0
                  do j1=1,3
                     dxkin(j1)=0.d0
                     dxtuf(j1)=0.d0
                  enddo
                  do i1=1,nope
                     xkin=xkin+shp(4,i1)*voldtul(1,i1)
                     xtuf=xtuf+shp(4,i1)*voldtul(2,i1)
                     do j1=1,3
                        dxkin(j1)=dxkin(j1)+shp(j1,i1)*voldtul(1,i1)
                        dxtuf(j1)=dxtuf(j1)+shp(j1,i1)*voldtul(2,i1)
                     enddo
                  enddo
                  xkin=xkin/rho
                  xtuf=xtuf/rho
                  do j1=1,3
                     dxkin(j1)=dxkin(j1)/rho
                     dxtuf(j1)=dxtuf(j1)/rho
                  enddo
               endif
!     
!     calculation of turbulent auxiliary variables
!     
!     vorticity
!     
               vort=dsqrt((vkl(3,2)-vkl(2,3))**2+
     &              (vkl(1,3)-vkl(3,1))**2+
     &              (vkl(2,1)-vkl(1,2))**2)
!     
!     kinematic viscosity
!     
               un=dvi/rho
!     
!     factor F2
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
!     
!     factor F1
!     
               if(turbulent.eq.1) then
!     
!     k-epsilon model
!     
                  f1=1.d0
               elseif(turbulent.eq.2) then
!     
!     q-omega model
!     
                  f1=0.d0
               else
!     
!     SST model
!     
                  cdktuf=max(2.d0*rho*stuf2*
     &                 (dxkin(1)*dxtuf(1)+dxkin(2)*dxtuf(2)+
     &                 dxkin(3)*dxtuf(3))/xtuf,1.d-20)
                  arg1=min(max(c1,c2),4.d0*rho*stuf2*xkin/(cdktuf*y*y))
                  f1=dtanh(arg1**4.d0)
               endif
               f1m=1.d0-f1
!     
!     interpolation of the constants
!     
               skin=f1*skin1+f1m*skin2
               stuf=f1*stuf1+f1m*stuf2
!     
!     auxiliary quantities
!     
               umsk=dvi+skin*umt
               umst=dvi+stuf*umt
!     
!     determining the stress and and stress x velocity + conductivity x
!     temperature gradient
!     
               do i1=1,3
                     tvk(i1)=umsk*dxkin(i1)
                     tvt(i1)=umsk*dxtuf(i1)
               enddo
!
               tvnk=tvk(1)*xsj2(1)+tvk(2)*xsj2(2)+tvk(3)*xsj2(3)
               tvnt=tvt(1)*xsj2(1)+tvt(2)*xsj2(2)+tvt(3)*xsj2(3)
!
               xsjmodk=tvnk*weight*dtime
               xsjmodt=tvnt*weight*dtime
               do k=1,nopes
                  if((nope.eq.20).or.(nope.eq.8)) then
                     ipointer=ifaceq(k,ig)
                  elseif((nope.eq.10).or.(nope.eq.4)) then
                     ipointer=ifacet(k,ig)
                  else
                     ipointer=ifacew(k,ig)
                  endif
                  ffk(ipointer)=ffk(ipointer)+shp2(4,k)*xsjmodk
                  fft(ipointer)=fft(ipointer)+shp2(4,k)*xsjmodt
               enddo
            enddo
            idf=idf-1
         enddo
      endif
      if(nelem.eq.1055) then
         write(*,*) 'e_c3d_krhs ',(ffk(j),j=1,nope)
         write(*,*) 'e_c3d_krhs ',(fft(j),j=1,nope)
      endif
!
      return
      end

