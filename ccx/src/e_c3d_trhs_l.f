!
!     CalculiX 3-dimensional finite element program
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
      subroutine e_c3d_trhs_l(co,nk,konl,lakonl,p1,p2,omx,bodyfx,
     &  nbody,ff,nelem,nmethod,rhcon,nrhcon,
     &  ielmat,ntmat_,vold,voldaux,nelemload,
     &  sideload,xload,nload,idist,dtl,matname,mi,
     &  ttime,time,istep,iinc,xloadold,reltimef,shcon,nshcon,cocon,
     &  ncocon,physcon,nelemface,sideface,nface,
     &  ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,compressible,v,
     &  voldtu,yy,turbulent)
!
!     computation of the energy element matrix and rhs for the element with
!     element with the topology in konl: step 4
!
!     ff: rhs 
!
      implicit none
!
      integer flux,compressible
!
      character*1 sideface(*)
      character*8 lakonl
      character*20 sideload(*)
      character*80 matname(*),amat
!
      integer konl(20),ifaceq(8,6),nelemload(2,*),nk,nbody,nelem,
     &  nload,idist,i,j,k,i1,i2,j1,ncocon(2,*),k1,node,nfield,
     &  nmethod,ii,jj,id,ipointer,ig,kk,nrhcon(*),ielmat(*),nshcon(*),
     &  ntmat_,nope,nopes,imat,mint2d,mint3d,mi(2),ifacet(6,4),nopev,
     &  ifacew(8,5),istep,iinc,layer,kspt,jltyp,iflag,nelemface(*),
     &  nface,igl,idf,ipompc(*),nodempc(3,*),nmpc,ikmpc(*),ilmpc(*),
     &  iscale,turbulent
!
      real*8 co(3,*),xl(3,20),shp(4,20),xs2(3,7),dvi,
     &  p1(3),p2(3),bodyf(3),bodyfx(3),ff(60),cond,enthalpy,
     &  bf(3),q(3),xsjmod,dtem(3),vkl(3,3),corio(3),sinktemp,
     &  rhcon(0:1,ntmat_,*),reltimef,t(3,3),tv(3),bfv,press,
     &  vel(3),div,shcon(0:3,ntmat_,*),pgauss(3),dxsj2,areaj,
     &  voldl(0:mi(2),20),xloadold(2,*),cocon(0:6,ntmat_,*),
     &  xl2(3,8),xsj2(3),shp2(7,8),vold(0:mi(2),*),xload(2,*),
     &  om,omx,xi,et,ze,const,xsj,field,physcon(*),tvn,
     &  temp,voldaux(0:4,*),voldauxl(0:4,20),rho,xi3d,et3d,ze3d,
     &  weight,shpv(20),xlocal20(3,9,6),coefmpc(*),v(0:mi(2),*),
     &  xlocal4(3,1,4),xlocal10(3,3,4),xlocal6(3,1,5),vl(0:mi(2),20),
     &  xlocal15(3,4,5),xlocal8(3,4,6),xlocal8r(3,1,6),omcor,
     &  shpvnithi(20),voldtu(2,*),voldtul(2,20),yy(*),yyl(20),
     &  y,xtuf,xkin,vort,un,unt,umt,f2,arg2,c1,c2,a1
!
      real*8 dtime,ttime,time,tvar(2),coords(3)
c
c        inhomogeneous dtlime
c
      real*8 dtl(*)
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
!     calculating a time increment based on the element length in
!     flow direction
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
     &         (lakonl(7:7).eq.'E')) then
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
         ff(i)=0.d0
      enddo
!
!     temperature, velocity and auxiliary variables
!     (rho*energy density, rho*velocity and rho)
!
         do i1=1,nope
            do i2=0,4
               voldl(i2,i1)=vold(i2,konl(i1))
               voldauxl(i2,i1)=voldaux(i2,konl(i1))
               vl(i2,i1)=v(i2,konl(i1))
            enddo
            if(turbulent.ne.0) then
               voldtul(1,i1)=voldtu(1,konl(i1))
               voldtul(2,i1)=voldtu(2,konl(i1))
               yyl(i1)=yy(konl(i1))
            endif
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
c         xsjmod=dtime*xsj*weight
         xsjmod=xsj*weight
!     
!        calculating of
!        the temperature temp
!        the total enthalpy
!        the velocity vel
!        the temperature gradient dtem
!        the velocity gradient vkl
!        the divergence of the velocity div 
!        the divergence of the shape function times the velocity shpv(*)
!              in the integration point
!     
         temp=0.d0
         enthalpy=0.d0
         do i1=1,3
            vel(i1)=0.d0
            dtem(i1)=0.d0
            do j1=1,3
               vkl(i1,j1)=0.d0
            enddo
         enddo
         do i1=1,nope
            temp=temp+shp(4,i1)*voldl(0,i1)
            do j1=1,3
               vel(j1)=vel(j1)+shp(4,i1)*voldl(j1,i1)
               dtem(j1)=dtem(j1)+shp(j1,i1)*voldl(0,i1)
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
     &                   shp(3,i1)*vel(3)+shp(4,i1)*div
                enthalpy=enthalpy+shpv(i1)*(voldauxl(0,i1)+voldl(4,i1))
         enddo
!     
!     material data (density, dynamic viscosity, heat capacity and
!     conductivity)
!     
         call materialdata_dvi(imat,ntmat_,temp,shcon,nshcon,dvi)
         call materialdata_cond(imat,ntmat_,temp,cocon,ncocon,cond)
!
!     determining the stress
!
         do i1=1,3
            do j1=i1,3
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
!     determining stress x velocity + conductivity x
!     temperature gradient
!
         do i1=1,3
            tv(i1)=dvi*(t(i1,1)*vel(1)+t(i1,2)*vel(2)+t(i1,3)*vel(3))+
     &             cond*dtem(i1)
         enddo
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
!     determination of the rhs of the energy equations
!     
         do jj=1,nope
            ff(jj)=ff(jj)-xsjmod*dtl(konl(jj))*(
     &           (shp(4,jj)+dtl(konl(jj))*shpv(jj)/2.d0)*enthalpy+
     &           shp(1,jj)*tv(1)+shp(2,jj)*tv(2)+shp(3,jj)*tv(3))
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
                  bf(i1)=bodyf(i1)+(q(i1)-const*p2(i1))*om+
     &                   corio(i1)*omcor
               enddo
            else
               do i1=1,3
                  bf(i1)=bodyf(i1)
               enddo
            endif
            bfv=bf(1)*vel(1)+bf(2)*vel(2)+bf(3)*vel(3)
            do jj=1,nope
               ff(jj)=ff(jj)+xsjmod*dtl(konl(jj))*(shp(4,jj)+
     &              dtl(konl(jj))*shpv(jj)/2.d0)*bfv
            enddo
         endif
!
!           distributed heat flux
!
         if(nload.gt.0) then
            call nident2(nelemload,nelem,nload,id)
            areaj=xsj*weight
            do
               if((id.eq.0).or.(nelemload(1,id).ne.nelem)) exit
               if(sideload(id)(1:2).ne.'BF') then
                  id=id-1
                  cycle
               endif
               if(sideload(id)(3:4).eq.'NU') then
                  do j=1,3
                     pgauss(j)=0.d0
                     do i1=1,nope
                        pgauss(j)=pgauss(j)+
     &                       shp(4,i1)*xl(j,i1)
                     enddo
                  enddo
                  jltyp=1
                  iscale=1
                  call dflux(xload(1,id),temp,istep,iinc,tvar,
     &                 nelem,kk,pgauss,jltyp,temp,press,sideload(id),
     &                 areaj,vold,co,lakonl,konl,ipompc,nodempc,coefmpc,
     &                 nmpc,ikmpc,ilmpc,iscale,mi)
                  if((nmethod.eq.1).and.(iscale.ne.0))
     &                  xload(1,id)=xloadold(1,id)+
     &                 (xload(1,id)-xloadold(1,id))*reltimef
               endif
               do jj=1,nope
                  ff(jj)=ff(jj)+xsjmod*dtl(konl(jj))*(shp(4,jj)+
     &              dtl(konl(jj))*shpv(jj)/2.d0)*xload(1,id)
               enddo
               exit
            enddo
         endif
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
            ig=ichar(sideface(idf)(1:1))-48
!
!           check for distributed flux
!
            flux=0
            call nident2(nelemload,nelem,nload,id)
            do
               if((id.eq.0).or.(nelemload(1,id).ne.nelem)) exit
               if((sideload(id)(1:1).ne.'F').and.
     &              (sideload(id)(1:1).ne.'R').and.
     &              (sideload(id)(1:1).ne.'S')) then
                  id=id-1
                  cycle
               endif
               igl=ichar(sideload(id)(2:2))-48
               if(igl.ne.ig) then
                  id=id-1
                  cycle
               endif
               flux=1
               exit
            enddo
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
!              calculating of
!              the temperature temp
!              the velocity vel
!              the temperature gradient dtem
!              the velocity gradient vkl
!              in the integration point
!     
               temp=0.d0
               do i1=1,3
                  vel(i1)=0.d0
                  dtem(i1)=0.d0
                  do j1=1,3
                     vkl(i1,j1)=0.d0
                  enddo
               enddo
               do i1=1,nope
                  temp=temp+shp(4,i1)*voldl(0,i1)
                  do j1=1,3
                     vel(j1)=vel(j1)+shp(4,i1)*voldl(j1,i1)
                     dtem(j1)=dtem(j1)+shp(j1,i1)*voldl(0,i1)
                     do k1=1,3
                        vkl(j1,k1)=vkl(j1,k1)+shp(k1,i1)*voldl(j1,i1)
                     enddo
                  enddo
               enddo
               if(compressible.eq.1) div=vkl(1,1)+vkl(2,2)+vkl(3,3)
!     
!     material data (density, dynamic viscosity, heat capacity and
!     conductivity)
!     
               call materialdata_dvi(imat,ntmat_,temp,shcon,nshcon,dvi)
               call materialdata_cond(imat,ntmat_,temp,cocon,ncocon,
     &                 cond)
!     
!     determining the stress and and stress x velocity + conductivity x
!     temperature gradient
!     
               do i1=1,3
                  do j1=1,3
                     t(i1,j1)=vkl(i1,j1)+vkl(j1,i1)
                  enddo
                  if(compressible.eq.1) t(i1,i1)=t(i1,i1)-2.d0*div/3.d0
                  tv(i1)=dvi*(t(i1,1)*vel(1)+t(i1,2)*vel(2)+
     &                 t(i1,3)*vel(3))
                  if(flux.eq.0) then
                     tv(i1)=tv(i1)+cond*dtem(i1)
                  endif
               enddo
!
               tvn=tv(1)*xsj2(1)+tv(2)*xsj2(2)+tv(3)*xsj2(3)
!
               if(flux.eq.1) then
                  dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &                 xsj2(3)*xsj2(3))
                  areaj=dxsj2*weight
                  sinktemp=xload(2,id)
!
!                 for nonuniform load: determine the coordinates of the
!                 point (transferred into the user subroutine)
!     
                  if((sideload(id)(3:4).eq.'NU').or.
     &                 (sideload(id)(5:6).eq.'NU')) then
                     do k=1,3
                        coords(k)=0.d0
                        do j=1,nopes
                           coords(k)=coords(k)+xl2(k,j)*shp2(4,j)
                        enddo
                     enddo
                     jltyp=ichar(sideload(id)(2:2))-48
                     jltyp=jltyp+10
                     if(sideload(id)(1:1).eq.'S') then
                        iscale=1
                        call dflux(xload(1,id),temp,istep,iinc,tvar,
     &                       nelem,i,coords,jltyp,temp,press,
     &                       sideload(id),areaj,vold,co,lakonl,konl,
     &                       ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,
     &                       iscale,mi)
                        if((nmethod.eq.1).and.(iscale.ne.0))
     &                        xload(1,id)=xloadold(1,id)+
     &                       (xload(1,id)-xloadold(1,id))*reltimef
                     elseif(sideload(id)(1:1).eq.'F') then
                        call film(xload(1,id),sinktemp,temp,istep,
     &                       iinc,tvar,nelem,i,coords,jltyp,field,
     &                       nfield,sideload(id),node,areaj,vold,mi)
                        if(nmethod.eq.1) xload(1,id)=xloadold(1,id)+
     &                       (xload(1,id)-xloadold(1,id))*reltimef
                     elseif(sideload(id)(1:1).eq.'R') then
                        call radiate(xload(1,id),xload(2,id),temp,istep,
     &                       iinc,tvar,nelem,i,coords,jltyp,field,
     &                       nfield,sideload(id),node,areaj,vold,mi)
                        if(nmethod.eq.1) xload(1,id)=xloadold(1,id)+
     &                       (xload(1,id)-xloadold(1,id))*reltimef
                     endif
                  endif
!
                  if(sideload(id)(1:1).eq.'S') then
!     
!     flux INTO the face is positive (input deck convention)
!     this is different from the convention in the theory
!     
                     tvn=tvn+xload(1,id)*dxsj2
                  elseif(sideload(id)(1:1).eq.'F') then
                     tvn=tvn-xload(1,id)*(temp-sinktemp)*dxsj2
                  elseif(sideload(id)(1:1).eq.'R') then
                     tvn=tvn-physcon(2)*
     &                    xload(1,id)*((temp-physcon(1))**4-
     &                    (xload(2,id)-physcon(1))**4)*dxsj2
                  endif
               endif
!     
c               xsjmod=tvn*weight*dtime
               xsjmod=tvn*weight
               do k=1,nopes
                  if((nope.eq.20).or.(nope.eq.8)) then
                     ipointer=ifaceq(k,ig)
                  elseif((nope.eq.10).or.(nope.eq.4)) then
                     ipointer=ifacet(k,ig)
                  else
                     ipointer=ifacew(k,ig)
                  endif
                  ff(ipointer)=ff(ipointer)+
     &                 shp2(4,k)*xsjmod*dtl(konl(ipointer))
               enddo
            enddo
            idf=idf-1
         enddo
      endif
!     
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

