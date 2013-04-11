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
      subroutine e_c3d_prhs(co,nk,konl,lakonl,sm,gg,nelem,nmethod,rhcon,
     &  nrhcon,ielmat,ntmat_,v,vold,voldaux,nelemface,sideface,nface,
     &  dtime,matname,mint_,shcon,nshcon,theta1,physcon)
!
!     computation of the pressure element matrix and rhs for the element with
!     element with the topology in konl: step 2
!
!     sm: lhs matrix
!     ff: rhs 
!
      implicit none
!
      logical liquid
!
      character*1 sideface(*)
      character*8 lakonl
      character*80 matname(*),amat
!
      integer konl(20),ifaceq(8,6),nelemface(*),nk,nelem,j1,
     &  nface,i,j,k,i1,i2,nmethod,ii,jj,id,ipointer,ig,kk,
     &  nrhcon(*),ielmat(*),nshcon(*),ntmat_,nope,nopes,imat,mint2d,
     &  mint3d,mint_,ifacet(6,4),nopev,ifacew(8,5),iflag
!
      real*8 co(3,*),xl(3,20),shp(4,20),xs2(3,2),r,cp,dvi,
     &  ff(60),theta1,c2i,xsjmod,rhcon(0:1,ntmat_,*),rhovel(3),
     &  shcon(0:3,ntmat_,*),vl(0:4,20),xl2(0:3,8),xsj2(3),shp2(4,8),
     &  v(0:4,*),xi,et,ze,xsj,sm(60,60),temp,voldaux(0:4,*),
     &  rho,weight,vold(0:4,*),delrhovel(3),aux(3),voldauxl(0:4,20),
     &  dpress(3),voldauxl2(0:4,8),vl2(0:4,8),voldl(0:4,20),
     &  physcon(3),gg(60),divrhovel,auxg(3)
!
      real*8 dtime
!
      include "gauss.f"
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
!     initialisation for distributed forces
!     
      do i=1,nope
         ff(i)=0.d0
         gg(i)=0.d0
      enddo
!     
!     temperature, velocity and auxiliary variables
!     (rho*energy density, rho*velocity and rho)
!     
      do i1=1,nope
         do i2=1,3
            vl(i2,i1)=v(i2,konl(i1))
         enddo
         do i2=1,3
            voldauxl(i2,i1)=voldaux(i2,konl(i1))
         enddo
         voldl(0,i1)=vold(0,konl(i1))
         voldl(4,i1)=vold(4,konl(i1))
      enddo
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
!     calculating of
!     the temperature temp
!     the density times velocity rhovel
!     auxiliary variables delrhovel and dpress
!     
         temp=0.d0
c         divrhovel=0.d0
         do i1=1,3
            rhovel(i1)=0.d0
            delrhovel(i1)=0.d0
            dpress(i1)=0.d0
         enddo
         do i1=1,nope
            temp=temp+shp(4,i1)*voldl(0,i1)
            do j1=1,3
               rhovel(j1)=rhovel(j1)+shp(4,i1)*voldauxl(j1,i1)
               delrhovel(j1)=delrhovel(j1)+shp(4,i1)*vl(j1,i1)
               dpress(j1)=dpress(j1)+shp(j1,i1)*voldl(4,i1)
c               divrhovel=divrhovel+shp(j1,i1)*voldauxl(j1,i1)
            enddo
         enddo
c         divrhovel=divrhovel*xsjmod
         do j1=1,3
            aux(j1)=xsjmod*(rhovel(j1)+
     &            theta1*(delrhovel(j1)-dtime*dpress(j1)))
            auxg(j1)=xsjmod*(
     &            theta1*(delrhovel(j1)-dtime*dpress(j1)))
         enddo
!     
!     determination of lhs and rhs
!     
         do jj=1,nope
!     
            ff(jj)=ff(jj)+
     &         shp(1,jj)*aux(1)+shp(2,jj)*aux(2)+shp(3,jj)*aux(3)
            gg(jj)=gg(jj)+
     &         shp(1,jj)*auxg(1)+shp(2,jj)*auxg(2)+shp(3,jj)*auxg(3)-
     &         shp(4,jj)*divrhovel
!
         enddo
!     
      enddo
!     
!     velocity normal to external surfaces
!
      call nident(nelemface,nelem,nface,id)
      do
         if((id.eq.0).or.(nelemface(id).ne.nelem)) exit
         read(sideface(id)(1:1),'(i1)') ig
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
                  vl2(j,i)=v(j,konl(ifaceq(i,ig)))
                  voldauxl2(j,i)=voldaux(j,konl(ifaceq(i,ig)))
               enddo
            enddo
         elseif((nope.eq.10).or.(nope.eq.4)) then
            do i=1,nopes
               do j=1,3
                  xl2(j,i)=co(j,konl(ifacet(i,ig)))
                  vl2(j,i)=v(j,konl(ifacet(i,ig)))
                  voldauxl2(j,i)=voldaux(j,konl(ifacet(i,ig)))
               enddo
            enddo
         else
            do i=1,nopes
               do j=1,3
                  xl2(j,i)=co(j,konl(ifacew(i,ig)))
                  vl2(j,i)=v(j,konl(ifacew(i,ig)))
                  voldauxl2(j,i)=voldaux(j,konl(ifacew(i,ig)))
               enddo
            enddo
         endif
!     
         do i=1,mint2d
            if((lakonl(4:5).eq.'8R').or.
     &           ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
               xi=gauss2d1(1,i)
               et=gauss2d1(2,i)
               weight=weight2d1(i)
            elseif((lakonl(4:4).eq.'8').or.
     &              (lakonl(4:6).eq.'20R').or.
     &              ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
               xi=gauss2d2(1,i)
               et=gauss2d2(2,i)
               weight=weight2d2(i)
            elseif(lakonl(4:4).eq.'2') then
               xi=gauss2d3(1,i)
               et=gauss2d3(2,i)
               weight=weight2d3(i)
            elseif((lakonl(4:5).eq.'10').or.
     &              ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
               xi=gauss2d5(1,i)
               et=gauss2d5(2,i)
               weight=weight2d5(i)
            elseif((lakonl(4:4).eq.'4').or.
     &              ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
               xi=gauss2d4(1,i)
               et=gauss2d4(2,i)
               weight=weight2d4(i)
            endif
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
            xsjmod=dtime*weight
!     
!     calculating rho*velocity after step 1 at the
!     integration point
!     
            do k=1,3
               rhovel(k)=0.d0
               do j=1,nopes
cc                  rhovel(k)=rhovel(k)+
cc     &                (voldauxl2(k,j)+vl2(k,j))*shp2(4,j)
c                  rhovel(k)=rhovel(k)+
c     &                (voldauxl2(k,j))*shp2(4,j)
                  rhovel(k)=rhovel(k)+
     &                (vl2(k,j))*shp2(4,j)
               enddo
            enddo
!     
            do k=1,nopes
               if((nope.eq.20).or.(nope.eq.8)) then
                  ipointer=ifaceq(k,ig)
               elseif((nope.eq.10).or.(nope.eq.4)) then
                  ipointer=ifacet(k,ig)
               else
                  ipointer=ifacew(k,ig)
               endif
               ff(ipointer)=ff(ipointer)-shp2(4,k)*
     &              (rhovel(1)*xsj2(1)+rhovel(2)*xsj2(2)+
     &               rhovel(3)*xsj2(3))*xsjmod
               gg(ipointer)=gg(ipointer)-shp2(4,k)*
     &              (rhovel(1)*xsj2(1)+rhovel(2)*xsj2(2)+
     &               rhovel(3)*xsj2(3))*xsjmod
            enddo
         enddo
         id=id-1
      enddo
!     
!     
!     for axially symmetric and plane stress/strain elements: 
!     complete s and sm
!     
c     if((lakonl(6:7).eq.'RA').or.(lakonl(6:7).eq.'RS').or.
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
c         do i=1,60
c            do j=i,60
c               k=abs(iperm(i))
c               l=abs(iperm(j))
c               if(k.gt.l) then
c                  m=k
c                  k=l
c                  l=m
c               endif
c               sax(i,j)=sm(k,l)*iperm(i)*iperm(j)/(k*l)
c            enddo
c         enddo
c         do i=1,60
c            do j=i,60
c               sm(i,j)=sm(i,j)+sax(i,j)
c            enddo
c         enddo
c      endif
c      if(nelem.eq.145) then
c         do jj=1,nope
c            write(*,*) 'ec3dprhs 2',jj,konl(jj),ff(jj),gg(jj)
c         enddo
c      endif
!
      return
      end

