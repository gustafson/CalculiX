!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine e_c3d_prhs(co,nk,konl,lakonl,sm,ff,nelem,nmethod,rhcon,
     &  nrhcon,ielmat,ntmat_,v,vold,voldcon,nelemface,sideface,nface,
     &  dtime,matname,mi,shcon,nshcon,theta1,physcon,
     &  iexplicit,ipvar,var,ipvarf,varf)
!
!     computation of the pressure element matrix and rhs for the element with
!     element with the topology in konl: step 2
!
!     sm: lhs matrix
!     ff: rhs 
!
      implicit none
!
      character*1 sideface(*)
      character*8 lakonl
      character*80 matname(*),amat
!
      integer konl(20),ifaceq(8,6),nelemface(*),nk,nelem,j1,index,
     &  nface,i,j,k,i1,i2,nmethod,ii,jj,id,ipointer,ig,kk,ipvar(*),
     &  mi(*),nrhcon(*),ielmat(mi(3),*),nshcon(*),ntmat_,nope,nopes,
     &  imat,mint2d,
     &  mint3d,ifacet(6,4),ifacew(8,5),iflag,iexplicit,
     &  ipvarf(*)
!
      real*8 co(3,*),shp(4,20),
     &  ff(60),theta1,c2i,xsjmod,rhcon(0:1,ntmat_,*),rhovel(3),
     &  shcon(0:3,ntmat_,*),vl(0:mi(2),20),xsj2(3),shp2(7,8),
     &  v(0:mi(2),*),xsj,sm(60,60),temp,voldcon(0:4,*),
     &  weight,vold(0:mi(2),*),delrhovel(3),aux(3),voldconl(0:4,20),
     &  dpress(3),voldconl2(0:4,8),voldl(0:mi(2),20),
     &  physcon(*),gg(60),divrhovel,auxg(3),var(*),varf(*)
!
      real*8 dtime
!
      include "gauss.f"
!
      ifaceq=reshape((/4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/),(/8,6/))
      ifacet=reshape((/1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/),(/6,4/))
      ifacew=reshape((/1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/),(/8,5/))
      iflag=3
!
      imat=ielmat(1,nelem)
      amat=matname(imat)
!     
      if(lakonl(4:4).eq.'2') then
         nope=20
         nopes=8
      elseif(lakonl(4:4).eq.'8') then
         nope=8
         nopes=4
      elseif(lakonl(4:5).eq.'10') then
         nope=10
         nopes=6
      elseif(lakonl(4:4).eq.'4') then
         nope=4
         nopes=3
      elseif(lakonl(4:5).eq.'15') then
         nope=15
      elseif(lakonl(4:4).eq.'6') then
         nope=6
      endif
!
      if(lakonl(4:5).eq.'8R') then
         mint2d=1
         mint3d=1
      elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) then
         mint2d=4
         mint3d=8
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
      elseif(lakonl(4:5).eq.'6 ') then
         mint3d=2
      elseif(lakonl(4:5).eq.'6R') then
         mint3d=1
      else
         mint3d=0
      endif
!     
!     initialisation for distributed forces
!     
      do i=1,nope
         ff(i)=0.d0
cc         gg(i)=0.d0
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
            voldconl(i2,i1)=voldcon(i2,konl(i1))
         enddo
c         voldl(0,i1)=vold(0,konl(i1))
         voldl(4,i1)=vold(4,konl(i1))
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
!         
c         index=index+nope+14
         index=index+nope+5
!
         xsjmod=dtime*xsj*weight
!     
!     calculating of
!     the temperature temp
!     the density times velocity rhovel
!     auxiliary variables delrhovel and dpress
!     
cc
cc         divrhovel=0.d0
cc
         do i1=1,3
            rhovel(i1)=0.d0
            delrhovel(i1)=0.d0
            dpress(i1)=0.d0
         enddo
         do i1=1,nope
            do j1=1,3
               rhovel(j1)=rhovel(j1)+shp(4,i1)*voldconl(j1,i1)
               delrhovel(j1)=delrhovel(j1)+shp(4,i1)*vl(j1,i1)
               dpress(j1)=dpress(j1)+shp(j1,i1)*voldl(4,i1)
cc
cc               divrhovel=divrhovel+shp(j1,i1)*voldconl(j1,i1)
cc
            enddo
         enddo
!
!        storing dpress
!
         do i1=1,3
            index=index+1
            var(index)=dpress(i1)
         enddo
         index=index+6
cc
cc         divrhovel=divrhovel*xsjmod
cc
         do j1=1,3
            aux(j1)=xsjmod*(rhovel(j1)+
     &            theta1*(delrhovel(j1)-dtime*dpress(j1)))
cc
cc            auxg(j1)=xsjmod*(
cc     &            theta1*(delrhovel(j1)-dtime*dpress(j1)))
cc
         enddo
!     
!     determination of lhs and rhs
!     
         do jj=1,nope
!     
            ff(jj)=ff(jj)+
     &         shp(1,jj)*aux(1)+shp(2,jj)*aux(2)+shp(3,jj)*aux(3)
ccd            gg(jj)=gg(jj)+
ccd     &         shp(1,jj)*auxg(1)+shp(2,jj)*auxg(2)+shp(3,jj)*auxg(3)
cc
cc            gg(jj)=gg(jj)+
cc     &         shp(1,jj)*auxg(1)+shp(2,jj)*auxg(2)+shp(3,jj)*auxg(3)-
cc     &         shp(4,jj)*divrhovel
cc
!
         enddo
!     
      enddo
!     
!     velocity normal to external surfaces
!
      if(nface.ne.0) then
         index=ipvarf(nelem)
         call nident(nelemface,nelem,nface,id)
         do
            if((id.eq.0).or.(nelemface(id).ne.nelem)) exit
c     read(sideface(id)(1:1),'(i1)') ig
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
                     voldconl2(j,i)=voldcon(j,konl(ifaceq(i,ig)))
                  enddo
               enddo
            elseif((nope.eq.10).or.(nope.eq.4)) then
               do i=1,nopes
                  do j=1,3
                     voldconl2(j,i)=voldcon(j,konl(ifacet(i,ig)))
                  enddo
               enddo
            else
               do i=1,nopes
                  do j=1,3
                     voldconl2(j,i)=voldcon(j,konl(ifacew(i,ig)))
                  enddo
               enddo
            endif
!     
            do i=1,mint2d
               if((lakonl(4:5).eq.'8R').or.
     &              ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
                  weight=weight2d1(i)
               elseif((lakonl(4:4).eq.'8').or.
     &                 (lakonl(4:6).eq.'20R').or.
     &                 ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
                  weight=weight2d2(i)
               elseif(lakonl(4:4).eq.'2') then
                  weight=weight2d3(i)
               elseif((lakonl(4:5).eq.'10').or.
     &                 ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
                  weight=weight2d5(i)
               elseif((lakonl(4:4).eq.'4').or.
     &                 ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
                  weight=weight2d4(i)
               endif
!     
               do i1=1,nopes
                  index=index+1
                  shp2(4,i1)=varf(index)
               enddo
               do i1=1,3
                  index=index+1
                  xsj2(i1)=varf(index)
               enddo
               index=index+4*nope+9
!
               xsjmod=dtime*weight
!
c            do i1=1,3
c               dpress(i1)=0.d0
c            enddo
c            do i1=1,nope
c               do j1=1,3
c                  dpress(j1)=dpress(j1)+shp(j1,i1)*voldl(4,i1)
c               enddo
c            enddo
!     
!     calculating rho*velocity after step 1 at the
!     integration point
!     
               do k=1,3
                  rhovel(k)=0.d0
                  do j=1,nopes
                     rhovel(k)=rhovel(k)+
     &                    (voldconl2(k,j))*shp2(4,j)
c     rhovel(k)=rhovel(k)+
c     &                (voldconl2(k,j)+theta1*vl2(k,j))*shp2(4,j)
                  enddo
c     rhovel(k)=rhovel(k)-dtime*theta1*dpress(k)
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
     &                 (rhovel(1)*xsj2(1)+rhovel(2)*xsj2(2)+
     &                 rhovel(3)*xsj2(3))*xsjmod
               enddo
            enddo
            id=id-1
         enddo
      endif
!     
      return
      end

