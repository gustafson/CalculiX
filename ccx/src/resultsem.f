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
      subroutine resultsem(co,kon,ipkon,lakon,v,elcon,nelcon,ielmat,
     &  ntmat_,vold,dtime,matname,mi,ncmat_,nea,neb,sti,eei,alcon,
     &  nalcon,h0)
!
!     calculates the heat flux and the material tangent at the integration
!     points and the internal concentrated flux at the nodes
!
      implicit none
!
      character*8 lakon(*)
      character*80 amat,matname(*)
!
      integer kon(*),konl(26),mi(*),nelcon(2,*),ielmat(mi(3),*),
     &  ntmat_,ipkon(*),null,three,iflag,mt,i,j,k,m1,kk,i1,m3,indexe,
     &  nope,imat,mint3d,ncmat_,nea,neb,nalcon(2,*)
!
      real*8 co(3,*),v(0:mi(2),*),shp(4,26),xl(3,26),vl(0:mi(2),26),
     &  elcon(0:ncmat_,ntmat_,*),vkl(0:mi(2),3),vold(0:mi(2),*),
     &  elconloc(21),xi,et,ze,xsj,t1l,dtime,weight,
     &  voldkl(0:mi(2),3),alpha(6),
     &  h0l(3),al(3),aoldl(3),sti(6,mi(1),*),eei(6,mi(1),*),
     &  um,alcon(0:6,ntmat_,*),h0(3,*),voldl(0:mi(2),26)
!
      include "gauss.f"
!
      iflag=3
      null=0
      three=3
!
      mt=mi(2)+1
!
!     calculation of temperatures and thermal flux
!
      do i=nea,neb
!
         if(ipkon(i).lt.0) cycle
!
         imat=ielmat(1,i)
         amat=matname(imat)
!
         indexe=ipkon(i)
         if(lakon(i)(4:5).eq.'20') then
            nope=20
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:5).eq.'15') then
            nope=15
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         else
            cycle
         endif
!
         if(lakon(i)(4:5).eq.'8R') then
            mint3d=1
         elseif((lakon(i)(4:4).eq.'8').or.
     &          (lakon(i)(4:6).eq.'20R')) then
            mint3d=8
         elseif(lakon(i)(4:4).eq.'2') then
            mint3d=27
         elseif(lakon(i)(4:5).eq.'10') then
            mint3d=4
         elseif(lakon(i)(4:4).eq.'4') then
            mint3d=1
         elseif(lakon(i)(4:5).eq.'15') then
            mint3d=9
         elseif(lakon(i)(4:4).eq.'6') then
            mint3d=2
         endif
!
         do j=1,nope
            konl(j)=kon(indexe+j)
            do k=1,3
               xl(k,j)=co(k,konl(j))
            enddo
            do k=0,5
               vl(k,j)=v(k,konl(j))
            enddo
            voldl(4,j)=vold(4,konl(j))
         enddo
!
         do kk=1,mint3d
            if(lakon(i)(4:5).eq.'8R') then
               xi=gauss3d1(1,kk)
               et=gauss3d1(2,kk)
               ze=gauss3d1(3,kk)
               weight=weight3d1(kk)
            elseif((lakon(i)(4:4).eq.'8').or.
     &             (lakon(i)(4:6).eq.'20R'))
     &        then
               xi=gauss3d2(1,kk)
               et=gauss3d2(2,kk)
               ze=gauss3d2(3,kk)
               weight=weight3d2(kk)
            elseif(lakon(i)(4:4).eq.'2') then
               xi=gauss3d3(1,kk)
               et=gauss3d3(2,kk)
               ze=gauss3d3(3,kk)
               weight=weight3d3(kk)
            elseif(lakon(i)(4:5).eq.'10') then
               xi=gauss3d5(1,kk)
               et=gauss3d5(2,kk)
               ze=gauss3d5(3,kk)
               weight=weight3d5(kk)
            elseif(lakon(i)(4:4).eq.'4') then
               xi=gauss3d4(1,kk)
               et=gauss3d4(2,kk)
               ze=gauss3d4(3,kk)
               weight=weight3d4(kk)
            elseif(lakon(i)(4:5).eq.'15') then
               xi=gauss3d8(1,kk)
               et=gauss3d8(2,kk)
               ze=gauss3d8(3,kk)
               weight=weight3d8(kk)
            elseif(lakon(i)(4:4).eq.'6') then
               xi=gauss3d7(1,kk)
               et=gauss3d7(2,kk)
               ze=gauss3d7(3,kk)
               weight=weight3d7(kk)
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
!                 vkl(m2,m3) contains the derivative of the m2-
!                 component of the displacement with respect to
!                 direction m3
!     
            do k=1,5
               do m3=1,3
                  vkl(k,m3)=0.d0
               enddo
!     
               do m1=1,nope
                  do m3=1,3
                     vkl(k,m3)=vkl(k,m3)+shp(m3,m1)*vl(k,m1)
                  enddo
               enddo
            enddo
!     
            do m1=1,nope
               do m3=1,3
                  voldkl(4,m3)=vkl(4,m3)+shp(m3,m1)*vl(4,m1)
               enddo
            enddo
!
!              calculating the temperature difference in
!              the integration point
!
            t1l=0.d0
            do j=1,3
               h0l(j)=0.d0
               al(j)=0.d0
               aoldl(j)=0.d0
            enddo
            if(lakon(i)(4:5).eq.'8 ') then
               do i1=1,8
                  do j=1,3
                     h0l(j)=h0l(j)+h0(j,konl(i1))/8.d0
                     al(j)=al(j)+v(j,konl(i1))/8.d0
                     aoldl(j)=aoldl(j)+vold(j,konl(i1))/8.d0
                  enddo
                  t1l=t1l+v(0,konl(i1))/8.d0
               enddo
            elseif(lakon(i)(4:6).eq.'20 ') then
               call linscal(v,konl,nope,kk,t1l,mi(2))
               call linvec(h0,konl,nope,kk,h0l,null,three)
               call linvec(v,konl,nope,kk,al,null,mi(2))
               call linvec(vold,konl,nope,kk,aoldl,null,mi(2))
            else
               do i1=1,nope
                  t1l=t1l+shp(4,i1)*v(0,konl(i1))
                  do j=1,3
                     h0l(j)=h0l(j)+shp(4,i1)*h0(j,konl(i1))
                     al(j)=al(j)+shp(4,i1)*v(j,konl(i1))
                     aoldl(j)=aoldl(j)+shp(4,i1)*vold(j,konl(i1))
                  enddo
               enddo
            endif
!
!                 material data (permeability)
!
            call materialdata_em(elcon,nelcon,alcon,nalcon,
     &           imat,ntmat_,t1l,elconloc,ncmat_,alpha)
!
            um=elconloc(1)
!
            if(int(elconloc(2)).eq.1) then
!
!              magnetic field in phi-domain
!
               do k=1,3
                  sti(k,kk,i)=um*(h0l(k)-vkl(5,k))
               enddo
            else
!
!              magnetic field in A and A-V domain
!
               sti(1,kk,i)=vkl(3,2)-vkl(2,3)
               sti(2,kk,i)=vkl(1,3)-vkl(3,1)
               sti(3,kk,i)=vkl(2,1)-vkl(1,2)
!
!              electric intensity in A-V domain
!
               if(int(elconloc(2)).eq.2) then
                  do k=1,3
                     eei(k,kk,i)=(aoldl(k)-al(k)+
     &                            voldl(4,k)-vl(4,k))/dtime
                  enddo
               endif
            endif
!
         enddo
      enddo
!
      return
      end
