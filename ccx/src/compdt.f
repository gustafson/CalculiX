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
      subroutine compdt(nk,dt,nshcon,shcon,nrhcon,rhcon,vold,ntmat_,
     &  iponoel,inoel,dtimef,iexplicit,ielmat,physcon,dh,cocon,
     &  ncocon,ithermal,mi,ipkon,kon,lakon,dtl,ne,v,co,turbulent,voldtu)
!
!     - determine the time step for each node (stored in field dt
!       and the minimum value across all nodes (dtimef)
!
      implicit none
!
      character*8 lakon(*),lakonl
!
      integer nk,i,j,k,iponoel(*),inoel(3,*),index,nelem,ithermal,mi(*),
     &  nshcon(*),nrhcon(*),ntmat_,ielmat(mi(3),*),imat,ncocon(2,*),
     &  ipkon(*),kon(*),ne,nope,indexe,iflag,iexplicit,turbulent
!
      real*8 dtimef,dt(*),dvi,r,cp,rho,shcon(0:3,ntmat_,*),
     &  rhcon(0:1,ntmat_,*),vold(0:mi(2),*),temp,vel,dtu,dtnu,
     &  physcon(*),dh(*),cocon(0:6,ntmat_,*),dtal,cond,voldl(3,20),
     &  xl(3,20),vertex6(3,6),vertex8(3,8),xi,et,ze,xsj,shp(4,20),
     &  dtl(*),h,v(0:mi(2),*),co(3,*),dd,voldtu(2,*)
!
      data vertex6 /0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,
     &              0.d0,1.d0,0.d0,0.d0,0.d0,1.d0,
     &              1.d0,0.d0,1.d0,0.d0,1.d0,1.d0/
      data vertex8 /-1.d0,-1.d0,-1.d0,1.d0,-1.d0,-1.d0,
     &              1.d0,1.d0,-1.d0,-1.d0,1.d0,-1.d0,
     &              -1.d0,-1.d0,1.d0,1.d0,-1.d0,1.d0,
     &              1.d0,1.d0,1.d0,-1.d0,1.d0,1.d0/
      data iflag /3/
c!
c!     determining the element height in flow direction
c!
c      if(iexplicit.eq.1) then
c         do i=1,ne
c            indexe=ipkon(i)
c            if(indexe.lt.0) cycle
c            lakonl(1:8)=lakon(i)(1:8)
c!     
c!     number of nodes in the element
c!     
c            if(lakonl(4:4).eq.'2') then
c               nope=20
c            elseif(lakonl(4:4).eq.'8') then
c               nope=8
c            elseif(lakonl(4:5).eq.'10') then
c               nope=10
c            elseif(lakonl(4:4).eq.'4') then
c               nope=4
c            elseif(lakonl(4:5).eq.'15') then
c               nope=15
c            elseif(lakonl(4:4).eq.'6') then
c               nope=6
c            else
c               cycle
c            endif
c!     
c!     velocity at the nodes
c!     
c            do j=1,nope
c               do k=1,3
c                  voldl(k,j)=vold(k,kon(indexe+j))
c                  xl(k,j)=co(k,kon(indexe+j))
c               enddo
c            enddo
c!     
c!     element height
c!     
c            h=0.d0
c            do j=1,nope
c               if(nope.eq.20) then
c                  call shape20h(xi,et,ze,xl,xsj,shp,iflag)
c               elseif(nope.eq.8) then
c                  xi=vertex8(1,j)
c                  et=vertex8(2,j)
c                  ze=vertex8(3,j)
c                  call shape8h(xi,et,ze,xl,xsj,shp,iflag)
c               elseif(nope.eq.10) then
c                  call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
c               elseif(nope.eq.4) then
c                  call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
c               elseif(nope.eq.15) then
c                  call shape15w(xi,et,ze,xl,xsj,shp,iflag)
c               elseif(nope.eq.6) then
c                  xi=vertex6(1,j)
c                  et=vertex6(2,j)
c                  ze=vertex6(3,j)
c                  call shape6w(xi,et,ze,xl,xsj,shp,iflag)
c               endif
c!     
c               dd=dsqrt(voldl(1,j)*voldl(1,j)+
c     &              voldl(2,j)*voldl(2,j)+voldl(3,j)*voldl(3,j))
c               if(dd.lt.1.d-10) then
c                  cycle
c               else
c                  h=h+dabs(shp(1,j)*voldl(1,j)+shp(2,j)*voldl(2,j)+
c     &              shp(3,j)*voldl(3,j))/dd
c               endif
c            enddo
c!     
cc            if(h.gt.0.d0) h=2.d0/h
c            if(h.gt.0.d0) h=nope/h
c!     
c!        height at the nodes of the elements is replaced by the
c!        element height of the latter is smaller
c!
c            do j=1,nope
c               if(dtl(kon(indexe+j)).gt.h) dtl(kon(indexe+j))=h
c            enddo
c         enddo
c      endif
!     
!     determining the time increment dt for each node.
!
!     edge nodes (fields iponoel and inoel are determined in precfd.f)
!
      dtimef=1.d30
!
      do i=1,nk
         index=iponoel(i)
         if(index.le.0) cycle
!
!        look at an element belonging to the edge node
!
         nelem=inoel(1,index)
!     
!        determining the time increment
!
         imat=ielmat(1,nelem)
         temp=vold(0,i)
!
!        density for gases
!         
         vel=dsqrt(vold(1,i)**2+vold(2,i)**2+vold(3,i)**2)
         if(iexplicit.eq.1) then
            call materialdata_cp(imat,ntmat_,temp,shcon,nshcon,cp)
            r=shcon(3,1,imat)
            dt(i)=dh(i)/(dsqrt(cp*r*temp/(cp-r))+vel)
c            dtl(i)=dtl(i)/(dsqrt(cp*r*temp/(cp-r))+vel)
cstart shallow
cccc            dt(i)=dh(i)/(dsqrt(10.d0*rho)+vel)
cend shallow
c            if(dtl(i).lt.dtimef) dtimef=dtl(i)
            if(dt(i).lt.dtimef) dtimef=dt(i)
         else
            call materialdata_dvi(imat,ntmat_,temp,shcon,nshcon,dvi)
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
            if(vel.lt.1.d-10) vel=1.d-10
            dtu=dh(i)/vel
            if(dvi.lt.1.d-10) dvi=1.d-10
            dtnu=dh(i)*dh(i)*rho/(2.d0*dvi)
            dt(i)=dtu*dtnu/(dtu+dtnu)
            if(ithermal.gt.1) then
               call materialdata_cond(imat,ntmat_,temp,cocon,ncocon,
     &                cond)
               call materialdata_cp(imat,ntmat_,temp,shcon,nshcon,cp)
               if(cond.lt.1.d-10) cond=1.d-10
               dtal=dh(i)*dh(i)*rho*cp/(2.d0*cond)
               dt(i)=(dt(i)*dtal)/(dt(i)+dtal)
            endif
            if(turbulent.ne.0) then
               dtal=dh(i)*dh(i)*rho/
     &              (2.d0*(dvi+dabs(voldtu(1,i)/voldtu(2,i))))
               dt(i)=(dt(i)*dtal)/(dt(i)+dtal)
            endif
            if(dt(i).lt.dtimef) dtimef=dt(i)
         endif
!
      enddo
!
!     middle nodes (interpolation between neighboring end nodes;
!     still to be done)
!      
      return
      end
