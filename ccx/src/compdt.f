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
      subroutine compdt(nk,dt,nshcon,shcon,nrhcon,rhcon,vold,ntmat_,
     &  iponoel,inoel,dtimef,iexplicit,ielmat,physcon,dh,cocon,
     &  ncocon,ithermal,mi,ipkon,kon,lakon,dtl,ne,v,co)
!
!     - determine the time step for each node (stored in field dt
!       and the minimum value across all nodes (dtimef)
!
      implicit none
!
      character*8 lakon(*),lakonl
!
      integer nk,i,j,k,iponoel(*),inoel(3,*),index,nelem,ithermal,
     &  nshcon(*),nrhcon(*),ntmat_,ielmat(*),imat,ncocon(2,*),mi(2),
     &  ipkon(*),kon(*),ne,nope,indexe,iflag,iexplicit
!
      real*8 dtimef,dt(*),dvi,r,cp,rho,shcon(0:3,ntmat_,*),
     &  rhcon(0:1,ntmat_,*),vold(0:mi(2),*),temp,vel,dtu,dtnu,
     &  physcon(*),dh(*),cocon(0:6,ntmat_,*),dtal,cond,voldl(3,20),
     &  xl(3,20),vertex6(3,6),vertex8(3,8),xi,et,ze,xsj,shp(4,20),
     &  dtl(*),h,v(0:mi(2),*),co(3,*),dd
!
      data vertex6 /0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,
     &              0.d0,1.d0,0.d0,0.d0,0.d0,1.d0,
     &              1.d0,0.d0,1.d0,0.d0,1.d0,1.d0/
      data vertex8 /-1.d0,-1.d0,-1.d0,1.d0,-1.d0,-1.d0,
     &              1.d0,1.d0,-1.d0,-1.d0,1.d0,-1.d0,
     &              -1.d0,-1.d0,1.d0,1.d0,-1.d0,1.d0,
     &              1.d0,1.d0,1.d0,-1.d0,1.d0,1.d0/
      data iflag /3/
!
!     determining the element height in flow direction
!
      if(iexplicit.eq.1) then
         do i=1,ne
            indexe=ipkon(i)
            lakonl(1:8)=lakon(i)(1:8)
!     
!     number of nodes in the element
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
!     velocity at the nodes
!     
            do j=1,nope
               do k=1,3
                  voldl(k,j)=vold(k,kon(indexe+j))
                  xl(k,j)=co(k,kon(indexe+j))
               enddo
            enddo
!     
!     element height
!     
            h=0.d0
            do j=1,nope
               if(nope.eq.20) then
                  call shape20h(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.8) then
                  xi=vertex8(1,j)
                  et=vertex8(2,j)
                  ze=vertex8(3,j)
                  call shape8h(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.10) then
                  call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.4) then
                  call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.15) then
                  call shape15w(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.6) then
                  xi=vertex6(1,j)
                  et=vertex6(2,j)
                  ze=vertex6(3,j)
                  call shape6w(xi,et,ze,xl,xsj,shp,iflag)
               endif
!     
               dd=dsqrt(voldl(1,j)*voldl(1,j)+
     &              voldl(2,j)*voldl(2,j)+voldl(3,j)*voldl(3,j))
               if(dd.lt.1.d-10) then
                  cycle
               else
                  h=h+dabs(shp(1,j)*voldl(1,j)+shp(2,j)*voldl(2,j)+
     &              shp(3,j)*voldl(3,j))/dd
               endif
            enddo
!     
            if(h.gt.0.d0) h=2.d0/h
!     
!        height at the nodes of the elements is replaced by the
!        element height of the latter is smaller
!
            do j=1,nope
c               if(kon(indexe+j).eq.9363) write(*,*) kon(indexe+j),
c     &                h,dtl(kon(indexe+j))
               if(dtl(kon(indexe+j)).gt.h) dtl(kon(indexe+j))=h
            enddo
         enddo
c      write(*,*) 'node ',kon(ipkon(4471)+1),dh(kon(ipkon(4471)+1)),
c     &           dtl(kon(ipkon(4471)+1))
      endif
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
         imat=ielmat(nelem)
         temp=vold(0,i)
!
!        density for gases
!         
         vel=dsqrt(vold(1,i)**2+vold(2,i)**2+vold(3,i)**2)
         if(iexplicit.eq.1) then
            call materialdata_cp(imat,ntmat_,temp,shcon,nshcon,cp)
            r=shcon(3,1,imat)
            if(dtl(i).lt.dh(i)) dtl(i)=dh(i)
            dt(i)=dh(i)/(dsqrt(cp*r*temp/(cp-r))+vel)
            dtl(i)=dt(i)*dtl(i)/dh(i)
            if(dtl(i).lt.dtimef) dtimef=dtl(i)
         else
            call materialdata_dvi(imat,ntmat_,temp,shcon,nshcon,dvi)
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_)
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
            if(dt(i).lt.dtimef) dtimef=dt(i)
         endif
!
      enddo
!
!     middle nodes (interpolation between neighboring end nodes;
!     still to be done)
!      
!
!     determining the minimum height across the complete fluid mesh            
!
c      dtimef=1.d30
c      if(iexplicit.eq.1) then
c         do i=1,nk
c            if(dtl(i).lt.dtimef) dtimef=dtl(i)
c         enddo
c      else
c         do i=1,nk
c            if(dt(i).lt.dtimef) dtimef=dt(i)
c         enddo
c      endif
!
      return
      end
