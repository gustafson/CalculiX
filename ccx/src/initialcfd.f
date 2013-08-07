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
      subroutine initialcfd(yy,nk,co,ne,ipkon,kon,lakon,x,y,z,xo,yo,zo,
     &  nx,ny,nz,isolidsurf,neighsolidsurf,xsolidsurf,dh,nshcon,shcon,
     &  nrhcon,rhcon,vold,vcon,ntmat_,iponoel,inoel,
     &  iexplicit,ielmat,nsolidsurf,turbulent,physcon,compressible,
     &  matname,inomat,vcontu,mi,euler,ithermal)
!
!     initial calculations for cfd applicatons:
!     - determine the distance from the nearest solid surface
!       (stored in yy)
!     - determing the distance from the nearest in-flow node
!       for solid surface nodes (stored in xsolidsurf)
!     - determine the adjacent element height for each node 
!       (stored in field dh)
!     - calculate the value of the conservative variables (vcon)
!
      implicit none
!
      character*8 lakon(*)
      character*80 matname(*)
!
      integer ne,ipkon(*),kon(*),indexe,ifaceq(8,6),ifacet(7,4),
     &  ifacew(8,5),kflag,isolidsurf(*),nsolidsurf,nope,node1,node2,
     &  nk,node,i,j,k,iponoel(*),inoel(3,*),nx(*),ny(*),index,nelem,
     &  nz(*),neighsolidsurf(*),kneigh,nodep(4),iplaneq(3,8),iplanet(4),
     &  iplanew(2,6),nshcon(*),nrhcon(*),ntmat_,neigh,nodel,ifacel,
     &  mi(*),ielmat(mi(3),*),imat,inomat(*),nonei20(3,12),nonei10(3,6),
     &  nonei15(3,9),euler,ithermal,iexplicit,turbulent,compressible
!
      real*8 x(*),y(*),z(*),xo(*),yo(*),zo(*),xsolidsurf(*),
     &  yy(*),co(3,*),dh(*),r,cp,rho,shcon(0:3,ntmat_,*),vcontu(2,*),
     &  rhcon(0:1,ntmat_,*),vold(0:mi(2),*),vcon(0:4,*),px,py,pz,
     &  a,b,c,d,temp,vel,dtu,dtnu,physcon(*),xtu,xkin,dvi
!
!     nodes belonging to the element faces
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,11,
     &             1,2,4,5,9,8,12,
     &             2,3,4,6,10,9,13,
     &             1,4,3,8,10,7,14/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
!
      data iplaneq /2,4,5,2,5,6,2,3,6,2,3,4,
     &              1,4,5,1,5,6,1,3,6,1,3,4/
      data iplanet /3,4,2,1/
      data iplanew /2,4,2,5,2,3,1,4,1,5,1,3/
!
      data nonei10 /5,1,2,6,2,3,7,3,1,8,1,4,9,2,4,10,3,4/
!
      data nonei15 /7,1,2,8,2,3,9,3,1,10,4,5,11,5,6,12,6,4,
     &  13,1,4,14,2,5,15,3,6/
!
      data nonei20 /9,1,2,10,2,3,11,3,4,12,4,1,
     &  13,5,6,14,6,7,15,7,8,16,8,5,
     &  17,1,5,18,2,6,19,3,7,20,4,8/
!
      if(turbulent.ne.0) then
!
!     determining the nearest solid boundary node
!
         kflag=2
         kneigh=1
!     
         do i=1,nk
            yy(i)=-1.d0
         enddo
!     
         do i=1,nsolidsurf
            node=isolidsurf(i)
            x(i)=co(1,node)
            y(i)=co(2,node)
            z(i)=co(3,node)
            xo(i)=x(i)
            yo(i)=y(i)
            zo(i)=z(i)
            nx(i)=i
            ny(i)=i
            nz(i)=i
         enddo
         call dsort(x,nx,nsolidsurf,kflag)
         call dsort(y,ny,nsolidsurf,kflag)
         call dsort(z,nz,nsolidsurf,kflag)
!     
         do node=1,nk
            index=iponoel(node)
            if(index.le.0) cycle
            px=co(1,node)
            py=co(2,node)
            pz=co(3,node)
!     
!     determining the neighboring solid surface node
!     
            call near3d(xo,yo,zo,x,y,z,nx,ny,nz,px,py,pz,
     &           nsolidsurf,neigh,kneigh)
!     
            neigh=isolidsurf(neigh)
            yy(node)=dsqrt((co(1,node)-co(1,neigh))**2+
     &           (co(2,node)-co(2,neigh))**2+
     &           (co(3,node)-co(3,neigh))**2)
         enddo
!     
!     determining the distance to the nearest in-flow node for
!     solid surface nodes (middle nodes do not count as valid
!     in-flow nodes)
!     
         do i=1,nsolidsurf
            node1=isolidsurf(i)
            node2=neighsolidsurf(i)
            xsolidsurf(i)=dsqrt((co(1,node1)-co(1,node2))**2+
     &           (co(2,node1)-co(2,node2))**2+
     &           (co(3,node1)-co(3,node2))**2)
c            write(*,*) 'xsolidsurf ',node1,node2,xsolidsurf(i)
         enddo
!
      endif
!     
!     determining the smallest element height dh for each node. This
!     is the minimum of the height of all elements to which the
!     node belongs
!
!     edge nodes (fields iponoel and inoel are determined in precfd.f)
!
      loop:do i=1,nk
         index=iponoel(i)
         if(index.le.0) cycle
         dh(i)=1.d30
!
!        loop over all elements belonging to the edge node:
!        determining the minimum height
!
         do
            if(index.le.0) exit
            nelem=inoel(1,index)
            nodel=inoel(2,index)
            indexe=ipkon(nelem)
            if((lakon(nelem)(4:4).eq.'2').or.
     &         (lakon(nelem)(4:4).eq.'8')) then
               do j=1,3
                  if(nodel.gt.8) cycle loop
                  ifacel=iplaneq(j,nodel)
                  do k=1,4
                     nodep(k)=kon(indexe+ifaceq(k,ifacel))
                  enddo
                  call plane4(co,i,nodep,a,b,c,d)
                  dh(i)=min(dh(i),-a*co(1,i)-b*co(2,i)-c*co(3,i)-d)
               enddo
            elseif((lakon(nelem)(4:4).eq.'4').or.
     &             (lakon(nelem)(4:5).eq.'10')) then
               if(nodel.gt.4) cycle loop
               ifacel=iplanet(nodel)
               do k=1,3
                  nodep(k)=kon(indexe+ifacet(k,ifacel))
               enddo
               call plane3(co,nodep,a,b,c,d)
               dh(i)=min(dh(i),-a*co(1,i)-b*co(2,i)-c*co(3,i)-d)
            else
               do j=1,2
                  if(nodel.gt.6) cycle loop
                  ifacel=iplanew(j,nodel)
                  if(ifacel.le.2) then
                     do k=1,3
                        nodep(k)=kon(indexe+ifacew(k,ifacel))
                     enddo
                     call plane3(co,nodep,a,b,c,d)
                  else
                     do k=1,4
                        nodep(k)=kon(indexe+ifacew(k,ifacel))
                     enddo
                     call plane4(co,i,nodep,a,b,c,d)
                  endif
                  dh(i)=min(dh(i),-a*co(1,i)-b*co(2,i)-c*co(3,i)-d)
               enddo
            endif
            index=inoel(3,index)
         enddo
!
      enddo loop
!
!     middle nodes (interpolation between neighboring end nodes)
!      
      do i=1,ne
         if(lakon(nelem)(1:1).ne.'F') cycle
         indexe=ipkon(nelem)
         if(indexe.lt.0) cycle
         if(lakon(nelem)(4:5).eq.'20') then
            do j=9,20
               node=kon(indexe+j)
               dh(node)=(dh(kon(indexe+nonei20(2,j-8)))+
     &              dh(kon(indexe+nonei20(3,j-8))))/2.d0
            enddo
         elseif(lakon(nelem)(4:5).eq.'10') then
            do j=5,10
               node=kon(indexe+j)
               dh(node)=(dh(kon(indexe+nonei10(2,j-4)))+
     &              dh(kon(indexe+nonei10(3,j-4))))/2.d0
            enddo
         elseif(lakon(nelem)(4:5).eq.'15') then
            do j=7,15
               node=kon(indexe+j)
               dh(node)=(dh(kon(indexe+nonei15(2,j-6)))+
     &              dh(kon(indexe+nonei15(3,j-6))))/2.d0
            enddo
         endif
      enddo
!
!     calculate auxiliary fields
!         
      do node=1,nk
         if(inomat(node).eq.0) cycle
         imat=inomat(node)
         temp=vold(0,node)
         call materialdata_cp_sec(imat,ntmat_,temp,shcon,nshcon,cp,
     &        physcon)
!     
!     different treatment for gases and liquids
!     
         if(compressible.eq.1) then
            r=shcon(3,1,imat)
            if(r.lt.1.d-10) then
               write(*,*) '*ERROR in initialcfd: specific gas '
               write(*,*) 'constant for material ',matname(imat)
               write(*,*) 'is close to zero; maybe it has'
               write(*,*) 'not been defined'
               stop
            endif
            if(vold(0,node)-physcon(1).le.1.d-10) then
               write(*,*) '*ERROR in initialcfd: absolute temperature '
               write(*,*) '       is nearly zero; maybe absolute zero '
               write(*,*) '       was wrongly defined or not defined'
               write(*,*) '       at all (*PHYSICAL CONSTANTS card)'
               stop
            endif
            if(vold(4,node).le.0.d0) then
               write(*,*) '*ERROR in inicialcfd: initial pressure'
               write(*,*) '       must be strictly positive'
               stop
            endif
            rho=vold(4,node)/(r*(vold(0,node)-physcon(1)))
            vcon(0,node)=rho*(cp*(temp-physcon(1))+
     &           (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &           /2.d0)-vold(4,node)
!
!           check for inviscous (= Euler) calculations
!
c            if(euler.eq.1) then
c               call materialdata_dvi(imat,ntmat_,temp,shcon,nshcon,dvi)
cc               if(dvi.gt.1.d-20) euler=0
c               euler=0
c            endif
cstart shallow
c            rho=dsqrt(vold(4,node)/5.d0+(0.005*co(1,node))**2)
c            vcon(0,node)=rho*(cp*(temp-physcon(1))+
c     &           (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
c     &           /2.d0)
cend shallow
         else
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
            if(ithermal.gt.1) then
               vcon(0,node)=rho*(cp*(temp-physcon(1))+
     &           (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &           /2.d0)
            endif
         endif
         vcon(4,node)=rho
         do k=1,3
            vcon(k,node)=rho*vold(k,node)
         enddo
      enddo
!
!     initial conditions for turbulence parameters:
!     freestream conditions
!     
      if(turbulent.ne.0) then
         if(dabs(physcon(8)).lt.1.d-20) then
            write(*,*) '*ERROR in initialcfd: typical length of the'
            write(*,*) '       computational domain is too small'
            write(*,*) '       use *VALUES AT INFINITY to define a'
            write(*,*) '       finite value'
            stop
         endif
         xtu=5.5d0*physcon(5)/physcon(8)
c         xkin=10.d0**(-3.5d0)*xtu
         xkin=10.d0**(-2.d0)*xtu
         do node=1,nk
            imat=inomat(node)
            if(imat.eq.0) cycle
            temp=vold(0,node)
            call materialdata_dvi(imat,ntmat_,temp,shcon,nshcon,dvi)
!     
!     density for gases
!     
            if(compressible.eq.1) then
               r=shcon(3,1,imat)
               rho=vold(4,node)/
     &              (r*(vold(0,node)-physcon(1)))
            else
               call materialdata_rho(rhcon,nrhcon,imat,rho,
     &              temp,ntmat_,ithermal)
            endif
!     
            vcontu(1,node)=xkin*dvi
            vcontu(2,node)=xtu*rho
         enddo
      endif
!     
      return
      end
