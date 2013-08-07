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
      subroutine applybounk(nodeboun,ndirboun,nboun,xbounact,
     &  iponoel,vold,ipompc,nodempc,coefmpc,nmpc,nfreestream,
     &  ifreestream,nsolidsurf,isolidsurf,xsolidsurf,
     &  inoel,physcon,compressible,ielmat,nshcon,shcon,nrhcon,
     &  rhcon,vcontu,ntmat_,labmpc,inomat,mi,ithermal)
!
!     applies turbulence boundary conditions
!
      implicit none
!
      character*20 labmpc(*)
!
      integer turbulent,
     &  nodeboun(*),ndirboun(*),i,j,nboun,node,
     &  index,nodei,ndiri,ist,ipompc(*),nodempc(3,*),nmpc,
     &  ndir,nfreestream,ifreestream(*),iponoel(*),mi(*),
     &  inoel,imat,ielmat(mi(3),*),ntmat_,nshcon(*),nrhcon(*),
     &  compressible,
     &  nsolidsurf,isolidsurf(*),inomat(*),ithermal
!
      real*8 vold(0:mi(2),*),xbounact(*),residuk,size,coefmpc(*),
     &  xtu,xkin,temp,r,dvi,rho,physcon(*),shcon(0:3,ntmat_,*),
     &  rhcon(0:1,ntmat_,*),xsolidsurf(*),vcontu(2,*),residut,
     &  correctionk,correctiont
!
!     freestream conditions
!      
      xtu=5.5d0*physcon(5)/physcon(8)
c      xkin=10.d0**(-3.5d0)*xtu
      xkin=10.d0**(-2.d0)*xtu
      do j=1,nfreestream
         node=ifreestream(j)
c         write(*,*) 'applybounk  freestream ',node
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
     &        (r*(vold(0,node)-physcon(1)))
         else
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
         endif
!
         vcontu(1,node)=xkin*dvi
         vcontu(2,node)=xtu*rho
      enddo
!
!     solid boundary conditions  
!
      do j=1,nsolidsurf
         node=isolidsurf(j)
c         write(*,*) 'applybounk solid  ',node
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
     &        (r*(vold(0,node)-physcon(1)))
         else
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
         endif
!
         vcontu(1,node)=0.d0
         vcontu(2,node)=800.d0*dvi/(xsolidsurf(j)**2)
c         write(*,*) 'applybounk ',node,xsolidsurf(j)
      enddo
!
!     taking fluid pressure MPC's into account: it is assumed
!     that cyclic fluid pressure MPC's also apply to the turbulent
!     parameters
!
      do i=1,nmpc
         if(labmpc(i)(1:6).ne.'CYCLIC') cycle
         ist=ipompc(i)
         ndir=nodempc(2,ist)
         if(ndir.ne.4) cycle
         node=nodempc(1,ist)
c         write(*,*) 'applybounk  mpc ',node
!     
!     check whether fluid MPC
!     
         imat=inomat(node)
         if(imat.eq.0) cycle
c         index=iponoel(node)
c         if(index.le.0) cycle
!     
         index=nodempc(3,ist)
         residuk=coefmpc(ist)*vcontu(1,node)
         residut=coefmpc(ist)*vcontu(2,node)
         size=coefmpc(ist)**2
         if(index.ne.0) then
            do
               nodei=nodempc(1,index)
c               ndiri=nodempc(2,index)
!     
               residuk=residuk+coefmpc(index)*vcontu(1,nodei)
               residut=residut+coefmpc(index)*vcontu(2,nodei)
               size=size+coefmpc(index)**2
               index=nodempc(3,index)
               if(index.eq.0) exit
            enddo
         endif
!     
!     correcting all terms of the MPC
!     
         residuk=residuk/size
         residut=residut/size
!     
         correctionk=-residuk*coefmpc(ist)
         correctiont=-residut*coefmpc(ist)
         vcontu(1,node)=vcontu(1,node)+correctionk
         vcontu(2,node)=vcontu(2,node)+correctiont
         index=nodempc(3,ist)
         if(index.ne.0) then
            do
               nodei=nodempc(1,index)
c               ndiri=nodempc(2,index)
!     
               correctionk=-residuk*coefmpc(index)
               correctiont=-residut*coefmpc(index)
               vcontu(1,nodei)=vcontu(1,nodei)+correctionk
               vcontu(2,nodei)=vcontu(2,nodei)+correctiont
               index=nodempc(3,index)
               if(index.eq.0) exit
            enddo
         endif
      enddo
!
      return
      end
