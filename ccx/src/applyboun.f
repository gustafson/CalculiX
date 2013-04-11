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
      subroutine applyboun(nodeboun,ndirboun,nboun,xbounact,
     &  ithermal,nk,iponoel,inoel,vold,vcontu,t1act,isolidsurf,
     &  nsolidsurf,xsolidsurf,nfreestream,ifreestream,turbulent,
     &  voldcon,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
     &  compressible,ismooth,nmpc,nodempc,ipompc,coefmpc,inomat,
     &  mi,ikboun,ilboun,ilmpc,labmpc)
!
!     applies temperature and velocity boundary conditions for
!     incompressible fluids (liquids)
!
      implicit none
!
      character*20 labmpc(*)
!
      integer turbulent,compressible
!
      integer nrhcon(*),mi(*),ielmat(mi(3),*),ntmat_,nodeboun(*),
     &  isolidsurf(*),j,ilboun(*),ilmpc(*),ikboun(*),
     &  ndirboun(*),nshcon(*),nk,i,nboun,node,imat,ithermal,iponoel(*),
     &  inoel(3,*),nsolidsurf,ifreenode,ifreestream(*),nfreestream,k,
     &  index,ismooth,indexi,nodei,nmpc,nodempc(3,*),ipompc(*),
     &  ist,ndir,ndiri,inomat(*),nref
!
      real*8 rhcon(0:1,ntmat_,*),rho,vold(0:mi(2),*),xbounact(*),shcon,
     &  vcontu(2,*),t1act(*),temp,xsolidsurf(*),reflength,sum,
     &  refkin,reftuf,refvel,voldcon(0:4,*),physcon(*),v(0:mi(2),*),
     &  rhoi,coefmpc(*),residu,size,correction,cp,xkin,xtu,sumk,sumt,
     &  dvi
!
!     inserting velocity and temperature SPC's for incompressible
!     materials (the pressure BC's have been applied in applybounp.f)
!
      nref=0
!
      do j=1,nboun
!
!        monotonically increasing DOF-order
!
         i=ilboun(j)
!
!        check whether temperature or velocity DOF
!
         ndir=ndirboun(i)
         if(ndir.gt.3) cycle
!
!        check whether fluid node
!
         node=nodeboun(i)
         imat=inomat(node)
         if(imat.eq.0) cycle
!
!        calculating the conservative variables for the previously
!        treated node, if any
!
         if(node.ne.nref) then
            if(nref.ne.0) then
               temp=vold(0,nref)
               call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
               call materialdata_cp_sec(imat,ntmat_,temp,shcon,nshcon,
     &              cp,physcon)
               if(ithermal.gt.1)
     &              voldcon(0,nref)=rho*(cp*(temp-physcon(1))+
     &              (vold(1,nref)**2+vold(2,nref)**2+vold(3,nref)**2)
     &              /2.d0)
               do k=1,3
                  voldcon(k,nref)=rho*vold(k,nref)
               enddo
c               voldcon(4,nref)=rho
            endif
            nref=node
         endif
!
!        calculating the physical variables for the node at stake
!
         vold(ndir,node)=xbounact(i)
      enddo
!
!     treating the remaining node
!
      if(nref.ne.0) then
         temp=vold(0,nref)
         call materialdata_rho(rhcon,nrhcon,imat,rho,
     &        temp,ntmat_,ithermal)
         if(ithermal.gt.1)
     &        voldcon(0,nref)=rho*(cp*(temp-physcon(1))+
     &        (vold(1,nref)**2+vold(2,nref)**2+vold(3,nref)**2)
     &        /2.d0)
         do k=1,3
            voldcon(k,nref)=rho*vold(k,nref)
         enddo
c         voldcon(4,nref)=rho
      endif
!
!     inserting velocity and temperature MPC's
!
      nref=0
!
      do j=1,nmpc
!
!        monotonically increasing DOF-order
!
         i=ilmpc(j)
         ist=ipompc(i)
!
!        check whether temperature or velocity DOF
!
         ndir=nodempc(2,ist)
         if(ndir.gt.3) cycle
!
!        check whether fluid node
!
         node=nodempc(1,ist)
         imat=inomat(node)
         if(imat.eq.0) cycle
!
!        calculating the value of the dependent DOF of the MPC
!
         index=nodempc(3,ist)
         if(index.eq.0) cycle
         sum=0.d0
         do
            nodei=nodempc(1,index)
            ndiri=nodempc(2,index)
            sum=sum+coefmpc(index)*vold(ndiri,nodei)
            index=nodempc(3,index)
            if(index.eq.0) exit
         enddo
!
!        calculating the conservative variables for the previously
!        treated node, if any
!
         if(node.ne.nref) then
            if(nref.ne.0) then
               temp=vold(0,nref)
               call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
               if(ithermal.gt.1)
     &              voldcon(0,nref)=rho*(cp*(temp-physcon(1))+
     &              (vold(1,nref)**2+vold(2,nref)**2+vold(3,nref)**2)
     &              /2.d0)
               do k=1,3
                  voldcon(k,nref)=rho*vold(k,nref)
               enddo
c               voldcon(4,nref)=rho
            endif
            nref=node
         endif
!
!        calculating the physical variables for the node at stake
!
         vold(ndir,node)=-sum/coefmpc(ist)
      enddo
!
!     treating the remaining node
!
      if(nref.ne.0) then
         temp=vold(0,nref)
         call materialdata_rho(rhcon,nrhcon,imat,rho,
     &        temp,ntmat_,ithermal)
         if(ithermal.gt.1)
     &        voldcon(0,nref)=rho*(cp*(temp-physcon(1))+
     &        (vold(1,nref)**2+vold(2,nref)**2+vold(3,nref)**2)
     &        /2.d0)
         do k=1,3
            voldcon(k,nref)=rho*vold(k,nref)
         enddo
c         voldcon(4,nref)=rho
      endif
!
!     freestream conditions for the turbulent variables
!      
      xtu=5.5d0*physcon(5)/physcon(8)
c      xkin=10.d0**(-3.5d0)*xtu
      xkin=10.d0**(-2.d0)*xtu
      do j=1,nfreestream
         node=ifreestream(j)
         imat=inomat(node)
         if(imat.eq.0) cycle
         temp=vold(0,node)
         call materialdata_dvi(imat,ntmat_,temp,shcon,nshcon,dvi)
!     
!        density 
!     
         if(compressible.eq.1) then
c            r=shcon(3,1,imat)
c            rho=vold(4,node)/
c     &        (r*(vold(0,node)-physcon(1)))
         else
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
         endif
!
         vcontu(1,node)=xkin*dvi
         vcontu(2,node)=xtu*rho
c         write(*,*) 'applyboun freestream ',node,xkin*dvi/rho,xtu
      enddo
!
!     solid boundary conditions for the turbulent variables 
!
      do j=1,nsolidsurf
         node=isolidsurf(j)
         imat=inomat(node)
         if(imat.eq.0) cycle
         temp=vold(0,node)
         call materialdata_dvi(imat,ntmat_,temp,shcon,nshcon,dvi)
!     
!        density
!     
         if(compressible.eq.1) then
c            r=shcon(3,1,imat)
c            rho=vold(4,node)/
c     &        (r*(vold(0,node)-physcon(1)))
         else
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
         endif
!
         vcontu(1,node)=0.d0
         vcontu(2,node)=800.d0*dvi/(xsolidsurf(j)**2)
c         write(*,*) 'applyboun solid boundary ',node,vcontu(1,node),
c     &            vcontu(2,node)/rho
      enddo
!
!     taking fluid pressure MPC's into account: it is assumed
!     that cyclic fluid pressure MPC's also apply to the turbulent
!     conservative variables (to be changed later -> physical
!     variables)
!
      do i=1,nmpc
         if(labmpc(i)(1:6).ne.'CYCLIC') cycle
         ist=ipompc(i)
         ndir=nodempc(2,ist)
         if(ndir.ne.4) cycle
         node=nodempc(1,ist)
!     
!     check whether fluid MPC
!     
         imat=inomat(node)
         if(imat.eq.0) cycle
!     
         index=nodempc(3,ist)
         sumk=0.d0
         sumt=0.d0
!
         do
            if(index.eq.0) exit
            sumk=sumk+coefmpc(index)*vcontu(1,nodei)
            sumt=sumt+coefmpc(index)*vcontu(2,nodei)
            index=nodempc(3,index)
         enddo
         vcontu(1,node)=-sumk/coefmpc(ist)
         vcontu(2,node)=-sumt/coefmpc(ist)
      enddo
!     
      return
      end
      
