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
      subroutine applyboun(nodeboun,ndirboun,nboun,xbounact,
     &  ithermal,nk,iponoel,inoel,vold,vcontu,t1act,isolidsurf,
     &  nsolidsurf,xsolidsurf,nfreestream,ifreestream,turbulent,
     &  vcon,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
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
      integer turbulent,compressible,
     &  nrhcon(*),mi(*),ielmat(mi(3),*),ntmat_,nodeboun(*),
     &  isolidsurf(*),j,ilboun(*),ilmpc(*),ikboun(*),
     &  ndirboun(*),nshcon(*),nk,i,nboun,node,imat,ithermal,iponoel(*),
     &  inoel(3,*),nsolidsurf,ifreenode,ifreestream(*),nfreestream,k,
     &  index,ismooth,indexi,nodei,nmpc,nodempc(3,*),ipompc(*),
     &  ist,ndir,ndiri,inomat(*),nref
!
      real*8 rhcon(0:1,ntmat_,*),rho,vold(0:mi(2),*),xbounact(*),
     &  shcon(0:3,ntmat_,*),
     &  vcontu(2,*),t1act(*),temp,xsolidsurf(*),reflength,sum,
     &  refkin,reftuf,refvel,vcon(0:4,*),physcon(*),v(0:mi(2),*),
     &  rhoi,coefmpc(*),residu,size,correction,cp,xkin,xtu,sumk,sumt,
     &  dvi,r
!
!     inserting velocity and temperature SPC's for incompressible
!     materials (the pressure BC's have been applied in applybounp.f)
!
!     inserting velocity, pressure and temperature SPC's for compressible
!     materials
!
      nref=0
!
      do j=1,nboun
!
!        monotonically increasing DOF-order
!
         i=ilboun(j)
!
!        pressure boundary conditions for incompressible materials
!        have already been treated
!
         ndir=ndirboun(i)
         if((compressible.eq.0).and.(ndir.gt.3)) cycle
!
!        check whether fluid node
!
         node=nodeboun(i)
         if(inomat(node).eq.0) cycle
!
!        calculating the conservative variables for the previously
!        treated node, if any
!
         if(node.ne.nref) then
            if(nref.ne.0) then
               call phys2con(inomat,nref,vold,ntmat_,shcon,nshcon,
     &              physcon,compressible,vcon,rhcon,nrhcon,ithermal,mi)
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
         call phys2con(inomat,nref,vold,ntmat_,shcon,nshcon,
     &        physcon,compressible,vcon,rhcon,nrhcon,ithermal,mi)
      endif
!
!     inserting velocity and temperature MPC's: incompressible fluids
!
      if((compressible.eq.0).or.(compressible.eq.1)) then
         nref=0
!     
         do j=1,nmpc
!     
!     monotonically increasing DOF-order
!     
            i=ilmpc(j)
            ist=ipompc(i)
!     
!     pressure multiple point constraints for incompressible materials
!     have already been treated
!     
            ndir=nodempc(2,ist)
            if((compressible.eq.0).and.(ndir.gt.3)) cycle
!     
!     check whether fluid node
!     
            node=nodempc(1,ist)
            if(inomat(node).eq.0) cycle
!     
!     calculating the value of the dependent DOF of the MPC
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
!     calculating the conservative variables for the previously
!     treated node, if any
!     
            if(node.ne.nref) then
               if(nref.ne.0) then
                  call phys2con(inomat,nref,vold,ntmat_,shcon,nshcon,
     &              physcon,compressible,vcon,rhcon,nrhcon,ithermal,mi)
               endif
               nref=node
            endif
!     
!     calculating the physical variables for the node at stake
!     
            vold(ndir,node)=-sum/coefmpc(ist)
         enddo
!
!     treating the remaining node
!     
         if(nref.ne.0) then
            call phys2con(inomat,nref,vold,ntmat_,shcon,nshcon,
     &           physcon,compressible,vcon,rhcon,nrhcon,ithermal,mi)
         endif
      else
!     
!     MPC's are treated by distributing the residual proportional to
!     the coefficients
!     
!     Right now it is assumed that the MPC's are independent of each 
!     other, i.e. degrees of freedom used in one MPC are not used in 
!     any other MPC
!     
         do i=1,nmpc
            index=ipompc(i)
!     
!     pressure multiple point constraints for incompressible materials
!     have already been treated
!     
            ndir=nodempc(2,index)
            if((compressible.eq.0).and.(ndir.gt.3)) cycle
!     
!     check whether fluid node
!     
            node=nodempc(1,index)
            if(inomat(node).eq.0) cycle
!     
!     calculating the value of the dependent DOF of the MPC
!     
            residu=coefmpc(index)*vold(ndir,node)
            size=coefmpc(index)**2
            if(index.eq.0) cycle
            do
               index=nodempc(3,index)
               if(index.eq.0) exit
               nodei=nodempc(1,index)
               ndiri=nodempc(2,index)
               residu=residu+coefmpc(index)*vold(ndiri,nodei)
               size=size+coefmpc(index)**2
            enddo
!     
!     correcting all terms of the MPC
!     
            residu=residu/size
!     
            index=ipompc(i)
            do
               nodei=nodempc(1,index)
               ndiri=nodempc(2,index)
!     
               correction=-residu*coefmpc(index)
               vold(ndiri,nodei)=vold(ndiri,nodei)+correction
               index=nodempc(3,index)
               if(index.eq.0) exit
            enddo
         enddo
!     
!     converting all nodes to conservative variables if mpc's were
!     active
!     
         if(nmpc.gt.0) then
            do i=1,nk
               if(inomat(i).eq.0) cycle
               call phys2con(inomat,nref,vold,ntmat_,shcon,nshcon,
     &              physcon,compressible,vcon,rhcon,nrhcon,ithermal,mi)
            enddo
         endif
!     
      endif
!     
      if(turbulent.eq.1) then
!     
!     freestream conditions for the turbulent variables
!     
         xtu=5.5d0*physcon(5)/physcon(8)
c     xkin=10.d0**(-3.5d0)*xtu
         xkin=10.d0**(-2.d0)*xtu
         do j=1,nfreestream
            node=ifreestream(j)
            imat=inomat(node)
            if(imat.eq.0) cycle
            temp=vold(0,node)
            call materialdata_dvi(imat,ntmat_,temp,shcon,nshcon,dvi)
!     
!     density 
!     
            if(compressible.eq.1) then
c     r=shcon(3,1,imat)
c     rho=vold(4,node)/
c     &        (r*(vold(0,node)-physcon(1)))
            else
               call materialdata_rho(rhcon,nrhcon,imat,rho,
     &              temp,ntmat_,ithermal)
            endif
!     
            vcontu(1,node)=xkin*dvi
            vcontu(2,node)=xtu*rho
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
!     density
!     
            if(compressible.eq.1) then
c     r=shcon(3,1,imat)
c     rho=vold(4,node)/
c     &        (r*(vold(0,node)-physcon(1)))
            else
               call materialdata_rho(rhcon,nrhcon,imat,rho,
     &              temp,ntmat_,ithermal)
            endif
!     
            vcontu(1,node)=0.d0
            vcontu(2,node)=800.d0*dvi/(xsolidsurf(j)**2)
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
      endif
!     
      return
      end
      
