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
      subroutine applybounv(nodeboun,ndirboun,nboun,xbounact,
     &  ithermal,nk,iponoel,inoel,vold,vcontu,t1act,isolidsurf,
     &  nsolidsurf,xsolidsurf,nfreestream,ifreestream,turbulent,
     &  vcon,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
     &  compressible,ismooth,nmpc,nodempc,ipompc,coefmpc,inomat,
     &  mi)
!
!     applies velocity boundary conditions
!
      implicit none
!
      integer turbulent,compressible,
     &  nrhcon(*),mi(*),ielmat(mi(3),*),ntmat_,nodeboun(*),
     &  isolidsurf(*),
     &  ndirboun(*),nshcon(*),nk,i,nboun,node,imat,ithermal,iponoel(*),
     &  inoel(3,*),nsolidsurf,ifreenode,ifreestream(*),nfreestream,k,
     &  index,ismooth,indexi,nodei,nmpc,nodempc(3,*),ipompc(*),
     &  ist,ndir,ndiri,inomat(*)
!
      real*8 rhcon(0:1,ntmat_,*),rho,vold(0:mi(2),*),xbounact(*),shcon,
     &  vcontu(2,*),t1act(*),temp,xsolidsurf(*),reflength,
     &  refkin,reftuf,refvel,vcon(0:4,*),physcon(*),v(0:mi(2),*),
     &  rhoi,coefmpc(*),residu,size,correction,sum
!
!     taking velocity MPC's into account
!
      do i=1,nmpc
         ist=ipompc(i)
         ndir=nodempc(2,ist)
         if((ndir.lt.1).or.(ndir.gt.3)) cycle
         node=nodempc(1,ist)
!     
!     check whether fluid MPC
!     
         imat=inomat(node)
         if(imat.eq.0) cycle
         if(compressible.eq.0) then
!     
!     determining rho from the material constants (for incompressible
!     fluids)
!     
            temp=vold(0,node)
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
         else
!     
!     determining rho from the solution field (for compressible fluids)
!     
            rho=vcon(4,node)+v(4,node)
         endif
!     
         index=nodempc(3,ist)
         sum=0.d0
         if(index.ne.0) then
            do
               nodei=nodempc(1,index)
               ndiri=nodempc(2,index)
!     
!     determining rho
!     
               if(compressible.eq.0) then
!     
!     determining rho from the material constants (for incompressible
!     fluids)
!     
                  imat=inomat(nodei)
                  if(imat.eq.0) cycle
                  temp=vold(0,nodei)
                  call materialdata_rho(rhcon,nrhcon,imat,rhoi,
     &                 temp,ntmat_,ithermal)
               else
!     
!     determining rho from the solution field (for compressible fluids)
!     
                  rhoi=vcon(4,nodei)+v(4,nodei)
               endif
!     
               sum=sum+coefmpc(index)*
     &              (vcon(ndiri,nodei)+v(ndiri,nodei))/rhoi
               index=nodempc(3,index)
               if(index.eq.0) exit
            enddo
         endif
         v(ndir,node)=-sum/coefmpc(ist)*rho-vold(ndir,node)
      enddo
!     
      return
      end
      
