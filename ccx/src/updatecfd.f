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
      subroutine updatecfd(vold,vcon,v,nk,
     &  ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,iout,
     &  nmethod,convergence,physcon,iponoel,inoel,ithermal,
     &  nactdoh,iit,compressible,ismooth,vcontu,vtu,turbulent,
     &  inomat,nodeboun,ndirboun,nboun,mi,co,factor)
!
!     calculates 
!       vold (temperature,velocity and pressure: for incompressible
!             fluids)
!       vcon (volumetric energy density, volumetric momentum
!                density and density)
!       at the nodes    
!
!     prints if iout=1
!
      implicit none
!
      integer convergence,compressible
!
      integer nrhcon(*),ntmat_,nactdoh(0:4,*),iit,turbulent,mi(*),
     &  nshcon(*),ielmat(mi(3),*),nk,ithermal,i,j,k,index,iout,
     &  nmethod,imat,nelem,iponoel(*),inoel(3,*),ismooth,
     &  inomat(*),node,nodeboun(*),ndirboun(*),nboun
!
      real*8 v(0:mi(2),*),vold(0:mi(2),*),vcon(0:4,*),
     &  rhcon(0:1,ntmat_,*),rho,c1,vmax(0:4),dummy,press,
     &  voldmax(0:4),cp,r,temp,temp0,c2,c3,tempnew,vel2,
     &  shcon(0:3,ntmat_,*),drho,dtemp,physcon(*),dpress,
     &  vcontu(2,*),vtu(2,*),co(3,*),factor
!
      if(ismooth.eq.0) then
!     
!     updates the volumetric energy density (only if ithermal>1),
!     the volumetric momentum density and the static pressure
!     
!     volumetric energy density
!     
         if(ithermal.gt.1) then
            do i=1,nk
               vcon(0,i)=vcon(0,i)+v(0,i)
            enddo
         endif
!     
!     volumetric momentum density
!     
         do i=1,nk
            do j=1,3
               vcon(j,i)=vcon(j,i)+v(j,i)
            enddo
         enddo
!     
!     volumetric turbulent density
!     
         if(turbulent.ne.0) then
            do i=1,nk
               vcontu(1,i)=vcontu(1,i)+vtu(1,i)
               vcontu(2,i)=vcontu(2,i)+vtu(2,i)
c               if(vcontu(1,i).lt.0.d0) vcontu(1,i)=0.d0
c               if(vcontu(2,i).lt.1.d0) vcontu(2,i)=1.d0
            enddo
         endif
      endif
!     
!     calculate the static temperature and the density
!     
      if(ithermal.gt.1) then
!     
         do i=1,nk
            if((compressible.eq.0).or.(ismooth.gt.0)) then
               if(inomat(i).eq.0) cycle
               imat=inomat(i)
               temp=vold(0,i)
            endif
!     
            if(compressible.eq.1) then
!     
!              gas: density was calculated
!     
               if(ismooth.eq.0) then
                  vcon(4,i)=vcon(4,i)+v(4,i)
                  cycle
               endif
!     
            else
!     
!              thermal liquid: pressure was calculated
!     
               vold(4,i)=vold(4,i)+v(4,i)
               c1=vcon(0,i)
               c2=(vcon(1,i)**2+vcon(2,i)**2+vcon(3,i)**2)/2.d0
               temp0=temp
               j=0
!     
!              iterating to find the temperature
!     
               do
                  call materialdata_cp_sec(imat,ntmat_,temp,shcon,
     &                 nshcon,cp,physcon)
                  call materialdata_rho(rhcon,nrhcon,imat,rho,
     &                 temp,ntmat_,ithermal)
                  temp=(c1-c2/rho)/(rho*cp)+physcon(1)
                  j=j+1
                  if((dabs(temp-temp0).lt.1.d-4*dabs(temp)).or.
     &               (dabs(temp-temp0).lt.1.d-10)) then
                     vold(0,i)=temp
                     exit
                  endif
                  if(j.gt.100) then
                     write(*,*) 
     &                 '*ERROR in updatecfd: too many iterations'
                     stop
                  endif
                  temp0=temp
               enddo
!     
!     calculating the velocity
!     
               do k=1,3
                  if(nactdoh(k,i).ne.0) then
                     vold(k,i)=vcon(k,i)/rho
                  endif
               enddo
            endif
         enddo
      else
!     
!     athermal liquid calculation
!     
         do i=1,nk
            if(inomat(i).eq.0) cycle
            imat=inomat(i)
            temp=vold(0,i)
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
            vold(4,i)=vold(4,i)+v(4,i)
            vcon(4,i)=rho
!     
!     storing the density
!     calculating the velocity
!     
            do k=1,3
               vold(k,i)=vcon(k,i)/rho
            enddo
         enddo
      endif
!     
      return
      end
      
