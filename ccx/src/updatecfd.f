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
      subroutine updatecfd(vold,voldcon,v,nk,
     &  ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,iout,
     &  nmethod,convergence,physcon,iponoel,inoel,ithermal,
     &  nactdoh,iit,compressible,ismooth,voldtu,vtu,turbulent,
     &  inomat,nodeboun,ndirboun,nboun,mi,co,factor)
!
!     calculates 
!       vold (temperature,velocity and pressure)
!       voldcon (volumetric energy density, volumetric momentum
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
      real*8 v(0:mi(2),*),vold(0:mi(2),*),voldcon(0:4,*),
     &  rhcon(0:1,ntmat_,*),rho,c1,vmax(0:4),dummy,press,
     &  voldmax(0:4),cp,r,temp,temp0,c2,c3,tempnew,vel2,
     &  shcon(0:3,ntmat_,*),drho,dtemp,physcon(*),dpress,
     &  voldtu(2,*),vtu(2,*),co(3,*),factor
!
      if(ismooth.eq.0) then
!     
!     updates the volumetric energy density (only if ithermal>1),
!     the volumetric momentum density and the static pressure
!     
c         do j=0,4
c            vmax(j)=0.d0
c            voldmax(j)=0.d0
c         enddo
!     
!     volumetric energy density
!     
         if(ithermal.gt.1) then
            do i=1,nk
c               vmax(0)=vmax(0)+v(0,i)**2
c               voldmax(0)=voldmax(0)+voldcon(0,i)**2
               voldcon(0,i)=voldcon(0,i)+v(0,i)
            enddo
!
!           subtracting the boundary conditions
!
c            do i=1,nboun
c               if(ndirboun(i).eq.0) then
c                  vmax(0)=vmax(0)-v(0,nodeboun(i))**2
c               endif
c            enddo
!
         endif
!     
!     volumetric momentum density
!     
         do i=1,nk
            do j=1,3
c               vmax(j)=vmax(j)+v(j,i)**2
c               voldmax(j)=voldmax(j)+voldcon(j,i)**2
               voldcon(j,i)=voldcon(j,i)+v(j,i)
            enddo
         enddo
!     
!     volumetric turbulent density
!     
         if(turbulent.ne.0) then
            do i=1,nk
               voldtu(1,i)=voldtu(1,i)+vtu(1,i)
               voldtu(2,i)=voldtu(2,i)+vtu(2,i)
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
c                  vmax(4)=vmax(4)+v(4,i)**2
c                  voldmax(4)=voldmax(4)+voldcon(4,i)**2
                  voldcon(4,i)=voldcon(4,i)+v(4,i)
                  cycle
               endif
!     
            else
!     
!              thermal liquid: pressure was calculated
!     
c               vmax(4)=vmax(4)+v(4,i)**2
c               voldmax(4)=voldmax(4)+vold(4,i)**2
               vold(4,i)=vold(4,i)+v(4,i)
               c1=voldcon(0,i)
               c2=(voldcon(1,i)**2+voldcon(2,i)**2+voldcon(3,i)**2)/2.d0
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
                     vold(k,i)=voldcon(k,i)/rho
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
!     
c            vmax(4)=vmax(4)+v(4,i)**2
c            voldmax(4)=voldmax(4)+vold(4,i)**2
c            if((i.eq.321).or.(i.eq.322)) then
c               write(*,*) 'updatecfd ',i,vold(4,i),v(4,i)
c            endif
            vold(4,i)=vold(4,i)+v(4,i)
c            if((i.eq.321).or.(i.eq.322)) then
c               write(*,*) 'updatecfd ',i,vold(4,i),v(4,i)
c            endif
            voldcon(4,i)=rho
!     
!     storing the density
!     calculating the velocity
!     
            do k=1,3
               vold(k,i)=voldcon(k,i)/rho
            enddo
         enddo
      endif
!     
!     for steady state calculations: check convergence
!     
      if(ismooth.eq.0) then
c         convergence=0
c         do i=0,4
c            vmax(i)=dsqrt(vmax(i))
c            voldmax(i)=dsqrt(voldmax(i))
c         enddo
c         if(nmethod.eq.1) then
c            if(((dabs(vmax(0)).lt.1.d-8*dabs(voldmax(0))).or.
c     &           (dabs(voldmax(0)).lt.1.d-10)).and.
c     &           ((dabs(vmax(1)).lt.1.d-8*dabs(voldmax(1))).or.
c     &           (dabs(voldmax(1)).lt.1.d-10)).and.
c     &           ((dabs(vmax(2)).lt.1.d-8*dabs(voldmax(2))).or.
c     &           (dabs(voldmax(2)).lt.1.d-10)).and.
c     &           ((dabs(vmax(3)).lt.1.d-8*dabs(voldmax(3))).or.
c     &           (dabs(voldmax(3)).lt.1.d-10)).and.
c     &           ((dabs(vmax(4)).lt.1.d-8*dabs(voldmax(4))).or.
c     &           (dabs(voldmax(4)).lt.1.d-10)).and.
c     &           (iit.gt.1)) convergence=1
c         endif
c         write(*,'(i10,10(1x,e11.4))') iit,vmax(0),voldmax(0),
c     &       vmax(1),voldmax(1),vmax(2),voldmax(2),
c     &       vmax(3),voldmax(3),vmax(4),voldmax(4)
c         write(*,*) 'convergence ',convergence

c         factor=min(1.d0,1.01d0*factor)
c         if(dabs(voldmax(0)).gt.1.d-3) then
c            factor=min(factor,voldmax(0)/vmax(0)*0.001)
c         endif
c         if(dabs(voldmax(1)).gt.1.d-3) then
c            factor=min(factor,voldmax(1)/vmax(1)*0.001)
c         endif
c         if(dabs(voldmax(2)).gt.1.d-3) then
c            factor=min(factor,voldmax(2)/vmax(2)*0.001)
c         endif
c         if(dabs(voldmax(3)).gt.1.d-3) then
c            factor=min(factor,voldmax(3)/vmax(3)*0.001)
c         endif
c         if(dabs(voldmax(4)).gt.1.d-3) then
c            factor=min(factor,voldmax(4)/vmax(4)*0.001)
c         endif
      endif
!     
      return
      end
      
