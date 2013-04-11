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
      subroutine updatecfd(vold,voldaux,v,nk,
     &  ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,iout,
     &  nmethod,convergence,physcon,iponoel,inoel,ithermal,
     &  nactdoh,iit,compressible,ismooth,voldtu,vtu,turbulent,
     &  inomat,nodeboun,ndirboun,nboun)
!
!     calculates 
!       vold (temperature,velocity and pressure)
!       voldaux (volumetric energy density, volumetric momentum
!                density and density)
!       at the nodes    
!
!     prints if iout=1
!
      implicit none
!
      integer convergence,compressible
!
      integer nrhcon(*),ntmat_,nactdoh(0:4,*),iit,turbulent,
     &  nshcon(*),ielmat(*),nk,ithermal,i,j,k,index,iout,
     &  nmethod,imat,nelem,iponoel(*),inoel(3,*),ismooth,
     &  inomat(*),node,nodeboun(*),ndirboun(*),nboun
!
      real*8 v(0:4,*),vold(0:4,*),voldaux(0:4,*),
     &  rhcon(0:1,ntmat_,*),rho,c1,vmax(0:4),dummy,press,
     &  voldmax(0:4),cp,r,temp,temp0,c2,c3,tempnew,vel2,
     &  shcon(0:3,ntmat_,*),drho,dtemp,physcon(*),dpress,
     &  voldtu(2,*),vtu(2,*)
!
      if(ismooth.eq.0) then
!     
!     updates the volumetric energy density (only if ithermal>1),
!     the volumetric momentum density and the static pressure
!     
         do j=0,4
            vmax(j)=0.d0
            voldmax(j)=0.d0
         enddo
!     
!     volumetric energy density
!     
         if(ithermal.gt.1) then
            do i=1,nk
               vmax(0)=vmax(0)+v(0,i)**2
               voldmax(0)=voldmax(0)+voldaux(0,i)**2
               voldaux(0,i)=voldaux(0,i)+v(0,i)
            enddo
!
!           subtracting the boundary conditions
!
            do i=1,nboun
               if(ndirboun(i).eq.0) then
                  vmax(0)=vmax(0)-v(0,nodeboun(i))**2
               endif
            enddo
!
         endif
!     
!     volumetric momentum density
!     
         do i=1,nk
            do j=1,3
               vmax(j)=vmax(j)+v(j,i)**2
               voldmax(j)=voldmax(j)+voldaux(j,i)**2
               voldaux(j,i)=voldaux(j,i)+v(j,i)
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
                  vmax(4)=vmax(4)+v(4,i)**2
                  voldmax(4)=voldmax(4)+voldaux(4,i)**2
                  voldaux(4,i)=voldaux(4,i)+v(4,i)
                  cycle
               endif
               rho=voldaux(4,i)
               c1=(voldaux(0,i)-(voldaux(1,i)**2+voldaux(2,i)**2+
     &              voldaux(3,i)**2)/(2.d0*rho))/rho
!     
!              temperature has to be calculated
!     
               temp0=temp
               j=0
               do
                  call materialdata_cp_sec(imat,ntmat_,temp,shcon,
     &                 nshcon,cp,physcon)
                  r=shcon(3,1,imat)
c                  call materialdata_tg_sec(imat,ntmat_,temp,
c     &                 shcon,nshcon,cp,r,dvi,rhcon,nrhcon,dummy,
c     &                 physcon)
                  temp=max(c1/(cp-r),1.d-2)+physcon(1)
                  j=j+1
                  if(dabs(temp-temp0).lt.1.d-4*temp) then
                     vold(0,i)=temp
                     exit
                  endif
                  if(j.gt.100) then
                     stop
                  endif
                  temp0=temp
               enddo
!     
!              determining the pressure (gas equation)
!     
               vold(4,i)=rho*r*(temp-physcon(1))
!
!              determining the Mach number
!
               if(ismooth.eq.2) then
                  vel2=vold(1,i)**2+vold(2,i)**2+vold(3,i)**2
                  v(0,i)=cp/(cp-r)
                  v(1,i)=dsqrt((vold(1,i)**2+vold(2,i)**2+vold(3,i)**2)
     &                      /(v(0,i)*r*(temp-physcon(1))))
               endif
!     
            else
!     
!              thermal liquid: pressure was calculated
!     
               vmax(4)=vmax(4)+v(4,i)**2
               voldmax(4)=voldmax(4)+vold(4,i)**2
               vold(4,i)=vold(4,i)+v(4,i)
               c1=voldaux(0,i)
               c2=(voldaux(1,i)**2+voldaux(2,i)**2+voldaux(3,i)**2)/2.d0
               temp0=temp
               j=0
!     
!              iterating to find the temperature
!     
               do
                  call materialdata_cp_sec(imat,ntmat_,temp,shcon,
     &                 nshcon,cp,physcon)
                  call materialdata_rho(rhcon,nrhcon,imat,rho,
     &                 temp,ntmat_)
c                  temp=max((c1-c2/rho)/(rho*cp),1.d-2)+physcon(1)
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
            endif
!     
!     calculating the velocity
!     
            do k=1,3
               if(nactdoh(k,i).ne.0) then
                  vold(k,i)=voldaux(k,i)/rho
               endif
            enddo
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
     &           temp,ntmat_)
!     
            vmax(4)=vmax(4)+v(4,i)**2
            voldmax(4)=voldmax(4)+vold(4,i)**2
            vold(4,i)=vold(4,i)+v(4,i)
            voldaux(4,i)=rho
!     
!     storing the density
!     calculating the velocity
!     
            do k=1,3
               vold(k,i)=voldaux(k,i)/rho
            enddo
         enddo
      endif
!     
!     for steady state calculations: check convergence
!     
      if(ismooth.eq.0) then
         convergence=0
         do i=0,4
            vmax(i)=dsqrt(vmax(i))
            voldmax(i)=dsqrt(voldmax(i))
         enddo
         if(nmethod.eq.1) then
            if(((dabs(vmax(0)).lt.1.d-8*dabs(voldmax(0))).or.
     &           (dabs(voldmax(0)).lt.1.d-10)).and.
     &           ((dabs(vmax(1)).lt.1.d-8*dabs(voldmax(1))).or.
     &           (dabs(voldmax(1)).lt.1.d-10)).and.
     &           ((dabs(vmax(2)).lt.1.d-8*dabs(voldmax(2))).or.
     &           (dabs(voldmax(2)).lt.1.d-10)).and.
     &           ((dabs(vmax(3)).lt.1.d-8*dabs(voldmax(3))).or.
     &           (dabs(voldmax(3)).lt.1.d-10)).and.
     &           ((dabs(vmax(4)).lt.1.d-8*dabs(voldmax(4))).or.
     &           (dabs(voldmax(4)).lt.1.d-10)).and.
     &           (iit.gt.1)) convergence=1
         endif
         write(*,*) 'rho*totenergy ',vmax(0),voldmax(0),iit
         write(*,*) 'rho*vx ',vmax(1),voldmax(1)
         write(*,*) 'rho*vy ',vmax(2),voldmax(2)
         write(*,*) 'rho*vz ',vmax(3),voldmax(3)
         write(*,*) 'pressure(fluids)/density(gas) ',vmax(4),voldmax(4)
      endif
!     
      return
      end
      
