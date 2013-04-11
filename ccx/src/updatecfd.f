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
     &  nactdoh)
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
      logical convergence
!
      integer nrhcon(*),ntmat_,nactdoh(0:4,*),
     &  nshcon(*),ielmat(*),nk,ithermal,i,j,k,index,iout,
     &  nmethod,imat,nelem,iponoel(*),inoel(3,*)
!
      real*8 v(0:4,*),vold(0:4,*),voldaux(0:4,*),rho0,rhov,
     &  rhcon(0:1,ntmat_,*),rho,c1,vmax(0:4),dummy,press,
     &  voldmax(0:4),cp,r,temp,dvi,temp0,c2,c3,tempnew,
     &  shcon(0:3,ntmat_,*),drho,dtemp,physcon(3),dpress
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
            vmax(0)=max(vmax(0),v(0,i))
            voldmax(0)=max(voldmax(0),voldaux(0,i))
            voldaux(0,i)=voldaux(0,i)+v(0,i)
         enddo
      endif
!
!     volumetric momentum density
!
      do i=1,nk
         do j=1,3
            vmax(j)=max(vmax(j),v(j,i))
            voldmax(j)=max(voldmax(j),voldaux(j,i))
            voldaux(j,i)=voldaux(j,i)+v(j,i)
         enddo
      enddo
!
!     calculate the static temperature and the density
!         
      if(ithermal.gt.1) then
!     
         do i=1,nk
            index=iponoel(i)
            if(index.le.0) cycle
            nelem=inoel(1,index)
            imat=ielmat(nelem)
            temp=vold(0,i)
            call materialdata_tg(imat,ntmat_,temp,
     &           shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rho)
            if(dabs(rho).lt.1.d-20) then
!
!              gas: density was calculated
!
c               vmax(4)=max(vmax(4),v(4,i))
c               voldmax(4)=max(voldmax(4),voldaux(4,i))
c               voldaux(4,i)=voldaux(4,i)+v(4,i)
c               rho=voldaux(4,i)
c               drho=v(4,i)
c               c1=v(0,i)+(voldaux(1,i)**2+voldaux(2,i)**2+
c     &                    voldaux(3,i)**2)*drho/(2.d0*rho**2)
c     &                  -(voldaux(1,i)*v(1,i)+voldaux(2,i)*v(2,i)+
c     &                    voldaux(3,i)*v(3,i))/rho
c               temp0=temp
c!     
c!     iterating to find the temperature
c!     
c               do
c                  call materialdata_tg_sec(imat,ntmat_,temp,
c     &               shcon,nshcon,cp,r,dvi,rhcon,nrhcon,dummy,physcon)
c                  dtemp=((r-cp)*(temp-physcon(1))*drho+c1)/(rho*(cp-r))
c                  tempnew=temp0+dtemp
c                  temp=tempnew
c                  exit
c               enddo
c               vold(0,i)=temp
cc               dpress=drho*r*(temp-physcon(1))+rho*r*dtemp
cc               vold(4,i)=vold(4,i)+dpress
c               vold(4,i)=rho*r*(temp-physcon(1))
c
               if(nactdoh(4,i).ne.0) then
!
!                 density was calculated
!                  
                  vmax(4)=max(vmax(4),v(4,i))
                  voldmax(4)=max(voldmax(4),voldaux(4,i))
                  voldaux(4,i)=voldaux(4,i)+v(4,i)
                  rho=voldaux(4,i)
                  c1=(voldaux(0,i)-(voldaux(1,i)**2+voldaux(2,i)**2+
     &                 voldaux(3,i)**2)/(2.d0*rho))/rho
!     
                  if(nactdoh(0,i).ne.0) then
!
!                    temperature has to be calculated
!
                     temp0=temp
                     do
                        call materialdata_tg_sec(imat,ntmat_,temp,
     &                       shcon,nshcon,cp,r,dvi,rhcon,nrhcon,dummy,
     &                       physcon)
                        temp=c1/(cp-r)+physcon(1)
                        if(i.eq.205) then
                        write(*,*) 'updatecfd...',i,rho,voldaux(0,i),
     &                      c1,temp0,temp
                        endif
                        if(dabs(temp-temp0).lt.1.d-4*temp) then
                           vold(0,i)=temp
                           exit
                        endif
                        temp0=temp
                     enddo
                  endif
!
!                 determining the pressure (gas equation)
!
                  vold(4,i)=rho*r*(temp-physcon(1))
                  if(i.eq.205) then
                  write(*,*) rho,r,temp,physcon(1),vold(4,i)
                  endif
               else
!
!                 the pressure is known; the density was not
!                 determined
!
                  if(nactdoh(0,i).ne.0) then
!
!                     the temperature has to be calculated
!
                     rhov=voldaux(1,i)**2+voldaux(2,i)**2+
     &                    voldaux(3,i)**2
                     if(dabs(rhov).lt.1.d-20) then
                        write(*,*) '*ERROR in updatecfd: node ',i
                        write(*,*) '       the pressure is known and'
                        write(*,*) '       the velocity is zero, yet'
                        write(*,*) '       the temperature is unknown'
                        write(*,*) '       and cannot be determined'
                        stop
                     endif
                     press=vold(4,i)
                     rhov=rhov*r
                     c1=(voldaux(0,i)+press)*press*2.d0/rhov
                     rhov=-2.d0*press*press/(rhov*r)
                     temp0=temp
                     do
                        call materialdata_tg_sec(imat,ntmat_,temp,
     &                       shcon,nshcon,cp,r,dvi,rhcon,nrhcon,dummy,
     &                       physcon)
                        temp=c1+rhov*cp+physcon(1)
                        if(dabs(temp-temp0).lt.1.d-4*temp) then
                           vold(0,i)=temp
                           exit
                        endif
                        temp0=temp
                     enddo
                  endif
                  voldaux(4,i)=press/(r*(temp-physcon(1)))
               endif
!                        
            else
!              
!              liquid: pressure was calculated
!
               vmax(4)=max(vmax(4),v(4,i))
               voldmax(4)=max(voldmax(4),vold(4,i))
               vold(4,i)=vold(4,i)+v(4,i)
               rho0=rho
               c1=v(0,i)
               c2=(voldaux(1,i)**2+voldaux(2,i)**2+voldaux(3,i)**2)/2.d0
               c3=voldaux(1,i)*v(1,i)+voldaux(2,i)*v(2,i)+
     &            voldaux(3,i)*v(3,i)
               temp0=temp
!     
!     iterating to find the temperature
!     
               do
                  call materialdata_tg_sec(imat,ntmat_,temp,
     &                 shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rho,physcon)
                  dtemp=(cp*temp*(rho-rho0)+c1+c2+(rho-rho0)/rho**2
     &                      -c3/rho)/(rho*(cp-r))
                  if(dabs(dtemp).lt.(1.d-4*temp)) exit
                  temp=temp0+dtemp
               enddo
               vold(0,i)=temp
            endif
!
!                 calculating the velocity
!                  
            do k=1,3
               if(nactdoh(k,i).ne.0) then
                  vold(k,i)=voldaux(k,i)/rho
               endif
            enddo
c               write(*,*) 'endofloop ',i,temp,vold(1,i),vold(2,i),
c     &             vold(3,i),vold(4,i)
         enddo
      else
!
         do i=1,nk
            index=iponoel(i)
            if(index.le.0) cycle
            nelem=inoel(1,index)
            imat=ielmat(nelem)
            temp=vold(0,i)
            call materialdata_tg(imat,ntmat_,temp,
     &           shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rho)
!     
!     gas
!     
            if(dabs(rho).lt.1.d-20) then
!
!              gas
!
!              gas: density was calculated
!
               vmax(4)=max(vmax(4),v(4,i))
               voldmax(4)=max(voldmax(4),voldaux(4,i))
               voldaux(4,i)=voldaux(4,i)+v(4,i)
               rho=voldaux(4,i)
               vold(4,i)=rho*r*temp
            else
!              
!              liquid: pressure was calculated
!
               vmax(4)=max(vmax(4),v(4,i))
               voldmax(4)=max(voldmax(4),vold(4,i))
               vold(4,i)=vold(4,i)+v(4,i)
               voldaux(4,i)=rho
            endif
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
      convergence=.false.
      if(nmethod.eq.1) then
         if((vmax(0).lt.1.d-5*voldmax(0)).and.
     &      (vmax(1).lt.1.d-5*voldmax(1)).and.
     &      (vmax(2).lt.1.d-5*voldmax(2)).and.
     &      (vmax(3).lt.1.d-5*voldmax(3)).and.
     &      (vmax(4).lt.1.d-5*voldmax(4))) convergence=.true.
      endif
!
!     no print requests
!
      if(iout.le.0) return
!
!     output in dat file (with *NODE PRINT or *EL PRINT)
!
      return
      end
