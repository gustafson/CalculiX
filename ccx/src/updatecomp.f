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
      subroutine updatecomp(vold,vcon,v,nk,
     &  ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,iout,
     &  nmethod,convergence,physcon,iponoel,inoel,ithermal,
     &  nactdoh,iit,compressible,ismooth,vcontu,vtu,turbulent,
     &  inomat,nodeboun,ndirboun,nboun,mi,co,factor)
!
!     calculates 
!       vold (temperature,velocity and pressure)
!       at the nodes from the conservative variables  
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
!     calculate the static temperature and the density
!     
      do i=1,nk
         if(inomat(i).eq.0) cycle
         imat=inomat(i)
         temp=vold(0,i)
!     
!              gas: density was calculated
!     
         rho=vcon(4,i)
         c1=(vcon(0,i)-(vcon(1,i)**2+vcon(2,i)**2+
     &        vcon(3,i)**2)/(2.d0*rho))/rho
!     
!              temperature has to be calculated
!     
         temp0=temp
         j=0
         do
            call materialdata_cp_sec(imat,ntmat_,temp,shcon,
     &           nshcon,cp,physcon)
            r=shcon(3,1,imat)
            temp=max(c1/(cp-r),1.d-2)+physcon(1)
cstart shallow
c                  temp=max(c1/(cp),1.d-2)+physcon(1)
cend shallow
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
cstart shallow
c               vold(4,i)=5.d0*(rho*rho-(0.005*co(1,i))**2)
cend shallow
!     
!     calculating the velocity
!     
         do k=1,3
            if(nactdoh(k,i).ne.0) then
               vold(k,i)=vcon(k,i)/rho
            endif
         enddo
      enddo
!     
      return
      end
      
