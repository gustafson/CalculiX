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
      subroutine compdt(nk,dt,nshcon,shcon,nrhcon,rhcon,vold,ntmat_,
     &  iponoel,inoel,dtimef,iexplicit,ielmat,physcon,dh)
!
!     - determine the time step for each node (stored in field dt
!       and the minimum value across all nodes (dtimef)
!
      implicit none
!
      logical iexplicit
!
      integer nk,i,iponoel(*),inoel(3,*),index,nelem,
     &  nshcon(*),nrhcon(*),ntmat_,ielmat(*),imat
!
      real*8 dtimef,dt(*),dvi,r,cp,rho,shcon(0:3,ntmat_,*),
     &  rhcon(0:1,ntmat_,*),vold(0:4,*),temp,vel,dtu,dtnu,physcon(3),
     &  dh(*)
!
!     
!     determining the time increment dt for each node.
!
!     edge nodes (fields iponoel and inoel are determined in precfd.f)
!
      do i=1,nk
         index=iponoel(i)
         if(index.le.0) cycle
!
!        look at an element belonging to the edge node
!
         nelem=inoel(1,index)
!     
!        determining the time increment
!
         imat=ielmat(nelem)
         temp=vold(0,i)
         call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,r,dvi,
     &        rhcon,nrhcon,rho)
!
!        density for gases
!         
         vel=dsqrt(vold(1,i)**2+vold(2,i)**2+vold(3,i)**2)
         if(vel.lt.1.d-10) vel=1.d-10
         if(iexplicit) then
            dt(i)=dh(i)/(dsqrt(cp*r*temp/(cp-r))+vel)
         else
            if(dabs(rho).lt.1.d-20) then
               rho=vold(4,i)/(r*(vold(0,i)-physcon(1)))
            endif
            dtu=dh(i)/vel
            dtnu=(dh(i)*dh(i)*rho)/(2.d0*dvi)
c            write(*,*) 'dtu ',vel,dtu,dtnu
            dt(i)=dtu*dtnu/(dtu+dtnu)
         endif
!
      enddo
!
!     middle nodes (interpolation between neighboring end nodes;
!     still to be done)
!      
!
!     determining the minimum height across the complete fluid mesh            
!
      dtimef=1.d30
      do i=1,nk
         if(dt(i).gt.0.d0) dtimef=min(dt(i),dtimef)
      enddo
!
!     the factor of 2 is introduced because dtimef is not 
!     sufficient for stability. Maybe a factor between 1 and 2
!     also works.
!
      dtimef=dtimef/2.d0
!
      return
      end
