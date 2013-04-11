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
     &  iponoel,inoel,dtimef,iexplicit,ielmat,physcon,dh,cocon,
     &  ncocon,ithermal)
!
!     - determine the time step for each node (stored in field dt
!       and the minimum value across all nodes (dtimef)
!
      implicit none
!
      logical iexplicit
!
      integer nk,i,iponoel(*),inoel(3,*),index,nelem,ithermal,
     &  nshcon(*),nrhcon(*),ntmat_,ielmat(*),imat,ncocon(2,*)
!
      real*8 dtimef,dt(*),dvi,r,cp,rho,shcon(0:3,ntmat_,*),
     &  rhcon(0:1,ntmat_,*),vold(0:4,*),temp,vel,dtu,dtnu,physcon(*),
     &  dh(*),cocon(0:6,ntmat_,*),dtal,cond
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
c         call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,r,dvi,
c     &        rhcon,nrhcon,rho)
!
!        density for gases
!         
         vel=dsqrt(vold(1,i)**2+vold(2,i)**2+vold(3,i)**2)
         if(iexplicit) then
            call materialdata_cp(imat,ntmat_,temp,shcon,nshcon,cp)
            r=shcon(3,1,imat)
            dt(i)=dh(i)/(dsqrt(cp*r*temp/(cp-r))+vel)
         else
            call materialdata_dvi(imat,ntmat_,temp,shcon,nshcon,dvi)
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_)
            if(vel.lt.1.d-10) vel=1.d-10
            dtu=dh(i)/vel
            if(dvi.lt.1.d-10) dvi=1.d-10
            dtnu=dh(i)*dh(i)*rho/(2.d0*dvi)
            dt(i)=dtu*dtnu/(dtu+dtnu)
            if(ithermal.gt.1) then
               call materialdata_cond(imat,ntmat_,temp,cocon,ncocon,
     &                cond)
               call materialdata_cp(imat,ntmat_,temp,shcon,nshcon,cp)
               if(cond.lt.1.d-10) cond=1.d-10
               dtal=dh(i)*dh(i)*rho*cp/(2.d0*cond)
               dt(i)=(dt(i)*dtal)/(dt(i)+dtal)
            endif
c            write(*,*) 'compdt ',i,dtu,dtnu,dt(i),dh(i),rho,dvi,vel
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
      return
      end
