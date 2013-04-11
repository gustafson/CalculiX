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
      subroutine cfdconv(vold,voldcon,v,nk,
     &  ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,iout,
     &  nmethod,convergence,physcon,iponoel,inoel,ithermal,
     &  nactdoh,iit,compressible,ismooth,voldtu,vtu,turbulent,
     &  inomat,nodeboun,ndirboun,nboun,mi,co,factor,voldconini,
     &  dtimef)
!
!     calculates the change in solution
!
      implicit none
!
      integer convergence,compressible
!
      integer nrhcon(*),ntmat_,nactdoh(0:4,*),iit,turbulent,
     &  nshcon(*),ielmat(*),nk,ithermal,i,j,k,index,iout,
     &  nmethod,imat,nelem,iponoel(*),inoel(3,*),ismooth,
     &  inomat(*),node,nodeboun(*),ndirboun(*),nboun,mi(2)
!
      real*8 v(0:mi(2),*),vold(0:mi(2),*),voldcon(0:4,*),
     &  rhcon(0:1,ntmat_,*),rho,c1,vmax(0:4),dummy,press,
     &  vconmax(0:4),cp,r,temp,temp0,c2,c3,tempnew,vel2,
     &  shcon(0:3,ntmat_,*),drho,dtemp,physcon(*),dpress,
     &  voldtu(2,*),vtu(2,*),co(3,*),factor,voldconini(0:4,*),
     &  dtimef
!     
      do j=0,4
         vmax(j)=0.d0
         vconmax(j)=0.d0
      enddo
!
      if(compressible.eq.1) then
         do i=1,nk
            do j=0,4
               vmax(j)=vmax(j)+(voldcon(j,i)-voldconini(j,i))**2
               vconmax(j)=vconmax(j)+voldconini(j,i)**2
               voldconini(j,i)=voldcon(j,i)
            enddo
         enddo
      else
         do i=1,nk
            do j=0,3
               vmax(j)=vmax(j)+(voldcon(j,i)-voldconini(j,i))**2
               vconmax(j)=vconmax(j)+voldconini(j,i)**2
               voldconini(j,i)=voldcon(j,i)
            enddo
!
!           for incompressible fluids the pressure is stored
!           in vold(4,*), the initial pressure in 
!           voldconini(4,*)
!
            do j=4,4
               vmax(j)=vmax(j)+(vold(j,i)-voldconini(j,i))**2
               vconmax(j)=vconmax(j)+voldconini(j,i)**2
               voldconini(j,i)=vold(j,i)
            enddo
         enddo
      endif
!     
!     for steady state calculations: check convergence
!     
      convergence=0
      do i=0,4
         vmax(i)=dsqrt(vmax(i))
         vconmax(i)=dsqrt(vconmax(i))
      enddo
      if(nmethod.eq.1) then
         if(((dabs(vmax(0)).lt.1.d-8*dabs(vconmax(0))).or.
     &        (dabs(vconmax(0)).lt.1.d-10)).and.
     &        ((dabs(vmax(1)).lt.1.d-8*dabs(vconmax(1))).or.
     &        (dabs(vconmax(1)).lt.1.d-10)).and.
     &        ((dabs(vmax(2)).lt.1.d-8*dabs(vconmax(2))).or.
     &        (dabs(vconmax(2)).lt.1.d-10)).and.
     &        ((dabs(vmax(3)).lt.1.d-8*dabs(vconmax(3))).or.
     &        (dabs(vconmax(3)).lt.1.d-10)).and.
     &        ((dabs(vmax(4)).lt.1.d-8*dabs(vconmax(4))).or.
     &        (dabs(vconmax(4)).lt.1.d-10)).and.
     &        (iit.gt.1)) convergence=1
      endif
      write(*,'(i10,11(1x,e11.4))') iit,vmax(0),vconmax(0),
     &     vmax(1),vconmax(1),vmax(2),vconmax(2),
     &     vmax(3),vconmax(3),vmax(4),vconmax(4),dtimef
      factor=min(1.d0,1.01d0*factor)
      if(dabs(vconmax(0)).gt.1.d-3) then
         factor=min(factor,vconmax(0)/vmax(0)*0.001)
      endif
      if(dabs(vconmax(1)).gt.1.d-3) then
         factor=min(factor,vconmax(1)/vmax(1)*0.001)
      endif
      if(dabs(vconmax(2)).gt.1.d-3) then
         factor=min(factor,vconmax(2)/vmax(2)*0.001)
      endif
      if(dabs(vconmax(3)).gt.1.d-3) then
         factor=min(factor,vconmax(3)/vmax(3)*0.001)
      endif
      if(dabs(vconmax(4)).gt.1.d-3) then
         factor=min(factor,vconmax(4)/vmax(4)*0.001)
      endif
!     
      return
      end
      
