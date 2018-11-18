!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2018 Guido Dhondt
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
      subroutine calcrhofacomp_ud(nface,vfa,shcon,ielmat,ntmat_,
     &  mi,ielfa,ipnei,vel,nef,flux,gradpel,gradtel,xxj,betam,
     &  xlet)
!
!     calculation of the density at the face centers
!     (compressible fluids)
!
      implicit none
!
      integer nface,i,j,imat,ntmat_,mi(*),ipnei(*),nef,iel1,iel2,
     &  ielmat(mi(3),*),ielfa(4,*),indexf
!
      real*8 t1l,vfa(0:7,*),shcon(0:3,ntmat_,*),vel(nef,0:7),flux(*),
     &  r,gradpel(3,*),gradtel(3,*),xxj(3,*),gamma,betam,phic,vud,
     &  vcd,xlet(*)
!     
c$omp parallel default(none)
c$omp& shared(nface,vfa,ielmat,ielfa,shcon,ipnei,flux,gradpel,gradtel,
c$omp&        xxj,betam,xlet,vel)
c$omp& private(i,j,t1l,imat,iel1,iel2,indexf,vcd,vud,r,gamma,phic)
c$omp do
      do i=1,nface
         t1l=vfa(0,i)
!
!        take the material of the first adjacent element
!
         imat=ielmat(1,ielfa(1,i))
         r=shcon(3,1,imat)
!
!        specific gas constant
!
         vfa(5,i)=vfa(4,i)/(r*t1l)
!
!        calculate gamma
!
         iel2=ielfa(2,i)
!
!        faces with only one neighbor need not be treated
!
         if(iel2.le.0) cycle
         iel1=ielfa(1,i)
         j=ielfa(4,i)
         indexf=ipnei(iel1)+j
!
         if(flux(indexf).ge.0.d0) then
!
c            vcd=vel(iel1,5)-vel(iel2,5)
c            if(dabs(vcd).lt.1.d-3*dabs(vel(iel1,5))) vcd=0.d0
c!
c            vud=2.d0*xlet(indexf)*
c     &           ((gradpel(1,iel1)-vel(iel1,5)*r*gradtel(1,iel1))
c     &             *xxj(1,indexf)+
c     &            (gradpel(2,iel1)-vel(iel1,5)*r*gradtel(2,iel1))
c     &             *xxj(2,indexf)+
c     &            (gradpel(3,iel1)-vel(iel1,5)*r*gradtel(3,iel1))
c     &             *xxj(3,indexf))/(r*vel(iel1,0))
c!
c            if(dabs(vud).lt.1.d-20) then
c!
c!           upwind difference
c!
c               vfa(5,i)=vel(iel1,5)
c               cycle
c            endif
c!     
c            phic=1.d0+vcd/vud
!
!              upwind difference
!
               vfa(5,i)=vel(iel1,5)
         else
!
c            vcd=vel(iel2,5)-vel(iel1,5)
c            if(dabs(vcd).lt.1.d-3*dabs(vel(iel2,5))) vcd=0.d0
c!
c            vud=-2.d0*xlet(indexf)*
c     &           ((gradpel(1,iel2)-vel(iel2,5)*r*gradtel(1,iel2))
c     &             *xxj(1,indexf)+
c     &            (gradpel(2,iel2)-vel(iel2,5)*r*gradtel(2,iel2))
c     &             *xxj(2,indexf)+
c     &            (gradpel(3,iel2)-vel(iel2,5)*r*gradtel(3,iel2))
c     &             *xxj(3,indexf))/(r*vel(iel2,0))
c!
c            if(dabs(vud).lt.1.d-20) then
c!
c!           upwind difference
c!
c               vfa(5,i)=vel(iel2,5)
c               cycle
c            endif
c!     
c            phic=1.d0+vcd/vud
!
!              upwind difference
!
            vfa(5,i)=vel(iel2,5)
         endif
!
      enddo
c$omp end do
c$omp end parallel
!            
      return
      end
