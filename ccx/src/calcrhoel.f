!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2014 Guido Dhondt
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
      subroutine calcrhoel(ne,lakon,vel,rhcon,nrhcon,ielmat,ntmat_,
     &  ithermal,mi)
!
!     calculation of rho in the element centers (incompressible
!     fluids)
!
      implicit none
!
      character*8 lakon(*)
!
      integer ne,i,nrhcon(*),imat,ithermal,ntmat_,mi(*),ielmat(mi(3),*)
!
      real*8 t1l,vel(0:5,*),rho,rhcon(0:1,ntmat_,*)
!     
      do i=1,ne
         if(lakon(i)(1:1).ne.'F') cycle
         t1l=vel(0,i)
         imat=ielmat(1,i)
         call materialdata_rho(rhcon,nrhcon,imat,rho,t1l,ntmat_,
     &            ithermal)
         vel(5,i)=rho
      enddo
!            
      return
      end
