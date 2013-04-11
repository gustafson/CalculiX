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
      subroutine storeresidual(nactdof,b,fn,filab,ithermal,nk)
!
!     This routine is called in case of divergence:
!     stores the residual forces in fn and changes the
!     file storage labels so that the independent
!     variables (displacements and/or temperatures) and
!     the corresponding residual forces are stored in the
!     frd file
!
      implicit none
!
      character*6 filab(*)
!
      integer nactdof(0:3,*),ithermal,i,j,nk
!
      real*8 b(*),fn(0:3,*)
!
!     storing the residual forces in field fn
!
      do i=1,nk
         do j=0,3
            if(nactdof(j,i).gt.0) then
               fn(j,i)=b(nactdof(j,i))
            else
               fn(j,i)=0.d0
            endif
         enddo
      enddo
!
!     adapting the storage labels
!
      if(ithermal.ne.2) then
         filab(1)='U     '
      else
         filab(1)='      '
      endif
      if(ithermal.gt.1) then
         filab(2)='NT    '
      else
         filab(2)='      '
      endif
      do i=3,10
         filab(i)='      '
      enddo
      if(ithermal.ne.2) then
         filab(13)='RFRES '
      else
         filab(13)='      '
      endif
      if(ithermal.gt.1) then
         filab(14)='RFLRES'
      else
         filab(14)='      '
      endif      
!
      return
      end


