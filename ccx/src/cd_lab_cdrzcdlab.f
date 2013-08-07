!
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2013 Guido Dhondt
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
!     This subroutine is a dummy subroutine , to the knowledge of the authors
!     no public source is available for such  data
!

      subroutine  cd_lab_cdrzcdlab (t,s,hst,x,p1,p2,cd_cdrzcdlab)
!
      implicit none
!     
      real*8 t,s,hst,x,p1,p2,cd_cdrzcdlab
!
      t=t
      s=s
      hst=hst
      x=x 
      p1=p1
      p2=p2
      cd_cdrzcdlab=1.d0
!
      write(*,*) '*WARNING while using subroutine cd_lab_cdrzcdlab.f'
      write(*,*) 'cd implicitely taken equal to 1'
!
      return
!
      end
