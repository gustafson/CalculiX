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
      subroutine dload(f,kstep,kinc,time,noel,npt,layer,kspt,
     &     coords,jltyp,loadtype)
!
!     user subroutine dload
!
!
!     INPUT:
!
!     kstep              step number
!     kinc               increment number
!     time(1)            current step time
!     time(2)            current total time
!     noel               element number
!     npt                integration point number
!     layer              currently not used
!     kspt               currently not used
!     coords(1..3)       global coordinates of the integration point
!     jltyp              loading face kode:
!                        21 = face 1 
!                        22 = face 2 
!                        23 = face 3 
!                        24 = face 4 
!                        25 = face 5 
!                        26 = face 6
!     loadtype           load type label
!
!     OUTPUT:
!
!     f                  magnitude of the distributed load
!           
      implicit none
!
      character*20 loadtype
      integer kstep,kinc,noel,npt,jltyp,layer,kspt
      real*8 f,time(2),coords(3)
!
c      f=100.d0*coords(2)
      f=(1.d0-coords(2))*(1.d0-coords(3))
c      f=1.d0
      write(*,100) jltyp-20,coords(1),coords(2),coords(3)
 100  format('     &',i5,3(',',d22.15))
!
      return
      end

