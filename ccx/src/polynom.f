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
      subroutine polynom(x,y,z,p)
!
!     calculates the polynomial terms for the Zienkiewicz-Zhu
!     stress recovery procedure
!
      implicit none
!
      real*8 p(20),x,y,z
!
      p(1)=1.d0
      p(2)=x
      p(3)=y
      p(4)=z
      p(5)=x*x
      p(6)=y*y
      p(7)=z*z
      p(8)=x*y
      p(9)=x*z
      p(10)=y*z
      p(11)=x*x*y
      p(12)=x*y*y
      p(13)=x*x*z
      p(14)=x*z*z
      p(15)=y*y*z
      p(16)=y*z*z
      p(17)=x*y*z
      p(18)=x*x*y*z
      p(19)=x*y*y*z
      p(20)=x*y*z*z
!
      return
      end
