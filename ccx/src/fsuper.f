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
      real*8 function fsuper(time,t,a,b,h1,h2,h3,h4,h5,h6)
!
      implicit none
!
      real*8 time,t,a,b,h1,h2,h3,h4,h5,h6,h7,h8
!
      h7=dexp(h1*t)
      h8=dexp(-h2*t)
!
      fsuper=(a+b*time)*(h7*h3+h8*h4)-b*(h7*(t*h3-h5)-h8*(-t*h4-h6))
!
      return
      end
