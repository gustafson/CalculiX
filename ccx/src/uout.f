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
      subroutine uout(v)
!
!     This routine allows the user to write user-defined output
!     to file. The output can be brought into the routine by commons
!     (FORTRAN77) or modules (FORTRAN90). The file management must
!     be taken care of by the user.
!
!     INPUT:
! 
!     v                  solution vector
!
!     OUTPUT: none
!           
      implicit none
!
      real*8 v(0:4,*)
!
      return
      end

