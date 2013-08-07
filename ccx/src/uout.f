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
      subroutine uout(v,mi,ithermal)
!
!     This routine allows the user to write user-defined output
!     to file. The output can be brought into the routine by commons
!     (FORTRAN77) or modules (FORTRAN90). The file management must
!     be taken care of by the user.
!
!     INPUT:
! 
!     v                  solution vector
!     mi(1)              max # of integration points per element (max
!                        over all elements)
!     mi(2)              max degree of freedomm per node (max over all
!                        nodes) in fields like v(0:mi(2))...
!     ithermal(1)        applies to the present step
!                        0: no thermal effects are taken into account
!                        1: thermal boundary conditions for mechanical
!                           calculations
!                        2: heat transfer calculation
!                        3: coupled temperature-displacement calculation
!     ithermal(2)        applies to the complete input deck:
!                        0: no thermal effects are taken into account
!                        1: only mechanical steps
!                        2: only heat transfer steps
!                        3: at least one mechanical and one heat transfer
!                           step, or at least one coupled temperature-
!                           displacement step
!
!     OUTPUT: none
!           
      implicit none
!
      integer mi(*),ithermal(*)
!
      real*8 v(0:mi(2),*)
!
      return
      end

