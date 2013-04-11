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
      subroutine uboun(boun,kstep,kinc,time,node,idof,coords,vold)
!
!     user subroutine uboun
!
!
!     INPUT:
!
!     kstep              step number
!     kinc               increment number
!     time(1)            current step time
!     time(2)            current total time
!     node               node number
!     idof               degree of freedom
!     coords  (1..3)     global coordinates of the node
!     vold(0..3,1..nk)   solution field in all nodes
!                        0: temperature
!                        1: displacement in global x-direction
!                           (or mass flow rate for fluid nodes)
!                        2: displacement in global y-direction
!                        3: displacement in global z-direction
!
!     OUTPUT:
!
!     boun               boundary value for degree of freedom idof
!                        in node "node"
!           
      implicit none
!
      integer kstep,kinc,node,idof 
      real*8 boun,time(2),coords(3),vold(0:4,*)
!
      boun=10.d0
!
      return
      end

