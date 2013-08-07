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
      subroutine ucreep(amat,iel,iint,t1l,epini,ep,dtime,svm,dsvm)
!
!     INPUT:
!
!     amat:  material name
!     iel:   element number
!     iint:  integration point number
!     t1l:   temperature
!     epini: equivalent creep strain at the start 
!            of the increment
!     ep:    present equivalent creep strain; values of ep < epini
!            are equivalent to ep=epini.
!     dtime: time increment
!
!     OUTPUT:
!
!     svm:   present Von Mises stress
!     dsvm:  derivative of the Von Mises true stress with respect
!            to the present equivalent creep strain.
!            Numerically: change the present equivalent
!            strain with a small amount, calculate the amount
!            of change this causes in the present Von Mises
!            true stress, and divide the latter amount through the
!            former amount.
!
      implicit none
!
      character*80 amat
      real*8 t1l,epini,ep,dtime,svm,dsvm
!
      integer iel,iint
      if(ep.le.epini) then
         svm=0.d0
         dsvm=1.d10
      else
         svm=((ep-epini)/(dtime*1.d-10))**0.2d0
         dsvm=((ep-epini)/(dtime*1.d-10))**(-0.8d0)/(5.d-10*dtime)
      endif
!
      RETURN
      end








