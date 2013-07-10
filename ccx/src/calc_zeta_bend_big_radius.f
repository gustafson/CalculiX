!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
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
!     Bend with a big radius (R/D>3,Idelchik) for ACC design
!     Written by Yavor Dobrev
!
      subroutine calc_zeta_bend_big_radius(zeta,R_D,reynolds,delta,deri)
!
!     Calculates the zeta-value of a bend with R/D>3 acc. to Idel'chik
!     This subroutine is also used in acctube.f
!
      implicit none
!
      integer iexp_acc(2),ier,deri,n1,n9,n10,n11,n13
!
      real*8
     &zeta,
     &R_D,
     &reynolds,
     &lambda_el,
     &delta
!
!     Bend for ACC-Design
!     Idelchik p.291
!
      REAL*8 lambda_el_idel(9,14)
      DATA lambda_el_idel /
     & 9.014,       3.1,  3.9,  4.4,   6.5,  12.5, 22.5, 40.d0, 50.d0,
     & 400.d0,      0.34, 0.3,  0.28,  0.26, 0.24, 0.22, 0.2,   0.18,
     & 600.d0,      0.26, 0.23, 0.22,  0.2,  0.18, 0.16, 0.15,  0.135,
     & 800.d0,      0.22, 0.19, 0.18,  0.16, 0.15, 0.14, 0.13,  0.105,
     & 1000.d0,     0.19, 0.17, 0.16,  0.14, 0.13, 0.12, 0.11,  0.09,
     & 2000.d0,     0.12, 0.11, 0.1,   0.09, 0.08, 0.075,0.07,  0.052,
     & 4000.d0,     0.078,0.07, 0.065, 0.06, 0.055,0.048,0.045, 0.04,
     & 6000.d0,     0.063,0.06, 0.056, 0.052,0.043,0.04, 0.038, 0.035,
     & 8000.d0,     0.058,0.055,0.052, 0.049,0.04, 0.037,0.035, 0.032,
     & 10000.d0,    0.055,0.053,0.049, 0.047,0.038,0.035,0.033, 0.03,
     & 20000.d0,    0.05, 0.047,0.045, 0.043,0.034,0.03, 0.028, 0.025,
     & 30000.d0,    0.048,0.045,0.043, 0.042,0.033,0.029,0.027, 0.023,
     & 50000.d0,    0.046,0.044,0.041, 0.040,0.03, 0.027,0.025, 0.022,
     & 100000.d0,   0.044,0.042,0.04,  0.038,0.028,0.026,0.023, 0.02/
!
      data n1 /1/
      data n9 /9/
      data n10 /10/
      data n11 /11/
      data n13 /13/
!
!     Calculation
!
!     Check if variables are in range
      if (R_D.lt.3.d0) then
         if (deri.eq.0)then
            write(*,*)'*WARNING in bend_big_radius:'
            write(*,*)'R/D outside valid range (R/D<3)'
            write(*,*)'R/D=',R_D
         endif
      endif
!
      if ((reynolds.lt.400.d0).or.(reynolds.gt.100000.d0)) then
         if (deri.eq.0)then
            write(*,*)'*WARNING in bend_big_radius:'
            write(*,*)'ACC-Bend:Reynolds outside valid range!'
            write(*,*)'400<Re<100000'
            write(*,*)'Re=',reynolds
            write(*,*)'Extrapolated values!'
         endif
      endif
!
!     Calculate lambda_el and zeta
!
      if (R_D.lt.50.d0) then
!        Extrapolation over R/D:constant
!        R/D<3, R/D>50 are handled elsewhere
!        Extrapolation over Re:linear for Re<400
!        constant for Re>100000
         DATA iexp_acc /0,10/
!
         call twodint(lambda_el_idel,n9,n11,R_D,
     &           reynolds,lambda_el,n1,IEXP_ACC,IER)
      else
!        Interpolation only about Re         
         call onedint(lambda_el_idel(1,2:14), 
     &      lambda_el_idel(9,2:14),n13,reynolds,lambda_el,n1,n1,n10,ier)       
      endif
!
      if ((reynolds.lt.400.d0).or.(reynolds.gt.100000.d0)) then
         if (deri.eq.0)then
            write(*,*)'lambda =',lambda_el
         endif
      endif
!
      zeta = 0.0175*lambda_el*R_D*delta          
!
      return
      end
