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
!     this subroutine solves iteratively the following equation
!     to determine the pressure for which section A2 is critical
!
      subroutine pt2_lim_calc (pt2,pt1,a2,a1,kappa,r,zeta,pt2_lim,
     &     isothermal)
!
      implicit none
!
      logical isothermal
!
      integer i 
!
      real*8 pt2,pt1,a2,a1,kappa,r,pt2_lim,x,zeta,f,df,expon1,
     &     expon2,expon3,cte,a2a1,kp1,km1,delta_x,fact1,fact2,term
!
      x=0.999
!
!     x belongs to interval [0;1]
!
      if(zeta.ge.0d0) then
         kp1=kappa+1d0
         km1=kappa-1d0
         a2a1=a2/a1
         expon1=-0.5d0*kp1/(zeta*kappa)
         expon2=-0.5d0*kp1/km1
         cte=a2a1*(0.5*kp1)**expon2
         expon3=-km1/(zeta*kappa)
         i=0
!
         do
!     
            f=x**(-1d0)-cte*x**(expon1)
     &           *(2d0/km1*(x**expon3-1.d0))**-0.5d0
!     
            df=-1.d0/X**2-cte*(x**expon1
     &           *(2d0/km1*(x**expon3-1.d0))**-0.5d0)
     &           *(expon1/X-1d0/km1*expon3*x**(expon3-1d0)
     &           *(2d0/km1*(x**expon3-1.d0))**(-1.d0))
            
            delta_x=-f/df
!     
            if(( dabs(delta_x/x).le.1.E-8)
     &           .or.(dabs(delta_x/1d0).le.1.E-10)) then
!
               pt2_lim=pt1*X
!
               exit
            endif
!     
            x=delta_x+x
!     
         enddo
!
      else
!
         do 
            kp1=kappa+1d0
            km1=kappa-1d0
            a2a1=a2/a1
            expon1=kp1/(zeta*kappa)
            expon2=km1/(zeta*kappa)
            expon3=kp1/km1
            cte=a2a1**2*(0.5*kp1)**-expon3*(2/km1)**-1
            fact1=x**-expon1
            fact2=x**-expon2
            term=fact2-1
!     
            f=x**-2-cte*fact1*term**-1
!     
            df=-2*x**-3-cte*(x**(-expon1-1)*term**-1)
     &           *(-expon1+expon2*(X**-expon2)*fact2*term**-1)
!     
            delta_x=-f/df
!     
            if(( dabs(delta_x/x).le.1.E-8)
     &           .or.(dabs(delta_x/1d0).le.1.E-10)) then
               pt2_lim=pt1*X
               exit
            endif
!     
            x=delta_x+x
!     
         enddo
!     
      endif
      
      return
      end
