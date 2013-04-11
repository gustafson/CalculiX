!
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2005 Guido Dhondt
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
      subroutine limit_case_calc(a1,a2,pt1,Tt1,Tt2,xflow,zeta,r,kappa,
     &     pt2_lim,M2)
!     
!     For restrictor elements A1<A2
!     if A1 is critical, pt2 can be as low as possible without any effect on the flow
!     it is necessary to compute pt2_lim which satisfies the element equation
!
!     xflow*dsqrt(R*Tt1)/(A1*Pt1*dsqrt(kappa)-dsqrt(2/(kappa-1)*(Pt1/Pt2)
!     **((kappa-1)-1))/(zeta*kappa)))/(Pt1/Pt2)**((kappa+1)/(2*zeta*kappa))=0
!
!     Once Pt2_lim is calculated it is possible in turn to calculate the corresponding
!     Mach number M2 satisfying
!     xflow*dsqrt(R*Tt2)/(A2*Pt2*dsqrt(kappa)-M2/(1+(kappa-1)/2*M2**2)
!     **(0.5*(kappa+1)/(kappa-1))
!
      implicit none
!
!      integer
!
      real*8 pt1,Tt1,Tt2,xflow,zeta,r,kappa,pt2_lim,M2,Qred,expon1,
     &     expon2,root,a1,a2,f,df,pt2,km1,kp1,pt1pt2
!
      Pt2=0.99*Pt1
      pt1pt2=pt1/pt2
      km1=kappa-1
      kp1=kappa+1
      expon1=-0.5*(kp1)/(zeta*kappa)
      expon2=(km1)/(zeta*kappa)
      root=2/(km1)*((pt1pt2)**expon2-1.d0)
!      Qred=dsqrt(kappa/R)*(pt1pt2)**(-0.5d0*kp1/(kappa*zeta))
!     &           *dsqrt(2.d0/km1*((pt1pt2)**(km1/(kappa*zeta))-1d0))
!xflow*dsqrt(Tt1)/(A1*Pt1)
!
!     pt2_lim calculation
!
      do
         Qred=dsqrt(kappa/R)*(pt1pt2)**(-0.5d0*kp1/(kappa*zeta))
     &           *dsqrt(2.d0/km1*((pt1pt2)**(km1/(kappa*zeta))-1d0))
         root=2/(km1)*((pt1pt2)**expon2-1.d0)
         f=Qred-dsqrt(kappa/R)*dsqrt(root)*(Pt1Pt2)**expon1
!
         df=dsqrt(kappa/R)*1/pt2*(Pt1Pt2)**expon1*dsqrt(root)*
     &        (1/(zeta*kappa)*(Pt1Pt2)**expon2*root**-1-expon1)
!
         if(dabs(-f/df)/pt2.le.1E-6) then
            pt2_lim=pt2-f/df
            write(*,*) 'Pt2_lim',pt2_lim
            exit
         endif
!
         pt2=pt2-f/df
         pt1pt2=pt1/pt2
      enddo
!
!     M2_lim calculation
!
      M2=0.5
      Qred=xflow*dsqrt(R*Tt2)/(A2*Pt2_lim*dsqrt(kappa))
      expon1=-0.5*(kp1)/(km1)
      
      do 
         
         root=(1+0.5d0*(km1)*M2**2)
         f=Qred-M2*root**(expon1)
!
         df=root**expon1*(-1d0+0.5*(kp1)*M2**2*root**-1)
!
         if(dabs(-f/df)/M2.le.1E-6) then
            M2=M2-f/df
            write(*,*) 'M2',M2
            exit
         endif
!
         M2=M2-f/df
      enddo
!
      return
!
      end
         
