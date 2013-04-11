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
!     calculate the maximal admissible pressure ratio pt2/pt1
!
!     1) assuming M2=1 for adiabatic respectively M2=1/dsqrt(kappa) for isotherm pipe choking 
!      M1 is calculated iteratively using a dichotomy scheme
!
!     2)the ratio of the critical pressure ratio  Qred_1/Qred_2crit=Pt2/Pt1
!     =D(M1)/D(M2_crit)is computed [D(M)=M*(1+0.5*(kappa-1)*M)**(-0.5*(kappa+1)/(kappa-1))]
!   
      subroutine pt2zpt1_crit(pt2,pt1,Tt1,lambda,kappa,r,l,d,A,iflag,
     &     xflow,inv,nactdog,voldgas,node1,node2,pt2zpt1_c,qred_crit,
     &     crit,case)
!     
      implicit none
!
      logical crit
!
      integer iflag,inv,nactdog(0:3,*),node2,node1,case
!     
      real*8 pt2,pt1,lambda,kappa,l,d,M1,pt2zpt1,pt2zpt1_c,
     &     km1,kp1,km1zk,kp1zk,e,xflow,Tt1,r,A,voldgas(0:3,*),
     &     xflow_crit,qred_crit,f1,f2,f3,m1_ac,m1_min,m1_max
!
!      write(*,*)  ''      
!      write(*,*) 'Pt2zPt1_crit.f'
!
!
!     useful variables and constants
!     
          km1=kappa-1.d0
          kp1=kappa+1.d0
          km1zk=km1/kappa
          kp1zk=kp1/kappa
          e=2.7182818d0
!     write(*,*) 'Pt1',Pt1
!     write(*,*) 'Pt2',Pt2
!
!     adiabatic case
!
       if(case.eq.0) then
!          write(*,*) 'adiabatic'
!     
!     computing M1 using dichotomy method
!     
          m1_max=1
          m1_min=0.001d0
          do
!     write(*,*) 'm1_min',m1_min
!     write(*,*) 'm1_max',m1_max
             m1_ac=(m1_min+m1_max)*0.5d0
!     write(*,*) 'm1_ac',m1_ac
!     
             f1=(1.d0-M1_min**2)/(kappa*M1_min**2)
     &            +0.5d0*kp1zk*log((1+0.5d0*km1)*M1_min**2
     &            /(1+0.5d0*km1*M1_min**2))/log(e)-lambda*l/d
!     write(*,*) 'f1',f1
!     
             f2=(1.d0-M1_ac**2)/(kappa*M1_ac**2)
     &            +0.5d0*kp1zk*log((1+0.5d0*km1)*M1_ac**2
     &            /(1+0.5d0*km1*M1_ac**2))/log(e)-lambda*l/d
!     write(*,*) 'f2',f2
!     
             f3=(1.d0-M1_max**2)/(kappa*M1_max**2)
     &            +0.5d0*kp1zk*log((1+0.5d0*km1)*M1_max**2
     &            /(1+0.5d0*km1*M1_max**2))/log(e)-lambda*l/d
!     write(*,*) 'f3',f3
!     
             if(abs(f3).le.1E-4) then
                M1=m1_ac
!     write(*,*) 'M1',M1
                exit
             endif
!     
             if((f3.gt.f2).and.(f2.ge.f1)) then
                if((f1.lt.0d0).and.(f2.lt.0d0)) then
                   m1_min=m1_ac
                else
                   m1_max=m1_ac
                endif
             elseif((f3.lt.f2).and.(f2.le.f1)) then  
                if((f3.lt.0d0).and.(f2.lt.0d0) )then
                   m1_max=m1_ac
                else
                   m1_min=m1_ac
                endif
             endif
          enddo
!
!     2) computing the critical pressure ratio in the adiabatic case
!     
          Pt2zPt1_c=M1*(2/kp1*(1+km1/2.d0*M1**2))**(-0.5d0*kp1/km1)
!     write(*,*) 'Pt2zPt1_crit=',Pt2zPt1_c
!     
!     isotherm case
!
       elseif (case.eq.1) then
!          write(*,*) 'isotherm'
!     write(*,*) 'Pt1',Pt1
!     write(*,*) 'Pt2',Pt2
!     
!     computing M1 using dichotomy method for choked conditions M2=1/dsqrt(kappa)
!     (1.d0-kappa*M1_min**2)/(kappa*M1_min**2)+log(kappa*M1_min**2)-lambda*l/d=0
!     
          m1_max=1/dsqrt(kappa)
          m1_min=0.1d0
          do
!     write(*,*) 'm1_min',m1_min
!     write(*,*) 'm1_max',m1_max
             m1_ac=(m1_min+m1_max)*0.5d0
!     write(*,*) 'm1_ac',m1_ac
!     
             f1=(1.d0-kappa*M1_min**2)/(kappa*M1_min**2)
     &            +log(kappa*M1_min**2)-lambda*l/d
!     write(*,*) 'f1',f1
!     
             f2=(1.d0-kappa*M1_ac**2)/(kappa*M1_ac**2)
     &            +log(kappa*M1_ac**2)-lambda*l/d
!     write(*,*) 'f2',f2
!     
             f3=(1.d0-kappa*M1_max**2)/(kappa*M1_max**2)
     &            +log(kappa*M1_max**2)-lambda*l/d
!     write(*,*) 'f3',f3
!     
             if(abs(f3).le.1E-5) then
                M1=m1_ac
!                write(*,*) 'M1_crit',M1
!                write(*,*) 'M2_crit',1/dsqrt(kappa)
                exit
             endif
!     
             if((f3.gt.f2).and.(f2.ge.f1)) then
                if((f1.lt.0d0).and.(f2.lt.0d0)) then
                   m1_min=m1_ac
                else
                   m1_max=m1_ac
                endif
             elseif((f3.lt.f2).and.(f2.le.f1)) then  
                if((f3.lt.0d0).and.(f2.lt.0d0) )then
                   m1_max=m1_ac
                else
                   m1_min=m1_ac
                endif
             endif
          enddo
!     
!     2) computing the critical pressure ratio in the isotherm case
!     
          Pt2zPt1_c=M1*dsqrt(kappa)*((1+0.5d0*km1/kappa)
     &         /(1+0.5d0*km1*M1**2))**(kappa/km1)
!          write(*,*) 'Pt2zPt1_crit=',Pt2zPt1_c
!     
       endif
!     
!     3)
!     
       pt2zpt1=pt2/pt1
!     write(*,*) 'Pt2/Pt1',Pt2zPt1
       if(Pt2zPt1.le.Pt2zPt1_c) then
!     write(*,*) 'Pipe is critical'
          crit=.true.
       endif
!     
       if (iflag.eq.1) then
          xflow_crit=inv*M1*Pt1*A/dsqrt(Tt1)*dsqrt(kappa/r)
     &         *(1+0.5d0*km1*M1**2)**(-0.5d0*kp1/km1) 
c          write(*,*) 'xflow_crit',xflow_crit
!         if(dabs(xflow).ge.dabs(xflow_crit))then
!             crit=.true.
!          endif
       endif
!     
          Qred_crit=M1*dsqrt(kappa/r)
     &         *(1+0.5d0*km1*M1**2)**(-0.5d0*kp1/km1)
!     
!          write(*,*) 'Qred_crit',qred_crit
!     
!       endif
!     
!     write(*,*) 'end of Pt2zPt1_crit.f'
!     write(*,*) ''
!     
      return
      end      
      
      
