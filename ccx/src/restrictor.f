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
      subroutine restrictor(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,iflag,voldgas,xflow,f,
     &     nodef,idirf,df,cp,r,physcon,dvi,numf)
!     
!     pressure loss element with partial total head loss 
!     
      implicit none
!     
      logical identity,crit,isothermal
      character*8 lakon(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(5),idirf(5),index,iflag,
     &     inv,ipkon(*),kon(*),kgas,case
!     
      real*8 prop(*),voldgas(0:3,*),xflow,f,df(5),kappa,R,d,
     &     Tt1,Tt2,pt1,pt2,cp,physcon(3),km1,dvi,
     &     kp1,kdkm1,reynolds,kdkp1,
     &     pt2pt1,pt1pt2,pt1pt2_crit,qred_crit,qred1,qred2,zeta,
     &     A1,A2,root, expon1,expon2,expon3,fact1,fact2,sqrt,pi,
     &     pt2_lim,Ts1,Ts2,M2,M1
!     
      numf=4
!     
      if (iflag.eq.0) then
         identity=.true.
!     
         if(nactdog(2,node1).ne.0)then
            identity=.false.
         elseif(nactdog(2,node2).ne.0)then
            identity=.false.
         elseif(nactdog(1,nodem).ne.0)then
            identity=.false.
         endif
!     
      elseif (iflag.eq.1)then
!     
         isothermal=.false.
         index=ielprop(nelem)
         kappa=(cp/(cp-R))
         kp1=kappa+1d0
         km1=kappa-1d0
!     
!     defining surfaces for branches elements
!     
         if(lakon(nelem)(2:6).eq.'REBRJ') then
            if(nelem.eq.int(prop(index+2))) then
               A1=prop(index+5)
               A2=A1
            elseif(nelem.eq.int(prop(index+3)))then
               A1=prop(index+6)
               A2=A1
            endif
         elseif(lakon(nelem)(2:6).eq.'REBRS') then
            if(nelem.eq.int(prop(index+2))) then
               A1=prop(index+5)
               A2=A1
            elseif(nelem.eq.int(prop(index+3)))then
               A1=prop(index+6)
               A2=A1
            endif
!     
!     for other Restrictor elements         
!     
         else
            A1=prop(index+1)
            A2=prop(index+2)
         endif
!     
         zeta=1.d0
!     
         pt1=voldgas(2,node1)
         pt2=voldgas(2,node2)
!     
         if(pt1.ge.pt2) then
            inv=1
            Tt1=voldgas(0,node1)+physcon(1)
            Tt2=voldgas(0,node2)+physcon(1)
         else
            inv=-1
            pt1=voldgas(2,node2)
            pt2=voldgas(2,node1)
            Tt1=voldgas(0,node2)+physcon(1)
            Tt2=voldgas(0,node1)+physcon(1)
         endif
!     
         pt1pt2=pt1/pt2
         pt2pt1=1/pt1pt2
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
!     
         if(.not.isothermal) then
            pt1pt2_crit=(0.5d0*kp1)**(zeta*kdkm1)
         else
            pt1pt2_crit=0.5d0*(3*kappa-1)**(zeta*kdkm1)
         endif
!     
         if(pt1pt2.gt.pt1pt2_crit) then
            crit=.true.
            pt1pt2=pt1pt2_crit
         endif
!     
         if(A1.le.A2) then
!     
            Qred1=dsqrt(kappa/R)*pt1pt2**(-0.5d0*kp1/(kappa*zeta))
     &           *dsqrt(2.d0/km1*(pt1pt2**(km1/(kappa*zeta))-1d0))
!     
            Qred2=pt1pt2*A1/A2*Qred1
!     
            if(.not.isothermal) then
               Qred_crit=dsqrt(kappa/R)*(1.d0+0.5d0*km1)
     &              **(-0.5d0*kp1/km1)
            else
               Qred_crit=dsqrt(1/R)*(1+0.5*km1/kappa)
     &              **(-0.5d0*kp1/km1)
            endif
!     
            if (Qred2.lt.Qred_crit) then
               if((Qred1.gt.Qred_crit).or.(pt1pt2.gt.pt1pt2_crit)) then
                  xflow=inv*A1*pt1*Qred_crit/dsqrt(Tt1)
               else
                  xflow=inv*A1*pt1*Qred1/dsqrt(Tt1)
               endif
            else
               call pt2_lim_calc(pt2,pt1,a2,a1,kappa,r,zeta,pt2_lim,
     &              isothermal)
!     
               xflow=inv*A2*pt2_lim*Qred_crit/dsqrt(Tt2)
!     
            endif
!     
         else
            Qred2=dsqrt(kappa/R)*pt1pt2**(-0.5d0*kp1/(kappa*zeta))
     &           *dsqrt(2.d0/km1*(pt1pt2**(km1/(kappa*zeta))-1d0))
            Qred1=pt2pt1*A2/A1*Qred2
            Qred_crit=(1.d0+0.5d0*km1)**(-0.5d0*kp1/km1)
!     
            if(Qred2.gt.Qred_crit) then
               xflow=inv*A2*pt2*Qred_crit/dsqrt(Tt2)
            else
               xflow=inv*A2*pt2*Qred2/dsqrt(Tt2)
            endif
         endif
!     
      elseif (iflag.eq.2)then
!     
         isothermal=.false.
         pi=4.d0*datan(1.d0)
         kappa=(cp/(cp-R))
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         kdkp1=kappa/kp1
!     
         pt1=voldgas(2,node1)
         pt2=voldgas(2,node2)
!     
         if(pt1.ge.pt2) then
            inv=1
            xflow=voldgas(1,nodem)
            Tt1=voldgas(0,node1)+physcon(1)
            Tt2=voldgas(0,node2)+physcon(1)
!     
            nodef(1)=node1
            nodef(2)=node1
            nodef(3)=nodem
            nodef(4)=node2
!     
         else
            inv=-1
            pt1=voldgas(2,node2)
            pt2=voldgas(2,node1)
            xflow=-voldgas(1,nodem)
            Tt1=voldgas(0,node2)+physcon(1)
            Tt2=voldgas(0,node1)+physcon(1)
!     
            nodef(1)=node2
            nodef(2)=node2
            nodef(3)=nodem
            nodef(4)=node1 
         endif
!     
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=2
!     
         index=ielprop(nelem)
!     
!     defining surfaces for branches elements
!     
         if(lakon(nelem)(2:6).eq.'REBRJ') then
            if(nelem.eq.int(prop(index+2))) then
               A1=prop(index+5)
               A2=A1
            elseif(nelem.eq.int(prop(index+3)))then
               A1=prop(index+6)
               A2=A1
            endif
         elseif(lakon(nelem)(2:6).eq.'REBRS') then
            if(nelem.eq.int(prop(index+2))) then
               A1=prop(index+5)
               A2=A1
            elseif(nelem.eq.int(prop(index+3)))then
               A1=prop(index+6)
               A2=A1
            endif
!     
!     for other Restrictor elements         
!     
         else
            A1=prop(index+1)
            A2=prop(index+2)
         endif
!     
!     calculation of the dynamic viscosity 
!     
         if( lakon(nelem)(2:3).eq.'RE') then
            case=0
         elseif(lakon(nelem)(2:5).eq.'REEX') then
            if(lakon(prop(index+4))(2:6).eq.'GAPIA') then
               case=0
            elseif(lakon(prop(index+4))(2:6).eq.'GAPII') then
               case=1
            endif
!     
         endif
!     
         call ts_calc(xflow,Tt1,Pt1,kappa,r,A1,Ts1,case)
!     
         call ts_calc(xflow,Tt2,Pt2,kappa,r,A2,Ts2,case)
!     
         if (A1.le.A2) then
            if(dabs(dvi).lt.1E-30) then
               kgas=0
               call dynamic_viscosity(kgas,Ts1,dvi)
            endif
         else
            if(dabs(dvi).lt.1E-30) then
               kgas=0
               call dynamic_viscosity(kgas,Ts2,dvi)
            endif
         endif
!     
!     Reynolds number calculation
!     
         if (lakon(nelem)(2:5).eq.'REBR') then
            d=dsqrt(4d0*A1/Pi)
            reynolds=dabs(xflow)*d/(dvi*A1)     
         else
            d=prop(index+3)
!     
            if(A1.le.A2) then
               reynolds=dabs(xflow)*d/(dvi*A1)
            else
               reynolds=dabs(xflow)*d/(dvi*A2)
            endif
         endif
!     
         call zeta_calc(nelem,prop,ielprop,lakon,reynolds,zeta,
     &        isothermal,kon,ipkon,R,Kappa,voldgas)
!     
         if(zeta.lt.0) then
            pt1=voldgas(2,node1)
            pt2=voldgas(2,node2)
            xflow=voldgas(1,nodem)
            Tt2=voldgas(0,node2)
            Tt1=voldgas(0,node1)
!     
            nodef(1)=node1
            nodef(2)=node1
            nodef(3)=nodem
            nodef(4)=node2
!     
         endif
!     
         if(.not.isothermal) then
            pt1pt2_crit=(0.5d0*kp1)**(zeta*kdkm1)
         else
            pt1pt2_crit=0.5d0*(3*kappa-1)**(zeta*kdkm1)
         endif
         pt1pt2=pt1/pt2
!     
!     Mach number caclulation
!     
         M1=dsqrt(2d0/km1*(Tt1/Ts1-1d0))
         if((1.d0-M1).le.1E-6) then
            if(zeta.gt.0d0) then
               call limit_case_calc(a1,a2,pt1,Tt1,Tt2,xflow,
     &              zeta,r,kappa,pt2_lim,M2)
!     
            endif
         else
            M2=dsqrt(2d0/km1*(Tt2/Ts2-1d0))
         endif
!     
!     Section A1 smaller than or equal to section A2
!     or for all BRANCHES ELEMENTS
!     
         if (A1.le.A2) then
!     
!     definition of the reduced mass flows
!     
            if(zeta.gt.0) then
!     
               Qred1=dsqrt(kappa/R)*pt1pt2**(-0.5d0*kp1/(kappa*zeta))
     &              *dsqrt(2.d0/km1*(pt1pt2**(km1/(kappa*zeta))-1d0))
!     
            elseif(zeta.lt.0d0) then
!     
               Qred1=dabs(xflow)*dsqrt(Tt1)/(pt1*A1)
!     
            endif
!     
            Qred2=pt1pt2*A1/A2*Qred1
!     
            if(.not.isothermal) then
               Qred_crit=dsqrt(kappa/R)*(1.d0+0.5d0*km1)
     &              **(-0.5d0*kp1/km1)
            else
               Qred_crit=dsqrt(1/R)*(1+0.5*km1/kappa)
     &              **(-0.5d0*kp1/km1)
            endif
!     
!     case zeta greater than zero
!     
            if(zeta.gt.0) then
!     
!     definition of the coefficients 
!     
               sqrt=dsqrt(R*Tt1/kappa)
               expon1=-0.5d0*kp1/(zeta*kappa)
               fact1=pt1pt2**expon1
               expon2=km1/(zeta*kappa)
               fact2=pt1pt2**expon2
               expon3=1d0/(zeta*kappa)
               root=2d0/km1*(fact2-1d0)
!     
               if(Qred2.lt.Qred_crit) then
!     
                  if((Qred1.lt.Qred_crit)
     &                 .and.(pt1pt2.lt.pt1pt2_crit))then
!     
!     
!     section 1 is not critical
!     
!     residual
!     
                     f=xflow*sqrt/(A1*Pt1)-fact1*dsqrt(root)
!     
!     
!     pressure node1
!     
                     df(1)=-xflow*sqrt/(A1*Pt1**2)+
     &                    fact1/pt1*dsqrt(root)
     &                    *(-expon1-expon3*fact2/root)
!     
!     temperature node1
!     
                     df(2)=0.5d0*xflow*dsqrt(R/(kappa*Tt1))/(A1*Pt1)
!     
!     mass flow
!     
                     df(3)=inv*sqrt/(A1*Pt1)
!     
!     pressure node2
!     
                     df(4)=fact1/pt2*dsqrt(root)*
     &                    (expon1+expon3*fact2/root)
!     
                  else
!     
!     section1 is critical
!     
                     f=xflow*dsqrt(Tt1)/pt1-A1*qred_crit
!     
!     pressure node1
!     
                     df(1)=-xflow*dsqrt(Tt1)/pt1**2
!     
!     temperature node1
!     
                     df(2)=0.5d0*xflow/(pt1*dsqrt(Tt1))
!     
!     mass flow
!     
                     df(3)=inv*dsqrt(Tt1)/pt1
!     
!     pressure node2
!     
                     df(4)=0.d0
!     
                  endif
!     
               else
!     
!     section A2 critical
!     
                  call pt2_lim_calc(pt2,pt1,a2,a1,kappa,r,zeta,
     &                 pt2_lim,isothermal)
                  pt1pt2=pt1/pt2_lim
!     
                  fact1=pt1pt2**expon1
!     
                  fact2=pt1pt2**expon2
!     
                  root=2d0/km1*(fact2-1d0)
!     
                  f=xflow*sqrt/(A1*Pt1)-fact1*dsqrt(root)
!     
!     pressure node1
!     
                  df(1)=-xflow*sqrt/(A1*Pt1**2)+
     &                 fact1/pt1*dsqrt(root)
     &                 *(-expon1-expon3*fact2/root)

!     
!     temperature node1
!     
                  df(2)=0.5d0*xflow*dsqrt(R/(kappa*Tt1))/(A1*Pt1)
!     
!     mass flow
!     
                  df(3)=inv*sqrt/(A1*Pt1)
!     
!     pressure node2
!     
                  df(4)=0
!     
               endif
!     
!     case zeta less than zero
!     
            elseif(zeta.lt.0) then
!     
               expon1=-kp1/(zeta*kappa)
               fact1=pt1pt2**expon1
               expon2=km1/(zeta*kappa)
               fact2=pt1pt2**expon2
               expon3=1d0/(zeta*kappa)
               root=2d0/km1*(fact2-1d0) 
!     
               if(Qred1.lt.Qred_crit) then
!     
!     section 1 is not critical
!     
!     residual
!     
                  f=xflow**2*R*Tt1/(A1**2*Pt1**2*Kappa)
     &                 -fact1*root
!     
!     pressure node1
!     
                  df(1)=-2*xflow**2*R*Tt1/(A1**2*Pt1**3*Kappa)
     &                 -1/pt1*fact1*(expon1*root
     &                 +2/(zeta*kappa)*fact2)
!     
!     temperature node1
!     
                  df(2)=xflow**2*R/(A1**2*Pt1**2*Kappa)
!     
!     mass flow
!     
                  df(3)=2*xflow*R*Tt1/(A1**2*Pt1**2*Kappa)
!     
!     pressure node2
!     
                  df(4)=-(1/Pt2*fact1)
     &                 *(-expon1*root-2/(zeta*kappa)*fact2)
!     
!     section1 is critical
!     
               else
!     
                  f=xflow**2*R*Tt1/(A1**2*Pt1**2*Kappa)
     &                 -R/kappa*qred_crit**2
!     
!     pressure node1
!     
                  df(1)=-2*xflow**2*R*Tt1/(A1**2*pt1**3*kappa)
!     
!     temperature node1
!     
                  df(2)=xflow**2*R/(A1**2*Pt1**2*Kappa)
!     
!     mass flow
!     
                  df(3)=2*xflow*R*Tt1/(A1**2*Pt1**2*Kappa)
!     
!     pressure node2
!     
                  df(4)=0.d0
!     
               endif
!     
            endif
!     
         else
!     
!     A1 greater than A2 or equal to A2
!     
            Qred2=dabs(xflow)*dsqrt(Tt2)/(A2*Pt2)
!     
            Qred1=1/pt1pt2*A2/A1*Qred2
!     
            Qred_crit=dsqrt(kappa/R)*(1.d0+0.5d0*km1)
     &           **(-0.5d0*kp1/km1)

!     definition of the coefficients 
!     
            if(zeta.gt.0d0) then
!     
               sqrt=dsqrt(R*Tt2/kappa)
               expon1=-0.5d0*kp1/(zeta*kappa)
               fact1=pt1pt2**expon1
               expon2=km1/(zeta*kappa)
               fact2=pt1pt2**expon2
               expon3=1d0/(zeta*kappa)
               root=2d0/km1*(fact2-1d0)
!     
               if(pt1pt2.ge.pt1pt2_crit) then
                  pt1pt2=pt1pt2_crit
                  pt2=pt1/pt1pt2_crit
               endif
!     
               if((Qred2.lt.Qred_crit)
     &              .and.(pt1/pt2.lt.pt1pt2_crit)) then
!     
!     section 2 is not critical
!     
!     residual
!     
                  f=xflow*sqrt/(A2*Pt2)-fact1*dsqrt(root)
!     
!     pressure node1
!     
                  df(1)=fact1/pt1*dsqrt(root)
     &                 *(-expon1-expon3*fact2/root)
!     
!     temperature node1
!     
                  df(2)=0.5d0*xflow*sqrt/(A2*Pt2*Tt2)
!     
!     mass flow
!     
                  df(3)=inv*sqrt/(A2*Pt2)
!     
!     pressure node2
!     
                  df(4)=-xflow*sqrt/(A2*Pt2**2)+fact1/pt2*dsqrt(root)*
     &                 (expon1+expon3*fact2/root)
!     
!     
               else
!     
!     section2 is critical
!     
                  pt2=pt1/pt1pt2_crit
!     
                  f=xflow*dsqrt(Tt2)/(pt2*A2)-qred_crit
!     
!     pressure node1
!     
                  df(1)=0
!     
!     temperature node1
!     
                  df(2)=0.5d0*xflow/(A2*pt2*dsqrt(Tt2))
!     
!     mass flow
!     
                  df(3)=inv*dsqrt(Tt2)/(A2*pt2)
!     
!     pressure node2
!     
                  df(4)=-xflow*dsqrt(Tt2)/(A2*pt2**2)
!     
               endif
!     
            elseif(zeta.eq.0d0) then
!     
               Qred1=dabs(xflow)*dsqrt(Tt1*kappa/R)/(A1*Pt1)
               Qred2=dabs(xflow)*dsqrt(Tt2*kappa/R)/(A2*Pt2)
               Qred_crit=dsqrt(kappa/R)*(1.d0+0.5d0*km1)
     &              **(-0.5d0*kp1/km1)
!     
               f=pt1/pt2-1.d0
!     
               df(1)=1/pt2
!     
               df(2)=0
!     
               df(3)=0
!     
               df(4)=-pt1/pt2**2
!     
            endif
         endif
!     
      endif
!     
      return
      end
      
      
 
