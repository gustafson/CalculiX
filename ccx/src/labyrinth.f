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
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!     
      subroutine labyrinth(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,iflag,voldgas,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,dvi,numf)
!     
!     labyrinth element
!     
      implicit none
!     
      logical identity
      character*8 lakon(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(4),idirf(4),index,iflag,
     &     inv,ipkon(*),kon(*),kgas,n
!
      real*8 prop(*),voldgas(0:3,*),xflow,f,df(4),kappa,R,a,d,
     &     p1,p2,T1,Aeff,C1,C2,C3,cd,cp,physcon(3),p2p1,km1,dvi,
     &     kp1,kdkm1,tdkp1,km1dk,x,y,ca1,cb1,ca2,cb2,dT1,alambda,
     &     rad,reynolds,pi,ppkrit,
     &     carry_over,lc,hst,e,szt,num,denom,t,s,b,h,cdu,
     &     cd_radius,cst,dh,cd_honeycomb,cd_lab,bdh,
     &     pt0zps1,cd_1spike,cdbragg,rzdh,
     &     cd_correction,p1p2
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
         index=ielprop(nelem)
         kappa=(cp/(cp-R))
         t=prop(index+1)
         s=prop(index+2)
         d=prop(index+3)
         n=int(prop(index+4))
         b=prop(index+5)
         h=prop(index+6)
!     hc=prop(index+7)
         lc=prop(index+7)
         rad=prop(index+8)
         X=prop(index+9)
         Hst=prop(index+10)
!     
         pi=4.d0*datan(1.d0)
         A=pi*D*s
         e=2.718281828459045d0
!     
         p1=voldgas(2,node1)
         p2=voldgas(2,node2)
         if(p1.ge.p2) then
            inv=1
            T1=voldgas(0,node1)+physcon(1)
         else
            inv=-1
            p1=voldgas(2,node2)
            p2=voldgas(2,node1)
            T1=voldgas(0,node2)+physcon(1)
         endif
!     
         cd=1.d0
         Aeff=A*cd
         p2p1=p2/p1
!     
!************************
!     one fin 
!*************************
         if(n.eq.1.d0) then
!     
            km1=kappa-1.d0
            kp1=kappa+1.d0
            kdkm1=kappa/km1
            tdkp1=2.d0/kp1
            C2=tdkp1**kdkm1
!     
!     subcritical
!     
            if(p2p1.gt.C2) then
               xflow=inv*p1*Aeff*dsqrt(2.d0*kdkm1*p2p1**(2.d0/kappa)
     &              *(1.d0-p2p1**(1.d0/kdkm1))/r)/dsqrt(T1)
!     
!     critical
!     
            else
               xflow=inv*p1*Aeff*dsqrt(kappa/r)*tdkp1**(kp1/(2.d0*km1))/
     &              dsqrt(T1)
            endif
c            write(*,*) 'xflow=',xflow
         endif
!     
!***********************
!     straight labyrinth and stepped labyrinth
!     method found in "Air system Correlations Part1 Labyrinth Seals"
!     H.Zimmermann and K.H. Wolff
!     ASME 98-GT-206
!**********************
!     
         if (n.ge.2) then
!     
            call lab_straight_ppkrit(n,ppkrit)
!     
!     subcritical case
!     
            if (p2p1.gt.ppkrit) then
               xflow=inv*p1*Aeff/dsqrt(T1)*dsqrt((1.d0-p2p1**2.d0)
     &              /(R*(n-log(p2p1)/log(e))))
!     
!     critical case
!     
            else
               xflow=inv*p1*Aeff/dsqrt(T1)*dsqrt(2.d0/R)*ppkrit
            endif
c            write(*,*) 'xflow=',xflow
         endif
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif (iflag.eq.2)then
!     
         alambda=10000.d0
!     
         p1=voldgas(2,node1)
         p2=voldgas(2,node2)
         if(p1.ge.p2) then
            inv=1
            xflow=voldgas(1,nodem)
            T1=voldgas(0,node1)+physcon(1)
            nodef(1)=node1
            nodef(2)=node1
            nodef(3)=nodem
            nodef(4)=node2
         else
            inv=-1
            p1=voldgas(2,node2)
            p2=voldgas(2,node1)
            xflow=-voldgas(1,nodem)
            T1=voldgas(0,node2)+physcon(1)
            nodef(1)=node2
            nodef(2)=node2
            nodef(3)=nodem
            nodef(4)=node1
         endif
!     
c     xflow=voldgas(1,nodem)
!     
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=2
!     
         index=ielprop(nelem)
         kappa=(cp/(cp-R))
         t=prop(index+1)
         s=prop(index+2)
         d=prop(index+3)
         n=int(prop(index+4))
         b=prop(index+5)
         h=prop(index+6)
!     hc=prop(index+7)
         lc=prop(index+7)
         rad=prop(index+8)
         X=prop(index+9)
         Hst=prop(index+10)
!     
         p2p1=p2/p1
         dT1=dsqrt(T1)
!     
         pi=4.d0*datan(1.d0)
         A=pi*D*s
         Aeff=A
         e=2.718281828459045d0
!     
!     honeycomb stator correction
!     
         cd_honeycomb=1.d0
         if (lc.ne.0.d0)then
            call cd_lab_honeycomb(s,lc,cd_honeycomb)
c            write(*,*) '1+cd_honeycomb/100',1+cd_honeycomb/100
            cd_honeycomb=1+cd_honeycomb/100
!            Aeff=Aeff*(1.d0+cd_honeycomb/100.d0)
         endif
!     
!     inlet radius correction
!     
         cd_radius=1.d0
         if((rad.ne.0.d0).and.(n.ne.1d0)) then
            call cd_lab_radius(rad,s,Hst,cd_radius)
c            write(*,*) 'cd_radius',cd_radius
!            Aeff=Aeff*cd_radius
         endif
!     
!     carry over factor (only for straight throught labyrinth)
!     
         if ((n.ge.2).and.(hst.eq.0d0)) then
            cst=n/(n-1.d0)
            szt=s/t
            carry_over=cst/dsqrt(cst-szt/(szt+0.02))
c            write(*,*) 'carry_over=',carry_over
            Aeff=Aeff*carry_over
         endif
!     
!     calculation of the dynamic viscosity 
!     
         if(dabs(dvi).lt.1E-30) then
            kgas=0
            call dynamic_viscosity(kgas,T1,dvi)
         endif     
!     
!     calculation of the number of reynolds for a gap
!     
         reynolds=dabs(xflow)*2.d0*s/(dvi*A*cd_honeycomb/cd_radius)
c         write(*,*) ''
c         write(*,*) 'dynamic viscosity=',dvi
c         write(*,*) 'reynolds=',reynolds
!     
!**************************************
!     single fin labyrinth 
!     the resolution procedure is the same as for the restrictor
!**************************************
!     
         if(n.eq.1)then
!     
!     single fin labyrinth
!     
!     incompressible basis cd , reynolds correction,and radius correction
!     
!     "Flow Characteristics of long orifices with rotation and corner radiusing"
!     W.F. Mcgreehan and M.J. Schotsch
!     ASME 87-GT-162
!     
            dh=2*s
            bdh=b/dh
            rzdh=rad/dh
!     
            call cd_Mcgreehan_Schotsch(rzdh,bdh,reynolds,cdu) 
c            write(*,*) 'McGreehan incompressible cd:cdu=',cdu
!     
!     compressibility correction factor
!     
!     S.L.Bragg
!     "Effect of conpressibility on the discharge coefficient of orifices and convergent nozzles"
!     Journal of Mechanical engineering vol 2 No 1 1960
!     
            call cd_bragg(cdu,p2p1,cdbragg)
            cd=cdbragg
c            write(*,*) 'Bragg correction cd=',cdbragg
c            write(*,*) ''
            Aeff=Aeff*cd 
!     
            km1=kappa-1.d0
            kp1=kappa+1.d0
            kdkm1=kappa/km1
            tdkp1=2.d0/kp1
            C2=tdkp1**kdkm1
!     
            if(p2p1.gt.C2) then
               C1=dsqrt(2.d0*kdkm1/r)*Aeff
               km1dk=1.d0/kdkm1
               y=p2p1**km1dk
               x=dsqrt(1.d0-y)
               ca1=-C1*x/(kappa*p1*y)
               cb1=C1*km1dk/(2.d0*p1)
               ca2=-ca1*p2p1-xflow*dT1/(p1*p1)
               cb2=-cb1*p2p1
               f=xflow*dT1/p1-C1*p2p1**(1.d0/kappa)*x
               if(cb2.le.-(alambda+ca2)*x) then
                  df(1)=-alambda
               elseif(cb2.ge.(alambda-ca2)*x) then
                  df(1)=alambda
               else
                  df(1)=ca2+cb2/x
               endif
               df(2)=xflow/(2.d0*p1*dT1)
               df(3)=inv*dT1/p1
               if(cb1.le.-(alambda+ca1)*x) then
                  df(4)=-alambda
               elseif(cb1.ge.(alambda-ca1)*x) then
                  df(4)=alambda
               else
                  df(4)=ca1+cb1/x
               endif
            else
               C3=dsqrt(kappa/r)*(tdkp1)**(kp1/(2.d0*km1))*Aeff
               f=xflow*dT1/p1-C3
               df(1)=-xflow*dT1/(p1)**2
               df(2)=xflow/(2*p1*dT1)
               df(3)=inv*dT1/p1
               df(4)=0.d0
            endif
         endif
!
!     
!****************************************
!     straight labyrinth & stepped labyrinth
!     method found in "Air system Correlations Part1 Labyrinth Seals"
!     H.Zimmermann and K.H. Wolff
!     ASME 98-GT-206
!****************************************
!     
         if(n.ge.2) then
            num=(1.d0-p2p1**2)
            denom=R*(n-log(p2p1)/log(e))
!     
!     straight labyrinth
!     
            if((hst.eq.0.d0).and.(n.ne.1)) then
               call cd_lab_straight(n,p2p1,s,b,reynolds,cd_lab)
c               write(*,*) 'cd_straight=',cd_lab
               Aeff=Aeff*cd_lab*cd_honeycomb*cd_radius
!     
!     Stepped Labyrinth
!     
            else 
!     corrective term for the first spike
               p1p2=p1/p2
               pt0zps1=(p1p2)**(1/prop(index+4))
c               write(*,*) 'p0/pn=',p1/p2
c               write(*,*) 'p0/p1=',pt0zps1
               call cd_lab_1spike (pt0zps1,s,b,cd_1spike)
c               write(*,*) 'cd1=',cd_1spike
!     
!     corrective term for cd_lab_1spike
!     
               call cd_lab_correction (p1p2,s,b,cd_correction)
c               write(*,*) 'cd_correction=', cd_correction
!     
!     calculation of the discharge coefficient of the stepped labyrinth
!     
               cd=cd_1spike*cd_correction
               cd_lab=cd
!     
c               write(*,*) 'cd_lab=',cd_lab
c               write(*,*) ''
!     
               Aeff=Aeff*cd_lab*cd_radius*cd_honeycomb
            endif
!     
            call lab_straight_ppkrit(n,ppkrit)
c            write(*,*) 'ppkrit=',ppkrit
!     
!     subcritical case
!     
            if (p2p1.gt.ppkrit) then
!     
               f=xflow*dT1/p1-dsqrt(num/denom)*Aeff
!     
               df(1)=xflow*dt1/p1**2.d0-Aeff/2.d0
     &              *dsqrt(denom/num)*(2.d0*(p2**2.d0/p1**3.d0)/denom)
     &              +num/denom**2.d0*r/p1
               df(2)=xflow/(2.d0*p1*dT1)
               df(3)=inv*dT1/p1
               df(4)=-Aeff/2.d0*dsqrt(denom/num)*(-2.d0*(p2/p1**2.d0)
     &              /denom)+num/denom**2.d0*r/p2
!     
!     critical case
!     
            else
               C2=dsqrt(2/R)*Aeff*ppkrit
!     
               f=xflow*dT1/p1-C2
               df(1)=-xflow*dT1/(p1**2)
               df(2)=xflow/(2.d0*p1*dT1)
               df(3)=inv*dT1/p1
               df(4)=0.d0
            endif
         endif
      endif
!         
      return
      end
      
      
