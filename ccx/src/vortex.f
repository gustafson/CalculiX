!
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine vortex(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,iflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,numf,set,mi)
!     
!     orifice element
!     
      implicit none
!     
      logical identity
      character*8 lakon(*)
      character*81 set(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(4),idirf(4),index,iflag,
     &     inv,ipkon(*),kon(*),t_chang,nelemswirl,mi(*)
!
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(4),kappa,r,cp,
     &     p1,p2,T1,T2,km1,pi,
     &     r2d,r1d,eta,U1,
     &     c1u,c2u, cinput, r1, r2, omega, K1, rpm,ciu,expon,
     &     Ui,Kr,cte1,cte2,qred_crit,A,xflow_oil
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
         kappa=(cp/(cp-R))
         pi=4.d0*datan(1.d0)
         index=ielprop(nelem)
         qred_crit=dsqrt(kappa/R)*
     &        (1+0.5d0*(kappa-1))**(-0.5*(kappa+1)/(kappa-1))
!
!     Because there is no explicit expression relating massflow
!     with to pressure loss for vortices
!     For  FREE as well as for FORCED VORTICES
!     initial mass flow is set to Qred_crit/2 = 0.02021518917
!     with consideration to flow direction
!
         node1=kon(ipkon(nelem)+1)
         node2=kon(ipkon(nelem)+3)
         p1=v(2,node1)
         p2=v(2,node2)
         T1=v(0,node1)
         T2=v(0,node2)
!
!     abstract cross section
         A=10E-6
!
         if(p1.gt.p2) then
            xflow=0.5/dsqrt(T1)*A*P1*qred_crit
         else
            xflow=-0.5/dsqrt(T1)*A*P1*qred_crit
         endif
!        
      elseif (iflag.eq.2)then
!     
         numf=4
         index=ielprop(nelem) 
         kappa=(cp/(cp-R))
         km1=kappa-1
         pi=4.d0*datan(1.d0)
!
!     radius downstream
         r2d=prop(index+1)
!
!     radius upstream
         r1d=prop(index+2)
!
!     pressure correction factor
         eta=prop(index+3)
!     
         p1=v(2,node1)
         p2=v(2,node2)
!     
         xflow=v(1,nodem)
!
         if(xflow.gt.0.d0) then
            inv=1.d0
            p1=v(2,node1)
            p2=v(2,node2)
            T1=v(0,node1)
            T2=v(0,node2)
            R1=r1d
            R2=r2d
!
            nodef(1)=node1
            nodef(2)=node1
            nodef(3)=nodem
            nodef(4)=node2
!
         elseif(xflow.lt.0.d0) then
            inv=-1.d0
            R1=r2d
            R2=r1d
            p1=v(2,node2)
            p2=v(2,node1)
            T1=v(0,node2)
            T2=v(0,node1)
            xflow=-v(1,nodem)
!            
            nodef(1)=node2
            nodef(2)=node2
            nodef(3)=nodem
            nodef(4)=node1
!
         endif
!
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=2
!
         kappa=(cp/(cp-R))
!
!     FREE VORTEX
!
         if(lakon(nelem)(4:5).eq.'FR')then
!
!     rotation induced loss (correction factor)
           K1= prop(index+4)
!
!     tangential velocity of the disk at vortex entry
             U1=prop(index+5)
!
!     number of the element generating the upstream swirl
             nelemswirl=int(prop(index+6))
!     
!     rotation speed (revolution per minutes)
             rpm=prop(index+7)
!
!     Temperature change
            t_chang=prop(index+8)
!
            if(rpm.gt.0) then
!
!     rotation speed is given (rpm) if the swirl comes from a rotating part
!     typically the blade of a coverplate
!
               omega=pi/30d0*rpm
               
!     C_u is given by radius r1d (see definition of the flow direction)
!     C_u related to radius r2d is a function of r1d
!     
               if(inv.gt.0) then
                  c1u=omega*r1
!
!     flow rotation at outlet 
                  c2u=c1u*r1/r2
!
               elseif(inv.lt.0) then
                  c2u=omega*r2
!     
                  c1u=c2u*r2/r1
               endif
!
            elseif(nelemswirl.gt.0) then
               if(lakon(nelemswirl)(2:5).eq.'ORPN') then
                  cinput=prop(ielprop(nelemswirl)+5)
               elseif(lakon(nelemswirl)(2:5).eq.'VOFR') then
                  cinput=prop(ielprop(nelemswirl)+9)
               elseif(lakon(nelemswirl)(2:5).eq.'VOFO') then
                  cinput=prop(ielprop(nelemswirl)+7)
               endif
!
               cinput=U1+K1*(cinput-U1)
!     
               if(inv.gt.0) then
                  c1u=cinput
                  c2u=c1u*R1/R2
               elseif(inv.lt.0) then
                  c2u=cinput
                  c1u=c2u*R2/R1
               endif
            endif
!
!     storing the tengential velocity for later use (wirbel cascade)
            if(inv.gt.0) then
               prop(index+9)=c2u
            elseif(inv.lt.0) then
               prop(index+9)=c1u
            endif  
!
!    inner rotation
!     
            if(R1.lt.R2) then
               ciu=c1u
            elseif(R1.ge.R2) then
               ciu=c2u
            endif
!     
            expon=kappa/km1
!     
            if(R2.ge.R1) then
!     
               cte1=c1u**2/(2*Cp*T1)
               cte2=1-(R1/R2)**2
               
               f=P2/P1-1d0-eta*((1+cte1*cte2)**expon-1d0)
!     
               df(1)=-p2/p1**2
!     
               df(2)=eta*expon*cte1/T1*cte2*
     &              (1+cte1*cte2)**(expon-1)
!     
               df(3)=0
!     
               df(4)=1/p1
!     
            elseif(R2.lt.R1) then 
!     
               cte1=c2u**2/(2*Cp*T2)
               cte2=1-(R2/R1)**2
!     
               f=P1/P2-1d0-eta*((1+cte1*cte2)**expon-1d0)
!     
               df(1)=1/p2
!     
               df(2)=eta*expon*cte1/T1*cte2*
     &           (1+cte1*cte2)**(expon-1)
!     
               df(3)=0
!     
               df(4)=-p1/p2**2
!     
            endif
!
!     FORCED VORTEX
!            
         elseif(lakon(nelem)(4:5).eq.'FO') then
!     
!    core swirl ratio
            Kr=prop(index+4)
!     
!     rotation speed (revolution per minutes) of the rotating part
!     responsible for the swirl
            rpm=prop(index+5)
!
!     Temperature change
            t_chang=prop(index+6)
!          
!     rotation speed
            omega=pi/30*rpm
!
            if(R2.ge.R1) then
               Ui=omega*R1
               c1u=Ui*kr
               c2u=c1u*R2/R1
            elseif(R2.lt.R1) then
               Ui=omega*R2
               c2u=Ui*kr
               c1u=c2u*R1/R2
            endif
!     
!     storing the tengential velocity for later use (wirbel cascade)
            if(inv.gt.0) then
               prop(index+7)=c2u
            elseif(inv.lt.0) then
               prop(index+7)=c1u
            endif   
!
            expon=kappa/km1
!
            if(((R2.ge.R1).and.(xflow.gt.0d0))
     &           .or.((R2.lt.R1).and.(xflow.lt.0d0)))then
!     
               cte1=(c1u)**2/(2*Cp*T1)
               cte2=(R2/R1)**2-1
!     
                f=p2/p1-1-eta*((1+cte1*cte2)**expon-1)
!     
!     pressure node1
               df(1)=-p2/p1**2
!     
!     temperature node1
               df(2)=eta*expon*cte1/T1*cte2*(1+cte1*cte2)**(expon-1)
!    
!     massflow nodem
               df(3)=0
!     
!     pressure node2
               df(4)=1/p1
!
            elseif(((R2.lt.R1).and.(xflow.gt.0d0))
     &              .or.((R2.gt.R1).and.(xflow.lt.0d0)))then
               cte1=(c2u)**2/(2*Cp*T2)
               cte2=(R1/R2)**2-1
!
               f=p1/p2-1-eta*((1+cte1*cte2)**expon-1)
!     
!     pressure node1
               df(1)=1/p2
!     
!     temperature node1
               df(2)=eta*expon*cte1/T2*cte2*(1+cte1*cte2)**(expon-1)
!     
!     massflow nodem
               df(3)=0
!     
!     pressure node2
               df(4)=-p1/p2**2
!
            endif
         endif
!
!     outpout
!
         elseif(iflag.eq.3) then
!
         index=ielprop(nelem) 
         kappa=(cp/(cp-R))
         km1=kappa-1
         pi=4.d0*datan(1.d0)
!
!     radius downstream
         r2d=prop(index+1)
!
!     radius upstream
         r1d=prop(index+2)
!
!     pressure correction factor
         eta=prop(index+3)
!     
         p1=v(2,node1)
         p2=v(2,node2)
!     
         xflow=v(1,nodem)
!
         if(xflow.gt.0.d0) then
            inv=1.d0
            p1=v(2,node1)
            p2=v(2,node2)
            T1=v(0,node1)
            T2=v(0,node2)
            R1=r1d
            R2=r2d
!
            nodef(1)=node1
            nodef(2)=node1
            nodef(3)=nodem
            nodef(4)=node2
!
         elseif(xflow.lt.0.d0) then
            inv=-1.d0
            R1=r2d
            R2=r1d
            p1=v(2,node2)
            p2=v(2,node1)
            T1=v(0,node2)
            T2=v(0,node1)
            xflow=v(1,nodem)
!            
            nodef(1)=node2
            nodef(2)=node2
            nodef(3)=nodem
            nodef(4)=node1
!
         endif
!
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=2
!
         kappa=(cp/(cp-R))
!
!     FREE VORTEX
!
         if(lakon(nelem)(4:5).eq.'FR')then
!
!     rotation induced loss (correction factor)
           K1= prop(index+4)
!
!     tengential velocity of the disk at vortex entry
             U1=prop(index+5)
!
!     number of the element generating the upstream swirl
             nelemswirl=int(prop(index+6))
!     
!     rotation speed (revolution per minutes)
             rpm=prop(index+7)
!
!     Temperature change
            t_chang=prop(index+8)
!
            if(rpm.gt.0) then
!
!     rotation speed is given (rpm) if the swirl comes from a rotating part
!     typically the blade of a coverplate
!
               omega=pi/30d0*rpm
               
!     C_u is given by radius r1d (see definition of the flow direction)
!     C_u related to radius r2d is a function of r1d
!     
               if(inv.gt.0) then
                  c1u=omega*r1
!
!     flow rotation at outlet 
                  c2u=c1u*r1/r2
!
               elseif(inv.lt.0) then
                  c2u=omega*r2
!     
                  c1u=c2u*r2/r1
               endif
!
            elseif(nelemswirl.gt.0) then
               if(lakon(nelemswirl)(2:5).eq.'ORPN') then
                  cinput=prop(ielprop(nelemswirl)+5)
               elseif(lakon(nelemswirl)(2:5).eq.'VOFR') then
                  cinput=prop(ielprop(nelemswirl)+9)
               elseif(lakon(nelemswirl)(2:5).eq.'VOFO') then
                  cinput=prop(ielprop(nelemswirl)+7)
               endif
!
               cinput=U1+K1*(cinput-U1)
!     
               if(inv.gt.0) then
                  c1u=cinput
                  c2u=c1u*R1/R2
               elseif(inv.lt.0) then
                  c2u=cinput
                  c1u=c2u*R2/R1
               endif
            endif
!
!     storing the tengential velocity for later use (wirbel cascade)
            if(inv.gt.0) then
               prop(index+9)=c2u
            elseif(inv.lt.0) then
               prop(index+9)=c1u
            endif   
!
!    inner rotation
!     
            if(R1.lt.R2) then
               ciu=c1u
            elseif(R1.ge.R2) then
               ciu=c2u
            endif
!     
            expon=kappa/km1
!     
            if(R2.ge.R1) then
!     
               cte1=c1u**2/(2*Cp*T1)
               cte2=1-(R1/R2)**2
               
               f=P2/P1-1d0-eta*((1+cte1*cte2)**expon-1d0)
!     
               df(1)=-p2/p1**2
!     
               df(2)=eta*expon*cte1/T1*cte2*
     &              (1+cte1*cte2)**(expon-1)
!     
               df(3)=0
!     
               df(4)=1/p1
!     
            elseif(R2.lt.R1) then 
!     
               cte1=c2u**2/(2*Cp*T2)
               cte2=1-(R2/R1)**2
!     
               f=P1/P2-1d0-eta*((1+cte1*cte2)**expon-1d0)
!     
               df(1)=1/p2
!     
               df(2)=eta*expon*cte1/T1*cte2*
     &           (1+cte1*cte2)**(expon-1)
!     
               df(3)=0
!     
               df(4)=-p1/p2**2
!     
            endif
!
!     FORCED VORTEX
!            
         elseif(lakon(nelem)(4:5).eq.'FO') then
!     
!    core swirl ratio
            Kr=prop(index+4)
!     
!     rotation speed (revolution per minutes) of the rotating part
!     responsible for the swirl
            rpm=prop(index+5)
!
!     Temperature change
            t_chang=prop(index+6)
!
!    no element generating the upstream swirl
             nelemswirl=0
!          
!     rotation speed
            omega=pi/30*rpm
!     
            if(R2.ge.R1) then
               Ui=omega*R1
               c1u=Ui*kr
               c2u=c1u*R2/R1
            elseif(R2.lt.R1) then
               Ui=omega*R2
               c2u=Ui*kr
               c1u=c2u*R1/R2
            endif
!     
!     storing the tengential velocity for later use (wirbel cascade)
            if(inv.gt.0) then
               prop(index+7)=c2u
            elseif(inv.lt.0) then
               prop(index+7)=c1u
            endif   
!
            expon=kappa/km1
         endif
!
         xflow_oil=0.d0
!
         write(1,*) ''
         write(1,55) 'In line',int(nodem/1000),' from node',node1,
     &' to node', node2,':   air massflow rate=',xflow,'kg/s',
     &', oil massflow rate=',xflow_oil,'kg/s'
 55      FORMAT(1X,A,I6.3,A,I6.3,A,I6.3,A,F9.5,A,A,F9.5,A)

         if(inv.eq.1) then
            write(1,56)'       Inlet node ',node1,':     Tt1= ',T1,
     &           'K, Ts1= ',T1,'K, Pt1= ',P1/1E5,
     &           'Bar'
            write(1,*)'             element V    ',set(numf)(1:20)
            write(1,57)'             C1u= ',C1u,'m/s ,C2u= ',C2u,'m/s'
            write(1,56)'       Outlet node ',node2,':    Tt2= ',T2,
     &           'K, Ts2= ',T2,'K, Pt2= ',P2/1e5,
     &           'Bar'
!     
         else if(inv.eq.-1) then
            write(1,56)'       Inlet node ',node2,':     Tt1= ',T1,
     &           'K, Ts1= ',T1,'K, Pt1= ',P1/1E5,
     &           'Bar'
            write(1,*)'             element V    ',set(numf)(1:20)
            write(1,57)'             C1u= ',C1u,'m/s ,C2u= ',C2u,'m/s'
            write(1,56)'       Outlet node ',node1,'     Tt2= ',
     &           T2,'K, Ts2= ',T2,'K, Pt2= ',P2/1e5,
     &           'Bar'
         endif
 56      FORMAT(1X,A,I6.3,A,f6.1,A,f6.1,A,f9.5,A,f9.5)  
 57      FORMAT(1X,A,f6.2,A,f6.2,A)
      endif
!      
      return
      end
      
























                  



