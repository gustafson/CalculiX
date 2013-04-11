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
      subroutine liquidchan(node1,node2,nodem,nelem,lakon,
     &     nactdog,identity,ielprop,prop,iflag,v,xflow,f,
     &     nodef,idirf,df,rho,g,co,dvi,numf,vold,mi)
!
!     pipe element for incompressible media
!     
      implicit none
!     
      logical identity
      character*8 lakon(*)
!      
      integer nelem,nactdog(0:3,*),node1,node2,nodem,iaxial,
     &     ielprop(*),nodef(4),idirf(4),index,iflag,mi(2),
     &     inv,ncoel,ndi,nbe,id,nen,ngv,numf,nodea,nodeb
!      
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(4),b,d,pi,radius,
     &     h1,h2,rho,dvi,friction,reynolds,vold(0:mi(2),*),dg,
     &     g(3),a1,a2,xn,xk,xk1,xk2,zeta,dl,rh,a0,alpha,a,
     &     coarseness,rd,xks,z1,z2,co(3,*),xcoel(3),yel(3),
     &     yco(3),xdi(10),ydi(10),xbe(7),ybe(7),zbe(7),ratio,
     &     xen(10),yen(10),xgv(8),ygv(8),ds,dd,dfriction,xkn,xkp,
     &     fjump,dfjumpdh1,dfjumpdh2,feq,dfeqdh1,dfeqdh2,xflow2,du1,
     &     dxidh1,angle,b1,b2,sgn,s0,sqrts0,xi,hk
!
      data ncoel /3/
      data xcoel /0.,14.,90./
      data yco /0.,0.06,0.5/
      data yel /0.,0.2,0.5/
!
      data ndi /10/
      data xdi /0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1./
      data ydi /226.,47.5,17.5,7.8,3.75,1.80,0.8,0.29,0.06,0./
!
      data nbe /7/
      data xbe /1.,1.5,2.,3.,4.,6.,10./
      data ybe /0.21,0.12,0.10,0.09,0.09,0.08,0.2/
      data zbe /0.51,0.32,0.29,0.26,0.26,0.17,0.31/
!
      data nen /10/
      data xen /0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1./
      data yen /232.,51.,18.8,9.6,5.26,3.08,1.88,1.17,0.734,0.46/
!
      data ngv /8/
      data xgv /0.125,0.25,0.375,0.5,0.625,0.75,0.875,1./
      data ygv /98.,17.,5.52,2.,0.81,0.26,0.15,0.12/
!
      numf=3
!
!     iflag=0: check whether all parameters in the element equation
!              are known => equation is not needed
!     iflag=1: calculation of the initial flux
!     iflag=2: evaluate the element equation and all derivatives
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
      elseif((iflag.eq.1).or.(iflag.eq.2))then
!     
         index=ielprop(nelem)
!     
         h1=v(2,node1)
         h2=v(2,node2)
!     
         z1=-g(1)*co(1,node1)-g(2)*co(2,node1)-g(3)*co(3,node1)
         z2=-g(1)*co(1,node2)-g(2)*co(2,node2)-g(3)*co(3,node2)
!
         dg=dsqrt(g(1)*g(1)+g(2)*g(2)+g(3)*g(3))
!     
         if(iflag.eq.1) then
            inv=0
         else
            xflow=v(1,nodem)
            if(xflow.ge.0.d0) then
               inv=1
            else
               inv=-1
            endif
            nodef(1)=node1
            nodef(2)=nodem
            nodef(3)=node2
            idirf(1)=2
            idirf(2)=1
            idirf(3)=2
         endif
!     
         if(lakon(nelem)(6:7).eq.'WC') then
!     
!     channel, White-Colebrook
!     
            b=prop(index+1)
            s0=prop(index+2)
            xks=prop(index+3)
            dl=dsqrt((co(1,node2)-co(1,node1))**2+
     &           (co(2,node2)-co(2,node1))**2+
     &           (co(3,node2)-co(3,node1))**2)
!
            if(s0.lt.0.d0) then
               s0=dasin((z1-z2)/dl)
            endif
            sqrts0=dsqrt(1.d0-s0*s0)
!
!           auxiliary variable
!
            xi=dg*(rho*b)**2*(h1+h2)**3
!
            if(iflag.eq.1) then
!
!              assuming large reynolds number
!
               friction=1.d0/(2.03*dlog10(xks/(d*3.7)))**2
            else
!
!              solving the implicit White-Colebrook equation
!
               reynolds=xflow*d/(a*dvi)
               friction=(200.d0*d/(xks*reynolds))**2
               do
                  ds=dsqrt(friction)
                  dd=2.51d0/(reynolds*ds)+xks/(2.7d0*d)
                  dfriction=(1.d0/ds+2.03*dlog10(dd))*2.d0*friction*ds/
     &               (1.d0+2.213d0/(reynolds*dd))
                  if(dfriction.le.friction*1.d-3) then
                     friction=friction+dfriction
                     exit
                  endif
                  friction=friction+dfriction
               enddo
            endif
            if(inv.ne.0) then
               xk=friction*dl
            else
               xkn=friction*dl
               xkp=xkn
            endif
!
!           check whether the critical depth lies in between h1 and h2
!
            if(iflag.gt.1) then
!
!              calculate HK (critical depth)
!               
               hk=((xflow/(rho*b))**2/(dg*sqrts0))**(1.d0/3.d0)
               if(inv.eq.1) then
                  if((h1.lt.hk).and.(h2.gt.hk)) then
                     fjump=((h2-h1)**3)/(4.d0*h1*h2)
                  else
                     fjump=0.d0
                     dfjumpdh1=0.d0
                     dfjumpdh2=0.d0
                  endif
               else
                  if((h1.gt.hk).and.(h2.lt.hk)) then
                     fjump=((h1-h2)**3)/(4.d0*h1*h2)
                  else
                     fjump=0.d0
                     dfjumpdh1=0.d0
                     dfjumpdh2=0.d0
                  endif
               endif
               feq=0.d0
               dfeqdh1=0.d0
               dfeqdh2=0.d0
            endif
         elseif(lakon(nelem)(6:7).eq.'EL') then
!     
!     pipe, sudden enlargement
!     
            b1=prop(index+1)
            b2=prop(index+2)
            angle=prop(index+3)
            call ident(xcoel,angle,ncoel,id)
            dl=0.d0
!
!           calculating the special head loss
!
            if(h2*b2.gt.h1*b1) then
               sgn=-1.d0
            else
               sgn=1.d0
            endif
            du1=((b1+b2)**2*(h1+h2)**3)/8.d0*sgn
            feq=du1*(1.d0/(h2*b2)**2-1.d0/(h1*b1)**2)
!
            if(inv.ne.0) then
               fjump=0.d0
               dfjumpdh1=0.d0
               dfjumpdh2=0.d0
               dfeqdh1=3.d0*feq/(h1+h2)
!
               dfeqdh2=dfeqdh1-2.d0*du1/(h2**3*b2**2)
               dfeqdh1=dfeqdh1+2.d0*du1/(h1**3*b1**2)
            endif
!
            if(inv.ge.0) then
               if(id.eq.0) then
                  zeta=yel(1)
               elseif(id.eq.ncoel) then
                  zeta=yel(ncoel)
               else
                  zeta=yel(id)+(yel(id+1)-yel(id))*(ratio-xcoel(id))/
     &                 (xcoel(id+1)-xcoel(id))
               endif
               if(inv.ne.0) then
                  xk=zeta*feq
               else
                  xkp=zeta*feq
               endif
            endif
            if(inv.le.0) then
               if(id.eq.0) then
                  zeta=yco(1)
               elseif(id.eq.ncoel) then
                  zeta=yco(ncoel)
               else
                  zeta=yco(id)+(yco(id+1)-yco(id))*(ratio-xcoel(id))/
     &                 (xcoel(id+1)-xcoel(id))
               endif
               if(inv.ne.0) then
                  xk=zeta*feq
               else
                  xkn=zeta*feq
               endif
            endif
            
         elseif(lakon(nelem)(6:7).eq.'CO') then
!     
!     pipe, sudden contraction
!     
            a1=prop(index+1)
            a2=prop(index+2)
            ratio=a2/a1
            call ident(xcoel,ratio,ncoel,id)
            if(inv.ge.0) then
               if(id.eq.0) then
                  zeta=yco(1)
               elseif(id.eq.ncoel) then
                  zeta=yco(ncoel)
               else
                  zeta=yco(id)+(yco(id+1)-yco(id))*(ratio-xcoel(id))/
     &                 (xcoel(id+1)-xcoel(id))
               endif
               if(inv.ne.0) then
                  xk=zeta/(a2*a2)
               else
                  xkp=zeta/(a2*a2)
               endif
            endif
            if(inv.le.0) then
               if(id.eq.0) then
                  zeta=yel(1)
               elseif(id.eq.ncoel) then
                  zeta=yel(ncoel)
               else
                  zeta=yel(id)+(yel(id+1)-yel(id))*(ratio-xcoel(id))/
     &                 (xcoel(id+1)-xcoel(id))
               endif
               if(inv.ne.0) then
                  xk=zeta/(a2*a2)
               else
                  xkn=zeta/(a2*a2)
               endif
            endif
         elseif(lakon(nelem)(6:7).eq.'DI') then
!     
!     pipe, diaphragm
!     
            a=prop(index+1)
            a0=prop(index+2)
            a1=a
            a2=a
            ratio=a0/a
            call ident(xdi,ratio,ndi,id)
            if(id.eq.0) then
               zeta=ydi(1)
            elseif(id.eq.ndi) then
               zeta=ydi(ndi)
            else
               zeta=ydi(id)+(ydi(id+1)-ydi(id))*(ratio-xdi(id))/
     &              (xdi(id+1)-xdi(id))
            endif
            if(inv.ne.0) then
               xk=zeta/(a*a)
            else
               xkn=zeta/(a*a)
               xkp=xkn
            endif
         elseif(lakon(nelem)(6:7).eq.'EN') then
!     
!     pipe, entrance
!     
            a=prop(index+1)
            a0=prop(index+2)
            a1=a*1.d10
            a2=a
            ratio=a0/a
            call ident(xen,ratio,nen,id)
            if(id.eq.0) then
               zeta=yen(1)
            elseif(id.eq.nen) then
               zeta=yen(nen)
            else
               zeta=yen(id)+(yen(id+1)-yen(id))*(ratio-xen(id))/
     &              (xen(id+1)-xen(id))
            endif
            if(inv.ne.0) then
               xk=zeta/(a*a)
            else
c
c              to be changed: entrance is different from exit
c
               xkn=zeta/(a*a)
               xkp=xkn
            endif
         elseif(lakon(nelem)(6:7).eq.'GV') then
!     
!     pipe, gate valve
!     
            a=prop(index+1)
            alpha=prop(index+2)
            a1=a
            a2=a
            call ident(xgv,alpha,ngv,id)
            if(id.eq.0) then
               zeta=ygv(1)
            elseif(id.eq.ngv) then
               zeta=ygv(ngv)
            else
               zeta=ygv(id)+(ygv(id+1)-ygv(id))*(alpha-xgv(id))/
     &              (xgv(id+1)-xgv(id))
            endif
            if(inv.ne.0) then
               xk=zeta/(a*a)
            else
               xkn=zeta/(a*a)
               xkp=xkn
            endif
         elseif(lakon(nelem)(6:7).eq.'BE') then
!     
!     pipe, bend
!     
            a=prop(index+1)
            rd=prop(index+2)
            alpha=prop(index+3)
            coarseness=prop(index+4)
            a1=a
            a2=a
            call ident(xbe,rd,nbe,id)
            if(id.eq.0) then
               zeta=ybe(1)+(zbe(1)-ybe(1))*coarseness
            elseif(id.eq.nbe) then
               zeta=ybe(nbe)+(zbe(nbe)-ybe(nbe))*coarseness
            else
               zeta=(1.d0-coarseness)*
     &              (ybe(id)+(ybe(id+1)-ybe(id))*(rd-xbe(id))/
     &              (xbe(id+1)-xbe(id)))
     &              +coarseness*
     &              (zbe(id)+(zbe(id+1)-zbe(id))*(rd-xbe(id))/
     &              (xbe(id+1)-xbe(id)))
            endif
            zeta=zeta*alpha/90.d0
            if(inv.ne.0) then
               xk=zeta/(a*a)
            else
               xkn=zeta/(a*a)
               xkp=xkn
            endif
         endif
!     
         if(iflag.eq.1) then
            xflow=xi*((h2-h1)*sqrts0-dl*s0)/
     &            (8.d0*(h2-h1)-xkp)
            if(xflow.lt.0.d0) then
               xflow=xi*((h2-h1)*sqrts0-dl*s0)/
     &              (8.d0*(h2-h1)+xkn)
               if(xflow.lt.0.d0) then
                  write(*,*) '*WARNING in liquidpipe:'
                  write(*,*) '         initial mass flow could'
                  write(*,*) '         not be determined'
                  write(*,*) '         1.d-10 is taken'
                  xflow=1.d-10
               else
                  xflow=-dsqrt(xflow)
               endif
            else
               xflow=dsqrt(xflow)
            endif
         else
            dxidh1=3.d0*xi/(h1+h2)
            df(1)=((h2-h1)*sqrts0-(dl*s0-inv*fjump))*dxidh1
            df(2)=df(1)
!
            xflow2=xflow*xflow
            du1=xi*sqrts0-8.d0*xflow2
!
            df(1)=df(1)-du1+(xi*dfjumpdh1+xflow2*dfeqdh1)*inv
            df(2)=df(2)+du1+(xi*dfjumpdh2+xflow2*dfeqdh2)*inv
            df(3)=(-8.d0*(h2-h1)+inv*(friction*dl+feq))*2.d0*xflow
            f=(h2-h1)*du1-xi*(dl*s0-inv*fjump)+
     &                    inv*(friction*dl+feq)*xflow2
         endif
!     
      endif
!     
      return
      end
      

