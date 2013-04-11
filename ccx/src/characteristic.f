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
      subroutine  characteristic(node1,node2,nodem,nelem,lakon,
     &     nactdog,identity,ielprop,prop,iflag,voldgas,xflow,f,
     &     nodef,idirf,df,cp,r,physcon,dvi,numf)
!     
!     This subroutine is used to enables the processing of empiric 
!     given under the form
!     massflow*dsqrt(T1)/Pt1=f((Pt1-Pt2)/Pt1) and T1=T2
!     characteristics the subroutine proceeds using 
!     linear interpolation to estimate the values for the whole characteristic
!     note that the characteristic is implicitely containing the point (0,0)
!     
      implicit none
!     
      logical identity
      character*8 lakon(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,
     &     ielprop(*),nodef(4),idirf(4),index,iflag,
     &     inv,id,numf,npu,i
!     
      real*8 prop(*),voldgas(0:3,*),xflow,f,df(4),R,
     &     p1,p2,cp,physcon(3),dvi,
     &     xpu(10),ypu(10),Qred,p1mp2zp1,T1,scal
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
      elseif ((iflag.eq.1).or.(iflag.eq.2)) then
!     
         index=ielprop(nelem)
!     
         npu=nint(prop(index+1))
         scal=prop(index+2)
!         scal=prop(index+npu)
         do i=1,npu
            xpu(i)=prop(index+2*i)
            ypu(i)=prop(index+2*i+1)
         enddo
!     
         p1=voldgas(2,node1)
         p2=voldgas(2,node2)       
!     
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
         p1mp2zp1=(P1-P2)/P1
!     
         if(iflag.eq.1) then
            
            call ident(xpu,p1mp2zp1,npu,id)
            if(id.le.2) then
               Qred=scal*ypu(2)/xpu(2)*p1mp2zp1
               xflow=inv*Qred*P1/dsqrt(T1)
            elseif(id.ge.npu) then
               Qred=scal*ypu(npu-2)
               xflow=inv*Qred*P1/dsqrt(T1)
            else
               Qred=scal*ypu(id)+(ypu(id+1)-ypu(id))
     &             *(p1mp2zp1-xpu(id))/(xpu(id+1)-xpu(id))
               xflow=inv*Qred*P1/dsqrt(T1)
            endif
!     
            write(*,*) ''
            write(*,*) 'Qred=',Qred
            write(*,*) 'xflow=',xflow
            write(*,*) ''
            
         elseif (iflag.eq.2) then
!     
!     index=ielprop(nelem)
!     
!     npu=nint(prop(index+1))
!     do i=1,npu
!     xpu(i)=prop(index+2*i)
!     ypu(i)=prop(index+2*i+1)
!     enddo
!     
            p1=voldgas(2,node1)
            p2=voldgas(2,node2) 
            xflow=voldgas(1,nodem)
!     
            if (p1.ge.p2) then
!     
               inv=1
               xflow=voldgas(1,nodem)
               T1=voldgas(0,node1)+physcon(1)
               nodef(1)=node1
               nodef(2)=node1
               nodef(3)=nodem
               nodef(4)=node2
!     
            else 
!     
               inv=-1
               p1=voldgas(2,node2)
               p2=voldgas(2,node1) 
               T1=voldgas(0,node2)+physcon(1)
               xflow=-voldgas(1,nodem)
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
            df(2)=xflow/(2.d0*P1*dsqrt(T1))
            df(3)=inv*dsqrt(T1)/P1
!     
            call ident(xpu,p1mp2zp1,npu,id)
            write(*,*)'id=',id
!     
            if(id.lt.2) then
               f=dabs(xflow)*dsqrt(T1)/p1-scal*ypu(2)/xpu(2)*p1mp2zp1
               df(4)=scal*ypu(2)/(xpu(2)*P1)
               df(1)=-xflow*dsqrt(T1)/(P1**2.d0)-(P2/P1**2.d0)
     &              *scal*ypu(2)/xpu(2)
!     
            elseif(id.ge.npu) then
               f=dabs(xflow)/P1*dsqrt(T1)-scal*ypu(npu)
               df(4)=0.01d0
               df(1)=-xflow*dsqrt(T1)/P1**2
!    
            else
               f=dabs(xflow)/P1*dsqrt(T1)-(scal*ypu(id)
     &             +scal*(ypu(id+1)-ypu(id))
     &              *(p1mp2zp1-xpu(id))/(xpu(id+1)-xpu(id)))
!
               df(4)=scal*(ypu(id+1)-ypu(id))/(xpu(id+1)-xpu(id))*1/p1
!
               df(1)=-xflow*dsqrt(T1)/P1**2-(P2/P1**2)
     &              *(scal*(ypu(id+1)-ypu(id))/(xpu(id+1)-xpu(id)))
            endif
         endif
!     
      endif
!     
      return
      end
      
      
      
