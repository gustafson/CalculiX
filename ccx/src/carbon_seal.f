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
      subroutine carbon_seal(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,iflag,voldgas,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,dvi,numf)
!     
!     carbon seal element calculated with Richter method
!      Richter "Rohrhydraulik", Springer ,1971,p. 175
!     
      implicit none
!     
      logical identity
      character*8 lakon(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(4),idirf(4),index,iflag,
     &     inv,ipkon(*),kon(*)
!
      real*8 prop(*),voldgas(0:3,*),xflow,f,df(4),R,d,l,
     &     p1,p2,T1,cp,physcon(3),dvi,pi,s
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
         d=prop(index+1)
         s=prop(index+2)
         l=prop(index+3)
         pi=4.d0*datan(1.d0)
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
         if(lakon(nelem)(2:6).eq.'CARBS') then
!     
!     gapflow
!     Richter "Rohrhydraulik", Springer ,1971,p. 175
!     
            xflow=inv*Pi*d*s**3*(P1**2-P2**2)/(24.d0*R*T1*dvi*l)
            
         elseif(lakon(nelem)(2:6).ne.'CARBS') then
            write(*,*) '*WARNING in Carbon_seal.f'
            write(*,*) 'unable to perform carbon seal calculation'
            write(*,*) 'check input file'
         endif
!
c         write(*,*) 'xflow',xflow
!     
      elseif (iflag.eq.2)then
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
         d=prop(index+1)
         s=prop(index+2)
         l=prop(index+3)
         pi=4.d0*datan(1.d0)
         
!     
         if (lakon(nelem)(2:8).eq.'CARBS') then
!     
            f=xflow*T1-pi*d*s**3*(P1**2-P2**2)/(24.d0*R*dvi*l)
!     
            df(1)=-(pi*d*s**3*P1)/(12.d0*R*dvi*l)
            df(2)=xflow
            df(3)=T1
            df(4)=(pi*d*s**3*P2)/(12.d0*R*dvi*l)
!       
         endif
      endif
!     
c      do i=1,4
c         write(*,*)'df(',i,')=',df(i)
c      enddo
!     
      return
      end
      

