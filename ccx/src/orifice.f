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
      subroutine orifice(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,iflag,voldgas,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,dvi,numf)
!     
!     orifice element
!     
      implicit none
!     
      logical identity
      character*8 lakon(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(4),idirf(4),index,iflag,
     &     inv,ipkon(*),kon(*),number
!
      real*4 ofvidg
!     
      real*8 prop(*),voldgas(0:3,*),xflow,f,df(4),kappa,R,a,d,xl,
     &     p1,p2,T1,Aeff,C1,C2,C3,cd,cp,physcon(3),p2p1,km1,dvi,
     &     kp1,kdkm1,tdkp1,km1dk,x,y,ca1,cb1,ca2,cb2,dT1,alambda,
     &     rad,beta,reynolds,theta,k_phi,c2u_new,u,pi,
     &     ps1pt1,uref,nelemref,cd_chamf,angle,vid
!
      external ofvidg
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
         a=prop(index+1)
         d=prop(index+2)
         xl=prop(index+3)
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
!     
         p2p1=p2/p1
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         tdkp1=2.d0/kp1
         C2=tdkp1**kdkm1
         Aeff=A*cd
         if(p2p1.gt.C2) then
            xflow=inv*p1*Aeff*dsqrt(2.d0*kdkm1*p2p1**(2.d0/kappa)
     &           *(1.d0-p2p1**(1.d0/kdkm1))/r)/dsqrt(T1)
         else
            xflow=inv*p1*Aeff*dsqrt(kappa/r)*tdkp1**(kp1/(2.d0*km1))/
     &           dsqrt(T1)
         endif
!     
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
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=2
!     
         index=ielprop(nelem)
         kappa=(cp/(cp-R))
         a=prop(index+1)
!     
         if ((lakon(nelem)(4:5).ne.'BT').and.
     &        (lakon(nelem)(4:5).ne.'PN').and.
     &        (lakon(nelem)(4:5).ne.'C1')) then
            d=prop(index+2)
            xl=prop(index+3)
            u=prop(index+7)
            nelemref=nint(prop(index+8))
            if (nelemref.eq.0) then
               uref=0.d0
            else
               if(lakon(nelemref)(2:5).eq.'ORPN') then
                  uref=voldgas(1,kon(ipkon(nelemref)+3))
               else
                  write(*,*) '*ERROR in orifice:'
                  write(*,*) ' element',nelemref
                  write(*,*) 'refered by element',nelem
                  write(*,*) 'is not a preswirl nozzle'
               endif
            endif
            u=u-uref
            angle=prop(index+5)
!     
         endif
!     
         if (lakon(nelem)(2:5).eq.'ORMA') then
!
!     calculate the discharge coefficient using own table data and 
!     using Dr.Albers method for rotating cavities
!     
            call cd_own_albers(p1,p2,xl,d,cd,u,T1,R,kappa)
!     
            write(*,*) 'cd OWN/Albers equals to',cd
!     
!     chamfer correction
!     
            if(angle.gt.0.d0)then
               call cd_chamfer(xl,d,p1,p2,angle,cd_chamf)
               cd=cd*cd_chamf
            endif
!     
         elseif (lakon(nelem)(2:5).eq.'ORMM') then
!
!     calculate the discharge coefficient using McGreehan and Schotsch method
!
            rad=prop(index+4)
!     
            reynolds=xflow*d/(dvi*a)
!     
            call cd_ms_ms(p1,p2,T1,rad,d,xl,kappa,r,reynolds,u,vid,cd)
!     
            if (cd.ge.1) then
               write(*,*) ''
               write(*,*) '**WARNING**'
               write(*,*) 'in RESTRICTOR ',nelem
               write(*,*) 'Calculation using' 
               write(*,*) ' McGreehan and Schotsch method:'
               write(*,*) ' Cd=',Cd,'>1 !'
               write(*,*) 'Calcultion will proceed will Cd=1'
               write(*,*) 'l/d=',xl/d,'r/d=',rad/d,'u/vid=',u/vid
               cd=1.d0
            endif
!     
c            write(*,*) 'cd McGreehan & Schotsch equals to',cd
!     
!     chamfer correction
!     
            if(angle.gt.0.d0) then
               call cd_chamfer(xl,d,p1,p2,angle,cd_chamf)
               cd=cd*cd_chamf
            endif
!     
         elseif  (lakon(nelem)(2:5).eq.'ORPA') then
!     
!     calculate the discharge coefficient using Parker and Kercher method
!     and using Dr. Albers method for rotating cavities
!     
            rad=prop(index+4)
!     
            beta=prop(index+6)
!     
            reynolds=xflow*d/(dvi*a)
!     
            call cd_pk_albers(rad,d,xl,reynolds,p2,p1,beta,kappa,
     &           cd,u,T1,R)
c            write(*,*) 'cd Parker & Kercher equals to',cd
!     
!     chamfer correction
!     
            if(angle.gt.0.d0) then
               call cd_chamfer(xl,d,p1,p2,angle,cd_chamf)
               cd=cd*cd_chamf
            endif
!     
         elseif (lakon(nelem)(2:5).eq.'ORPM') then
!     
!     calculate the discharge coefficient using Parker and Kercher method
!     and using Mac Grehan and Schotsch method for rotating cavities
!     
            rad=prop(index+4)
!     
            beta=prop(index+6)
            reynolds=xflow*d/(dvi*a)
!     
            call cd_pk_ms(rad,d,xl,reynolds,p2,p1,beta,kappa,cd,
     &           u,T1,R)
!     
!     chamfer correction
!     
            if(angle.gt.0.d0) then
               call cd_chamfer(xl,d,p1,p2,angle,cd_chamf)
               cd=cd*cd_chamf
            endif
!     
         elseif (lakon(nelem)(2:5).eq.'ORC1') then
!     
!     cd=1.d0
!     
            cd=1.d0
!     
         elseif (lakon(nelem)(2:5).eq.'ORBT') then
!     
!     calculate the discharge coefficient of bleed tappings (OWN tables)
!     
            ps1pt1=prop(index+2)
            number=nint(prop(index+3))
!     
            call cd_bleedtapping(p2,p1,ps1pt1,number,cd)
!     
         elseif (lakon(nelem)(2:5).eq.'ORPN') then
!     
!     calculate the discharge coefficient of preswirl nozzle (OWN tables)
!     
            number=nint(prop(index+3))
            call cd_preswirlnozzle(p2,p1,number,cd)
            theta=prop(index+2)
            k_phi=prop(index+3)
!     
            if(p2/p1.gt.(2/(kappa+1.d0))**(kappa/(kappa-1.d0))) then
               c2u_new=k_phi*cd*sin(theta*Pi/180.d0)*r*
     &              dsqrt(2.d0*kappa/(r*(kappa-1)))*
     &              dsqrt(T1*(1.d0-(p2/p1)**((kappa-1)/kappa)))
!     
            else
               c2u_new=k_phi*cd*sin(theta*Pi/180.d0)*r*
     &              dsqrt(2.d0*kappa/(r*(kappa-1)))*
     &              dsqrt(T1*(1.d0-2/(kappa+1)))
            endif
!     
c            write(*,*) 'tangential velocity is',c2u_new
c            write(*,*) 'at pre swirl nozzle exit'
         endif
!     
c         write(*,*)
c         write(*,*) 'cd',cd
c         write(*,*) 'dvi',dvi
c         write(*,*) 'reynolds=',reynolds
!     
         if (cd.gt.1.d0) then
            Write(*,*) '*WARNING:'
            Write(*,*) 'In RESTRICTOR',nelem
            write(*,*) 'Cd greater than 1'
            write (*,*) 'Calculation will proceed using Cd=1'
            cd=1.d0
         endif
         pi=4.d0*datan(1.d0)   
         p2p1=p2/p1
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         tdkp1=2.d0/kp1
         C2=tdkp1**kdkm1
         Aeff=A*cd
         dT1=dsqrt(T1)
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
      return
      end
      

