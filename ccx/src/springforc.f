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
      subroutine springforc(xl,konl,vl,imat,elcon,nelcon,
     &  elas,fnl,ncmat_,ntmat_,nope,lakonl,t0l,t1l,kode,elconloc,
     &  plicon,nplicon,npmat_,veoldl,senergy,iener,cstr,mi,
     &  springarea,nmethod)
!
!     calculates the force of the spring
!
      implicit none
!
      character*8 lakonl
!
      integer konl(9),i,j,imat,ncmat_,ntmat_,nope,nterms,iflag,mi(2),
     &  kode,niso,id,nplicon(0:ntmat_,*),npmat_,nelcon(2,*),iener,
     &  nmethod
!
      real*8 xl(3,9),elas(21),ratio(9),t0l,t1l,al(3),vl(0:mi(2),9),
     &  pl(3,9),xn(3),dm,alpha,beta,fnl(3,9),
     &  veoldl(0:mi(2),9),dist,c2,c3,t(3),dt,
     &  elcon(0:ncmat_,ntmat_,*),pproj(3),xsj2(3),xs2(3,7),val,
     &  shp2(7,8),xi,et,elconloc(21),plconloc(82),xk,fk,dd,
     &  xiso(20),yiso(20),dd0,plicon(0:2*npmat_,ntmat_,*),
     &  um,eps,pi,senergy,cstr(6),
     &  springarea
!
      data iflag /2/
!
!     actual positions of the nodes belonging to the contact spring
!
      do i=1,nope
         do j=1,3
            pl(j,i)=xl(j,i)+vl(j,i)
         enddo
      enddo
!
      if(lakonl(7:7).eq.'A') then
         dd0=dsqrt((xl(1,2)-xl(1,1))**2
     &           +(xl(2,2)-xl(2,1))**2
     &           +(xl(3,2)-xl(3,1))**2)
         dd=dsqrt((pl(1,2)-pl(1,1))**2
     &           +(pl(2,2)-pl(2,1))**2
     &           +(pl(3,2)-pl(3,1))**2)
         do i=1,3
            xn(i)=(pl(i,2)-pl(i,1))/dd
         enddo
         val=dd-dd0
!
!        interpolating the material data
!
         call materialdata_sp(elcon,nelcon,imat,ntmat_,i,t0l,t1l,
     &     elconloc,kode,plicon,nplicon,npmat_,plconloc,ncmat_)
!
!        calculating the spring force and the spring constant
!
         if(kode.eq.2)then
            xk=elconloc(1)
            fk=xk*val
            if(iener.eq.1) then
               senergy=fk*val/2.d0
            endif
         else
            niso=int(plconloc(81))
            do i=1,niso
               xiso(i)=plconloc(2*i-1)
               yiso(i)=plconloc(2*i)
            enddo
            call ident(xiso,val,niso,id)
            if(id.eq.0) then
               xk=0.d0
               fk=yiso(1)
               if(iener.eq.1) then
                  senergy=fk*val;
               endif
            elseif(id.eq.niso) then
               xk=0.d0
               fk=yiso(niso)
               if(iener.eq.1) then
                  senergy=yiso(1)*xiso(1)
                  do i=2,niso
                     senergy=senergy+(xiso(i)-xiso(i-1))*(yiso(i)+yiso(
     &               i-1))/2.d0
                  enddo
                  senergy=senergy+(val-xiso(niso))*yiso(niso)
               endif
            else
               xk=(yiso(id+1)-yiso(id))/(xiso(id+1)-xiso(id))
               fk=yiso(id)+xk*(val-xiso(id))
               if(iener.eq.1) then
                  senergy=yiso(1)*xiso(1)
                  do i=2, id
                     senergy=senergy+(xiso(i)-xiso(i-1))*
     &                    (yiso(i)+yiso(i-1))/2.d0
                  enddo
                  senergy=senergy+(val-xiso(id))*(fk+yiso(id))/2.d0
               endif
            endif
         endif
!
         do i=1,3
            fnl(i,1)=-fk*xn(i)
            fnl(i,2)=fk*xn(i)
         enddo
         return
      endif
!
      nterms=nope-1
!
!     vector vr connects the dependent node with its projection
!     on the independent face
!
      do i=1,3
         pproj(i)=pl(i,nope)
      enddo
c      write(*,*) 'springforc ',(pproj(i),i=1,3)
      call attach(pl,pproj,nterms,ratio,dist,xi,et)
      do i=1,3
         al(i)=pl(i,nope)-pproj(i)
      enddo
!
!     determining the jacobian vector on the surface 
!
      if(nterms.eq.8) then
         call shape8q(xi,et,pl,xsj2,xs2,shp2,iflag)
      elseif(nterms.eq.4) then
         call shape4q(xi,et,pl,xsj2,xs2,shp2,iflag)
      elseif(nterms.eq.6) then
         call shape6tri(xi,et,pl,xsj2,xs2,shp2,iflag)
      else
         call shape3tri(xi,et,pl,xsj2,xs2,shp2,iflag)
      endif
!
!     normal on the surface
!
      dm=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+xsj2(3)*xsj2(3))
      do i=1,3
         xn(i)=xsj2(i)/dm
      enddo
!
!     distance from surface along normal
!
      val=al(1)*xn(1)+al(2)*xn(2)+al(3)*xn(3)
      if(val.le.0.d0) cstr(1)=val
!
!     representative area: usually the slave surface stored in
!     springarea; however, if no area was assigned because the
!     node does not belong to any element, the master surface
!     is used
!
      if(springarea.le.0.d0) then
         if(nterms.eq.3) then
            springarea=dm/2.d0
         else
            springarea=dm*4.d0
         endif
      endif
!
      if(elcon(1,1,imat).gt.0.d0) then
!
!        exponential overclosure
!
         if(dabs(elcon(2,1,imat)).lt.1.d-30) then
            elas(1)=0.d0
            beta=1.d0
         else
!     
            alpha=elcon(2,1,imat)*springarea
            beta=elcon(1,1,imat)
            if(-beta*val.gt.23.d0-dlog(alpha)) then
               beta=(dlog(alpha)-23.d0)/val
            endif
            elas(1)=dexp(-beta*val+dlog(alpha))
         endif
      else
!     
!        linear overclosure
!
c         eps=1.d-6
         pi=4.d0*datan(1.d0)
         eps=-elcon(1,1,imat)*pi/elcon(2,1,imat)
c         elas(1)=-elcon(1,1,imat)-springarea*elcon(2,1,imat)*val*
c     &            (0.5d0+datan(-val/eps)/pi)
         elas(1)=-springarea*elcon(2,1,imat)*val*
     &            (0.5d0+datan(-val/eps)/pi)
      endif
!
!     forces in the nodes of the contact element
!
      do i=1,3
         fnl(i,nope)=-elas(1)*xn(i)
      enddo
      if(iener.eq.1) then
         senergy=elas(1)/beta;
      endif
c      write(*,*) 'springforc ',konl(nope),val,(-fnl(i,nope),i=1,3)
      cstr(4)=elas(1)/springarea
c      write(*,*) 'springforc ',konl(nope),cstr(4)
!
!     Coulomb friction for static calculations
!
      if(ncmat_.ge.7) then
         if(nmethod.eq.1) then
            um=elcon(6,1,imat)
            if(um.gt.0.d0) then
               eps=elcon(7,1,imat)
               pi=4.d0*datan(1.d0)
               do i=1,3
                  al(i)=al(i)-vl(i,nope)
               enddo
               val=al(1)*xn(1)+al(2)*xn(2)+al(3)*xn(3)
!
!              t is the vector connecting the undeformed position
!              of the slave node with the projection on the master
!              face of its deformed position
!
               do i=1,3
                  t(i)=al(i)-val*xn(i)
               enddo
               dt=dsqrt(t(1)*t(1)+t(2)*t(2)+t(3)*t(3))
               if(dt.lt.1.d-20) then
                  c2=1.d0/eps
               else
                  c2=datan(dt/eps)/dt
               endif
               c2=-um*2.d0*c2/pi
               c3=c2*elas(1)
               do i=1,3
                  fnl(i,nope)=fnl(i,nope)+c3*t(i)
               enddo
            endif
            cstr(2)=t(1)
            cstr(3)=t(2)
            cstr(5)=fnl(1,nope)/springarea
            cstr(6)=fnl(2,nope)/springarea
         endif
      endif
!
!     force in the master nodes
!
      do i=1,3
         do j=1,nterms
            fnl(i,j)=-ratio(j)*fnl(i,nope)
         enddo
      enddo
!
      return
      end

