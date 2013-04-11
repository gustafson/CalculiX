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
     &  plicon,nplicon,npmat_,veoldl)
!
!     calculates the force of the spring
!
      implicit none
!
      character*8 lakonl
!
      integer konl(9),i,j,imat,ncmat_,ntmat_,nope,nterms,iflag,
     &  kode,niso,id,nplicon(0:ntmat_,*),npmat_,nelcon(2,*)
!
      real*8 xl(3,9),elas(21),ratio(9),t0l,t1l,vr(3),vl(0:3,9),
     &  pl(0:3,9),xn(3),al,area,alpha,beta,fnl(3,9),veoldl(3,9),
     &  elcon(0:ncmat_,ntmat_,*),pproj(3),xsj2(3),xs2(3,2),dist,
     &  shp2(4,8),xi,et,elconloc(21),plconloc(82),xk,fk,dd,
     &  xiso(20),yiso(20),dd0,plicon(0:2*npmat_,ntmat_,*),
     &  damp,c0,eta,um,eps,fnd(3,9),fnv(3,9),ver(3),dvernor,
     &  dampforc,vertan(3),dvertan,fricforc,pi
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
         al=dd-dd0
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
            fk=xk*al
         else
            niso=int(plconloc(81))
            do i=1,niso
               xiso(i)=plconloc(2*i-1)
               yiso(i)=plconloc(2*i)
            enddo
            call ident(xiso,al,niso,id)
            if(id.eq.0) then
               xk=0.d0
               fk=yiso(1)
            elseif(id.eq.niso) then
               xk=0.d0
               fk=yiso(niso)
            else
               xk=(yiso(id+1)-yiso(id))/(xiso(id+1)-xiso(id))
               fk=yiso(id)+xk*(al-xiso(id))
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
      call attach(pl,pproj,nterms,ratio,dist,xi,et)
      do i=1,3
         vr(i)=pl(i,nope)-pproj(i)
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
      area=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+xsj2(3)*xsj2(3))
      do i=1,3
         xn(i)=xsj2(i)/area
      enddo
!
!     distance from surface along normal
!
      dist=vr(1)*xn(1)+vr(2)*xn(2)+vr(3)*xn(3)
!
!     representative area
!
      if((nterms.eq.8).or.(nterms.eq.4)) then
         area=area*4.d0
c         area=area*4.d0/konl(nope+1)
      else
         area=area/2.d0
c         area=area/2.d0/konl(nope+1)
      endif
!
      alpha=elcon(2,1,imat)*area
      beta=elcon(1,1,imat)
      if(-beta*dist.gt.23.d0-dlog(alpha)) then
         beta=(dlog(alpha)-23.d0)/dist
      endif
      elas(1)=alpha*dexp(-beta*dist)
      elas(2)=-beta*elas(1)
!
!     forces in the nodes of the contact element
!
      do i=1,3
         do j=1,nterms
            fnl(i,j)=ratio(j)*elas(1)*xn(i)
         enddo
         fnl(i,nope)=-elas(1)*xn(i)
      enddo
c      write(*,*) 'dist,-fn ',konl(nope),dist,(-fnl(i,nope),i=1,3)
!
!     contact damping
!
      if(ncmat_.ge.5) then
         damp=elcon(3,1,imat)
         if(damp.gt.0.d0) then
!
!           calculate the relative velocity
!
            do i=1,3
               ver(i)=0.d0
               do j=1,nterms
                  ver(i)=ver(i)+ratio(j)*veoldl(i,j)
               enddo
               ver(i)=veoldl(i,nope)-ver(i)
            enddo
            dvernor=ver(1)*xn(1)+ver(2)*xn(2)+ver(3)*xn(3)
!
            c0=elcon(4,1,imat)
            eta=elcon(5,1,imat)
!
            if(dist.gt.c0) then
               dampforc=0.d0
            elseif(dist.gt.eta*c0) then
               dampforc=dvernor*(c0-dist)/(c0*(1.d0-eta))*damp*area
            else
               dampforc=dvernor*damp*area
            endif
!
            do i=1,3
               do j=1,nterms
                  fnd(i,j)=ratio(j)*dampforc*xn(i)
               enddo
               fnd(i,nope)=-dampforc*xn(i)
            enddo
         endif
      endif
!
!     friction
!
      if(ncmat_.ge.7) then
         if(damp.le.0.d0) then
!
!           calculate the relative velocity
!
            do i=1,3
               ver(i)=0.d0
               do j=1,nterms
                  ver(i)=ver(i)+ratio(j)*veoldl(i,j)
               enddo
               ver(i)=veoldl(i,nope)-ver(i)
            enddo
            dvernor=ver(1)*xn(1)+ver(2)*xn(2)+ver(3)*xn(3)
         endif
!
         pi=4.d0*datan(1.d0)
!
!        calculate the tangential relative velocity
!
         do i=1,3
            vertan(i)=ver(i)-dvernor*xn(i)
         enddo
         dvertan=dsqrt(vertan(1)**2+vertan(2)**2+vertan(3)**2)
!
!        normalizing the tangent vector
!
         do i=1,3
            vertan(i)=vertan(i)/dvertan
         enddo
!
!        friction constants
!
         um=elcon(6,1,imat)
         eps=elcon(7,1,imat)
!
         fricforc=2.d0*um*datan(dvertan/eps)/pi
!
!        modify the friction force in case of contact damping
!
         if(damp.gt.0.d0) fricforc=fricforc+dampforc*dvertan/dvernor
!
         do i=1,3
            do j=1,nterms
               fnv(i,j)=ratio(j)*fricforc*vertan(i)
            enddo
            fnv(i,nope)=-fricforc*vertan(i)
         enddo
      endif
!
!     summing all forces
!
      if(ncmat_.ge.5) then
         if(damp.gt.0.d0) then
            do j=1,nope
               do i=1,3
                  fnl(i,j)=fnl(i,j)+fnd(i,j)
               enddo
            enddo
         endif
      endif
      if(ncmat_.ge.7) then
         do j=1,nope
            do i=1,3
               fnl(i,j)=fnl(i,j)+fnv(i,j)
            enddo
         enddo
      endif
!
      return
      end

