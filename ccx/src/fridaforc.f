!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2013 Guido Dhondt
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
      subroutine fridaforc(xl,konl,vl,imat,elcon,nelcon,
     &  elas,fnl,ncmat_,ntmat_,nope,lakonl,t0l,t1l,kode,elconloc,
     &  plicon,nplicon,npmat_,veoldl,senergy,iener,cstr,mi,springarea)
!
!     calculates the force of the spring
!
      implicit none
!
      character*8 lakonl
!
      integer konl(9),i,j,imat,ncmat_,ntmat_,nope,nterms,iflag,mi(*),
     &  kode,nplicon(0:ntmat_,*),npmat_,nelcon(2,*),iener
!
      real*8 xl(3,9),elas(21),ratio(9),t0l,t1l,al(3),vl(0:mi(2),9),
     &  pl(3,9),xn(3),areamaster,alpha,beta,fnl(3,9),
     &  veoldl(0:mi(2),9),dist,springarea,
     &  elcon(0:ncmat_,ntmat_,*),pproj(3),xsj2(3),xs2(3,7),val,
     &  shp2(7,8),xi,et,elconloc(21),plconloc(802),
     &  plicon(0:2*npmat_,ntmat_,*),fn,
     &  damp,c0,eta,um,eps,fnd(3,9),fnv(3,9),ver(3),dvernor,
     &  dampforc,vertan(3),dvertan,fricforc,pi,senergy,cstr(6)
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
      areamaster=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+xsj2(3)*xsj2(3))
      do i=1,3
         xn(i)=xsj2(i)/areamaster
      enddo
!
!     distance from surface along normal
!
      val=al(1)*xn(1)+al(2)*xn(2)+al(3)*xn(3)
!
!     representative area: usually the slave surface stored in
!     springarea; however, if no area was assigned because the
!     node does not belong to any element, the master surface
!     is used
!
      if(springarea.le.0.d0) then
         if(nterms.eq.3) then
            springarea=areamaster/2.d0
         else
            springarea=areamaster*4.d0
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
         pi=4.d0*datan(1.d0)
         eps=-elcon(1,1,imat)*pi/elcon(2,1,imat)
         elas(1)=-springarea*elcon(2,1,imat)*val*
     &            (0.5d0+datan(-val/eps)/pi)
c     &            -elcon(1,1,imat)*springarea
      endif
!
!     forces in the nodes of the contact element
!
c      do i=1,3
c         do j=1,nterms
c            fnl(i,j)=0.
c         enddo
c         fnl(i,nope)=0.
c      enddo
c      if(iener.eq.1) then
c         senergy=elas(1)/beta;
c      endif
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
            if(val.gt.c0) then
               dampforc=0.d0
            elseif(val.gt.eta*c0) then
               dampforc=dvernor*(c0-val)/(c0*(1.d0-eta))*damp*springarea
            else
               dampforc=dvernor*damp*springarea
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
         um=elcon(6,1,imat)
         if(um.gt.0.d0) then
            if(damp.le.0.d0) then
!     
!     calculate the relative velocity
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
!     calculate the tangential relative velocity
!     
            do i=1,3
               vertan(i)=ver(i)-dvernor*xn(i)
            enddo
            dvertan=dsqrt(vertan(1)**2+vertan(2)**2+vertan(3)**2)
c    write(*,*) 'dvertan ',dvertan
!     
!     normalizing the tangent vector
!
            if(dvertan.gt.0.d0)then      
               do i=1,3
                  vertan(i)=vertan(i)/dvertan
               enddo
            endif
!     
!     friction constants
!     
            eps=elcon(7,1,imat)
!
!     normal force
!
            fn=elas(1)
!
!     modify the friction force in case of contact damping
!
            if(damp.gt.0.d0) fn=fn+dampforc
!     
            fricforc=2.d0*um*datan(dvertan/eps)*fn/pi
!     
            do i=1,3
               do j=1,nterms
c                  fnv(i,j)=ratio(j)*fricforc*vertan(i)
                  fnv(i,j)=-ratio(j)*fricforc*vertan(i)
               enddo
c               fnv(i,nope)=-fricforc*vertan(i)
               fnv(i,nope)=fricforc*vertan(i)
            enddo
         endif
      endif
!
!     summing all forces
!
      if(ncmat_.ge.5) then
         if(damp.gt.0.d0) then
            do j=1,nope
               do i=1,3
                  fnl(i,j)=fnd(i,j)
               enddo
            enddo
         endif
      endif
      if(ncmat_.ge.7) then
         if((um.gt.0.d0).and.(val.lt.0.d0)) then
            if(damp.gt.0.d0) then
               do j=1,nope
                  do i=1,3
                     fnl(i,j)=fnd(i,j)+fnv(i,j)
                  enddo
               enddo
            else
               do j=1,nope
                  do i=1,3
                     fnl(i,j)=fnv(i,j)
                  enddo
               enddo
c       write(*,*) 'fnl(2,nope) ',fnl(2,nope)       
            endif
         else
            do j=1,nope
                do i=1,3
                   fnl(i,j)=0.d0
                enddo
            enddo 
         endif
      endif
!
      return
      end

