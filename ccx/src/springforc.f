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
     &  plicon,nplicon,npmat_,veoldl,senergy,iener,cstr,mi,ifricdamp)
!
!     calculates the force of the spring
!
      implicit none
!
      character*8 lakonl
!
      integer konl(9),i,j,imat,ncmat_,ntmat_,nope,nterms,iflag,mi(2),
     &  kode,niso,id,nplicon(0:ntmat_,*),npmat_,nelcon(2,*),iener,
     &  ifricdamp
!
      real*8 xl(3,9),elas(21),ratio(9),t0l,t1l,vr(3),vl(0:mi(2),9),
     &  pl(3,9),xn(3),al,area,alpha,beta,fnl(3,9),veoldl(0:mi(2),9),
     &  elcon(0:ncmat_,ntmat_,*),pproj(3),xsj2(3),xs2(3,7),dist,
     &  shp2(7,8),xi,et,elconloc(21),plconloc(82),xk,fk,dd,
     &  xiso(20),yiso(20),dd0,plicon(0:2*npmat_,ntmat_,*),fn,
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
            if(iener.eq.1) then
               senergy=fk*al/2.d0
            endif
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
               if(iener.eq.1) then
                  senergy=fk*al;
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
                  senergy=senergy+(al-xiso(niso))*yiso(niso)
               endif
            else
               xk=(yiso(id+1)-yiso(id))/(xiso(id+1)-xiso(id))
               fk=yiso(id)+xk*(al-xiso(id))
               if(iener.eq.1) then
                  senergy=yiso(1)*xiso(1)
                  do i=2, id
                     senergy=senergy+(xiso(i)-xiso(i-1))*
     &                    (yiso(i)+yiso(i-1))/2.d0
                  enddo
                  senergy=senergy+(al-xiso(id))*(fk+yiso(id))/2.d0
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
      if(dist.le.0.d0) cstr(1)=-dist
!
!     representative area
!
      if(elcon(1,1,imat).gt.0.d0) then
!
!        exponential overclosure
!
         if(dabs(elcon(2,1,imat)).lt.1.d-30) then
            elas(1)=0.d0
            elas(2)=0.d0
         else
            if((nterms.eq.8).or.(nterms.eq.4)) then
               area=area*4.d0
c     area=area*4.d0/konl(nope+1)
            else
               area=area/2.d0
c     area=area/2.d0/konl(nope+1)
            endif
!     
            alpha=elcon(2,1,imat)*area
            beta=elcon(1,1,imat)
            if(-beta*dist.gt.23.d0-dlog(alpha)) then
               beta=(dlog(alpha)-23.d0)/dist
            endif
            elas(1)=dexp(-beta*dist+dlog(alpha))
            elas(2)=-beta*elas(1)
         endif
      else
!
!        linear overclosure
!
         elas(1)=-area*elcon(2,1,imat)*dist
         elas(2)=-area*elcon(2,1,imat)
      endif
!
!     forces in the nodes of the contact element
!
      do i=1,3
         do j=1,nterms
            fnl(i,j)=ratio(j)*elas(1)*xn(i)
         enddo
         fnl(i,nope)=-elas(1)*xn(i)
      enddo
      if(iener.eq.1) then
         senergy=elas(1)/beta;
      endif
c      write(*,*) 'springforc ',konl(nope),dist,(-fnl(i,nope),i=1,3)
      cstr(4)=elas(1)/area
c      write(*,*) 'springforc ',konl(nope),cstr(4)
!
!     contact damping
!
      if(ncmat_.ge.5) then
         damp=elcon(3,1,imat)
         if(damp.gt.0.d0) then
            ifricdamp=1
         endif
      endif
!
!     friction
!
      if(ncmat_.ge.7) then
         um=elcon(6,1,imat)
         if(um.gt.0.d0) then
            ifricdamp=1
         endif
      endif
!
      return
      end

