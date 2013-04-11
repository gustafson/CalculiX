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
     &  elas,fn,ncmat_,ntmat_,nope,lakonl,t0l,t1l,kode,elconloc,
     &  plicon,nplicon,npmat_)
!
!     calculates the force of the spring
!
      implicit none
!
      character*8 lakonl
!
      integer konl(20),i,j,imat,ncmat_,ntmat_,nope,nterms,iflag,
     &  kode,niso,id,nplicon(0:ntmat_,*),npmat_,nelcon(2,*)
!
      real*8 xl(3,20),elas(21),ratio(9),t0l,t1l,dc(3),vl(0:3,20),
     &  pl(0:3,9),xn(3),al,dd,alpha,beta,fn(0:3,*),
     &  elcon(0:ncmat_,ntmat_,*),pproj(3),xsj2(3),xs2(3,2),dist,
     &  shp2(4,8),xi,et,elconloc(21),plconloc(82),xk,fk,
     &  xiso(20),yiso(20),dd0,plicon(0:2*npmat_,ntmat_,*)
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
            fn(i,konl(1))=-fk*xn(i)
            fn(i,konl(2))=fk*xn(i)
         enddo
         return
      endif
!
      nterms=nope-1
!
!     vector dc connects the dependent node with its projection
!     on the independent face
!
      do i=1,3
         pproj(i)=pl(i,nope)
      enddo
      call attach(pl,pproj,nterms,ratio,dist,xi,et)
      do i=1,3
         dc(i)=pl(i,nope)-pproj(i)
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
      dd=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+xsj2(3)*xsj2(3))
      do i=1,3
         xn(i)=xsj2(i)/dd
      enddo
!
!     distance from surface along normal
!
      al=dc(1)*xn(1)+dc(2)*xn(2)+dc(3)*xn(3)
!
!     representative area
!
      if((nterms.eq.8).or.(nterms.eq.4)) then
         dd=dd*4.d0
      else
         dd=dd/2.d0
      endif
!
      alpha=elcon(2,1,imat)*dd
      beta=elcon(1,1,imat)
      if(-beta*al.gt.23.d0-dlog(alpha)) then
         beta=(dlog(alpha)-23.d0)/al
      endif
      elas(1)=alpha*dexp(-beta*al)
      elas(2)=-beta*elas(1)
!
!     forces in the nodes of the contact element
!
      do i=1,3
         do j=1,nterms
            fn(i,konl(j))=fn(i,konl(j))+ratio(j)*elas(1)*xn(i)
         enddo
         fn(i,konl(nope))=fn(i,konl(nope))-elas(1)*xn(i)
      enddo
!
      return
      end

