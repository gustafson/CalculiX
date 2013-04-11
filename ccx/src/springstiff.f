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
      subroutine springstiff(xl,elas,konl,voldl,s,imat,elcon,nelcon,
     &  ncmat_,ntmat_,nope,lakonl,t0l,t1l,kode,elconloc,plicon,
     &  nplicon,npmat_,iperturb)
!
!     calculates the stiffness of a spring
!
      implicit none
!
      character*8 lakonl
!
      integer konl(20),i,j,imat,ncmat_,ntmat_,k,l,nope,nterms,iflag,
     &  i1,kode,niso,id,nplicon(0:ntmat_,*),npmat_,nelcon(2,*),
     &  iperturb
!
      real*8 xl(3,9),elas(21),ratio(9),q(3),dist,shp2(7,9),
     &  dc(3),s(60,60),voldl(3,9),pl(3,9),xn(3),al,dd,
     &  c1,c2,c3,c4,alpha,beta,elcon(0:ncmat_,ntmat_,*),xsj2(3),
     &  xsju(3,3,9),dxsju(3,9),h(3,9),fpu(3,3,9),xi,et,
     &  xs2(3,7),t0l,t1l,elconloc(21),plconloc(82),xk,fk,
     &  xiso(20),yiso(20),dd0,plicon(0:2*npmat_,ntmat_,*),
     &  a11,a12,a22,b1(3,9),b2(3,9),dal(3,3,9),qxxy(3),fnl(3),
     &  qxyy(3),dxi(3,9),det(3,9),determinant,c11,c12,c22,
     &  qxyx(3),qyxy(3)
!
      data iflag /4/
!
!     actual positions of the nodes belonging to the contact spring
!
      if(iperturb.eq.0) then
         do i=1,nope
            do j=1,3
               pl(j,i)=xl(j,i)
            enddo
         enddo
      else
         do i=1,nope
            do j=1,3
               pl(j,i)=xl(j,i)+voldl(j,i)
            enddo
         enddo
      endif
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
         c1=fk/dd
         c2=xk-c1
         do i=1,3
            do j=1,3
               s(i,j)=c2*xn(i)*xn(j)
            enddo
            s(i,i)=s(i,i)+c1
         enddo
         do i=1,3
            do j=1,3
               s(i+3,j)=-s(i,j)
               s(i,j+3)=-s(i,j)
               s(i+3,j+3)=s(i,j)
            enddo
         enddo
         return
      endif
!
!     contact springs
!
      nterms=nope-1
!
!     vector dc connects the dependent node with its projection
!     on the independent face
!
      do i=1,3
         q(i)=pl(i,nope)
      enddo
      call attach(pl,q,nterms,ratio,dist,xi,et)
      do i=1,3
         dc(i)=pl(i,nope)-q(i)
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
!     dxi(i,j) is the derivative of xi w.r.t. pl(i,j),
!     det(i,j) is the derivative of eta w.r.t. pl(i,j)
!
!     dxi and det are determined from the orthogonality
!     condition
!
      a11=-(xs2(1,1)*xs2(1,1)+xs2(2,1)*xs2(2,1)+xs2(3,1)*xs2(3,1))
     &    +dc(1)*xs2(1,5)+dc(2)*xs2(2,5)+dc(3)*xs2(3,5)
      a12=-(xs2(1,1)*xs2(1,2)+xs2(2,1)*xs2(2,2)+xs2(3,1)*xs2(3,2))
     &    +dc(1)*xs2(1,6)+dc(2)*xs2(2,6)+dc(3)*xs2(3,6)
      a22=-(xs2(1,2)*xs2(1,2)+xs2(2,2)*xs2(2,2)+xs2(3,2)*xs2(3,2))
     &    +dc(1)*xs2(1,7)+dc(2)*xs2(2,7)+dc(3)*xs2(3,7)
!
      do i=1,3
         do j=1,nterms
            b1(i,j)=shp2(4,j)*xs2(i,1)-shp2(1,j)*dc(i)
            b2(i,j)=shp2(4,j)*xs2(i,2)-shp2(2,j)*dc(i)
         enddo
         b1(i,nope)=-xs2(i,1)
         b2(i,nope)=-xs2(i,2)
      enddo
!
      determinant=a11*a22-a12*a12
      c11=a22/determinant
      c12=-a12/determinant
      c22=a11/determinant
!
      do i=1,3
         do j=1,nope
            dxi(i,j)=c11*b1(i,j)+c12*b2(i,j)
            det(i,j)=c12*b1(i,j)+c22*b2(i,j)
         enddo
      enddo
!
!     dal(i,j,k) is the derivative of dc(i) w.r.t pl(j,k)
!
      do i=1,nope
         do j=1,3
            do k=1,3
               dal(j,k,i)=-xs2(j,1)*dxi(k,i)-xs2(j,2)*det(k,i)
            enddo
         enddo
      enddo
      do i=1,nterms
         do j=1,3
            dal(j,j,i)=dal(j,j,i)-shp2(4,i)
         enddo
      enddo
      do j=1,3
         dal(j,j,nope)=dal(j,j,nope)+1.d0
      enddo
!
!     d2q/dxx x dq/dy
!
      qxxy(1)=xs2(2,5)*xs2(3,2)-xs2(3,5)*xs2(2,2)
      qxxy(2)=xs2(3,5)*xs2(1,2)-xs2(1,5)*xs2(3,2)
      qxxy(3)=xs2(1,5)*xs2(2,2)-xs2(2,5)*xs2(1,2)
!
!     dq/dx x d2q/dyy
!
      qxyy(1)=xs2(2,1)*xs2(3,7)-xs2(3,1)*xs2(2,7)
      qxyy(2)=xs2(3,1)*xs2(1,7)-xs2(1,1)*xs2(3,7)
      qxyy(3)=xs2(1,1)*xs2(2,7)-xs2(2,1)*xs2(1,7)
!
!     Modified by Stefan Sicklinger
!
!     dq/dx x d2q/dxy
!
      qxyx(1)=xs2(2,1)*xs2(3,6)-xs2(3,1)*xs2(2,6)
      qxyx(2)=xs2(3,1)*xs2(1,6)-xs2(1,1)*xs2(3,6)
      qxyx(3)=xs2(1,1)*xs2(2,6)-xs2(2,1)*xs2(1,6)
!
!
!     d2q/dxy x dq/dy
!
      qyxy(1)=xs2(2,6)*xs2(3,2)-xs2(3,6)*xs2(2,2)
      qyxy(2)=xs2(3,6)*xs2(1,2)-xs2(1,6)*xs2(3,2)
      qyxy(3)=xs2(1,6)*xs2(2,2)-xs2(2,6)*xs2(1,2)
!
!
!     End modifications
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
c      write(*,*) 'springstiff ',al
!
!     alpha and beta, taking the representative area into account
!     (conversion of pressure into force)
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
               alpha=elcon(2,1,imat)*dd*4.d0
c     alpha=elcon(2,1,imat)*dd*4.d0/konl(nope+1)
            else
               alpha=elcon(2,1,imat)*dd/2.d0
c     alpha=elcon(2,1,imat)*dd/2.d0/konl(nope+1)
            endif
            beta=elcon(1,1,imat)
            if(-beta*al.gt.23.d0-dlog(alpha)) then
               beta=(dlog(alpha)-23.d0)/al
            endif
            elas(1)=dexp(-beta*al+dlog(alpha))
            elas(2)=-beta*elas(1)
         endif
      else
!     
!        linear overclosure
!
         elas(1)=-dd*elcon(2,1,imat)*dist
         elas(2)=-dd*elcon(2,1,imat)
      endif
!
!     contact force
!
      do i=1,3
         fnl(i)=-elas(1)*xn(i)
      enddo
!
!     derivatives of the jacobian vector w.r.t. the displacement
!     vectors
!
      do k=1,nterms
         xsju(1,1,k)=0.d0
         xsju(2,2,k)=0.d0
         xsju(3,3,k)=0.d0
         xsju(1,2,k)=shp2(1,k)*xs2(3,2)-shp2(2,k)*xs2(3,1)
         xsju(2,3,k)=shp2(1,k)*xs2(1,2)-shp2(2,k)*xs2(1,1)
         xsju(3,1,k)=shp2(1,k)*xs2(2,2)-shp2(2,k)*xs2(2,1)
         xsju(1,3,k)=-xsju(3,1,k)
         xsju(2,1,k)=-xsju(1,2,k)
         xsju(3,2,k)=-xsju(2,3,k)
      enddo
      do i=1,3
         do j=1,3
            xsju(i,j,nope)=0.d0
         enddo
      enddo
!
!     correction due to change of xi and eta
!
      do k=1,nope
         do i=1,3
            do j=1,3
!
!     modified by Stefan Sicklinger
!
               xsju(i,j,k)=xsju(i,j,k)+(qxxy(i)+qxyx(i))*dxi(j,k)
     &                                +(qxyy(i)+qyxy(i))*det(j,k)
c               xsju(i,j,k)=xsju(i,j,k)+qxxy(i)*dxi(j,k)
c     &                                +qxyy(i)*det(j,k)
            enddo
         enddo
      enddo
!     
!     derivatives of the size of the jacobian vector w.r.t. the
!     displacement vectors
!
      do k=1,nope
         do i=1,3
            dxsju(i,k)=xn(1)*xsju(1,i,k)+xn(2)*xsju(2,i,k)+
     &           xn(3)*xsju(3,i,k)
         enddo
!
!        auxiliary variables
!
         do i=1,3
            h(i,k)=dc(1)*xsju(1,i,k)+dc(2)*xsju(2,i,k)+
     &             dc(3)*xsju(3,i,k)-al*dxsju(i,k)
         enddo
!
      enddo
!
      c1=1.d0/dd
      c2=c1*c1
      c3=elas(2)*c2
      c4=elas(1)*c1
!
!     derivatives of the forces w.r.t. the displacement vectors
!
      do k=1,nope
         do i=1,3
            do j=1,3
               fpu(i,j,k)=-c3*xsj2(i)*(h(j,k)
     &      +(xsj2(1)*dal(1,j,k)+xsj2(2)*dal(2,j,k)+xsj2(3)*dal(3,j,k)))
     &                    +c4*(xn(i)*dxsju(j,k)-xsju(i,j,k))
            enddo
         enddo
      enddo
!
!     determining the stiffness matrix contributions
!
c      do k=1,nterms
c         ratio(k)=-ratio(k)
c      enddo
c      ratio(nope)=1.d0
!
!     complete field shp2 
!
      shp2(1,nope)=0.d0
      shp2(2,nope)=0.d0
      shp2(4,nope)=-1.d0
!
      do k=1,nope
         do l=1,nope
            do i=1,3
               i1=i+(k-1)*3
               do j=1,3
c                  s(i1,j+(l-1)*3)=ratio(k)*fpu(i,j,l)
                  s(i1,j+(l-1)*3)=-shp2(4,k)*fpu(i,j,l)
     &                            -shp2(1,k)*fnl(i)*dxi(j,l)
     &                            -shp2(2,k)*fnl(i)*det(j,l)
               enddo
            enddo
         enddo
      enddo
!
!     symmetrizing the matrix
!
      do j=1,3*nope
         do i=1,j-1
            s(i,j)=(s(i,j)+s(j,i))/2.d0
c            write(*,*) 'springstiff stiffness ',i,j,s(i,j)
         enddo
      enddo
!
      return
      end

