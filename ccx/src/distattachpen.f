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
      subroutine distattachpen(xig,etg,pneigh,pnode,a,p,ratio,nterms,
     &                         xnormastface)
!
!     calculates the distance between the node with coordinates
!     in "pnode" and the node with local coordinates xig and etg
!     in a face described by "nterms" nodes with coordinates
!     in pneigh. The distance vector is perpendicular to the local
!     normal "xn(*)" going through the node with the local coordinates xig
!     and etg and is called "xa(*)".
!
      implicit none
!
      integer nterms,i,j
!
      real*8 ratio(8),pneigh(3,*),pnode(3),a,xi,et,xig,etg,p(3),
     &  dummy,xn(3),xnormastface(3,8),xa(3),scal1
!
      if(nterms.eq.3) then
         xi=(xig+1.d0)/2.d0
         et=(etg+1.d0)/2.d0
         if(xi+et.gt.1.d0) then
            dummy=xi
            xi=1.d0-et
            et=1.d0-dummy
         endif
         ratio(1)=1.d0-xi-et
         ratio(2)=xi
         ratio(3)=et
      elseif(nterms.eq.4) then
         xi=xig
         et=etg
         ratio(1)=(1.d0-xi)*(1.d0-et)/4.d0
         ratio(2)=(1.d0+xi)*(1.d0-et)/4.d0
         ratio(3)=(1.d0+xi)*(1.d0+et)/4.d0
         ratio(4)=(1.d0-xi)*(1.d0+et)/4.d0
      elseif(nterms.eq.6) then
         xi=(xig+1.d0)/2.d0
         et=(etg+1.d0)/2.d0
         if(xi+et.gt.1.d0) then
            dummy=xi
            xi=1.d0-et
            et=1.d0-dummy
         endif
         ratio(1)=2.d0*(0.5d0-xi-et)*(1.d0-xi-et)
         ratio(2)=xi*(2.d0*xi-1.d0)
         ratio(3)=et*(2.d0*et-1.d0)
         ratio(4)=4.d0*xi*(1.d0-xi-et)
         ratio(5)=4.d0*xi*et
         ratio(6)=4.d0*et*(1.d0-xi-et)  
      elseif(nterms.eq.8) then
         xi=xig
         et=etg
         ratio(1)=(1.d0-xi)*(1.d0-et)*(-xi-et-1.d0)/4.d0
         ratio(2)=(1.d0+xi)*(1.d0-et)*(xi-et-1.d0)/4.d0
         ratio(3)=(1.d0+xi)*(1.d0+et)*(xi+et-1.d0)/4.d0
         ratio(4)=(1.d0-xi)*(1.d0+et)*(-xi+et-1.d0)/4.d0
         ratio(5)=(1.d0-xi*xi)*(1.d0-et)/2.d0
         ratio(6)=(1.d0+xi)*(1.d0-et*et)/2.d0
         ratio(7)=(1.d0-xi*xi)*(1.d0+et)/2.d0
         ratio(8)=(1.d0-xi)*(1.d0-et*et)/2.d0
      else
         write(*,*) '*ERROR in distattach: case with ',nterms
         write(*,*) '       terms is not covered'
         stop
      endif
!
!     calculating the position in the face and the corresponding normal
!      
      do i=1,3
         p(i)=0.d0
         xn(i)=0.d0
         xa(i)=0.d0
         do j=1,nterms
            p(i)=p(i)+ratio(j)*pneigh(i,j)
            xn(i)=xn(i)+ratio(j)*xnormastface(i,j)
         enddo
      enddo     
!
!     calculating the distance "a" which is perpendicular to the normal
!     through the projection point xig,etg
!
      scal1=((pnode(1)-p(1))*xn(1)+(pnode(2)-p(2))*xn(2)
     &      +(pnode(3)-p(3))*xn(3))/
     &       (xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))      
      do i=1,3
         xa(i)=(pnode(i)-p(i))-scal1*xn(i)
      enddo      
            
      a=xa(1)*xa(1)+xa(2)*xa(2)+xa(3)*xa(3)
!
      return
      end
      
