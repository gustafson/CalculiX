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
      subroutine rectcylvi(co,v,fn,stn,qfn,een,cs,n,icntrl,t,filab,
     &  imag)
!
!     cf. subroutine rectcyl
!     In the present routine, the imaginary part of the 
!     displacements and stresses for all sectors are transformed 
!     from the cylindrical into the rectangular system
!
      implicit none
!
      character*6 filab(*)
      integer i,j,n,icntrl,imag
      real*8 co(3,*),v(0:4,*),fn(0:3,*),stn(6,*),een(6,*),a(3,3),
     &  xr,xt,xz,b(3,3),cs(17,*),t(3),u(3),qfn(3,*),csab(7),
     &  xn(3),r(3),z,theta,rr,c(3,3),ctm,ct,st,ddx,ddy,dd
!
      do i=1,7
         csab(i)=cs(5+i,1)
      enddo
!
      do i=1,n
         j=i
         call transformatrix(csab,co(1,i),a)
!     
         if(filab(11)(1:4).eq.'PU  ') then
            xr=v(1,j)*a(1,1)+v(2,j)*a(1,2)+v(3,j)*a(1,3)
            xt=v(1,j)*a(2,1)+v(2,j)*a(2,2)+v(3,j)*a(2,3)
            xz=v(1,j)*a(3,1)+v(2,j)*a(3,2)+v(3,j)*a(3,3)
            v(1,j)=xr
            v(2,j)=xt
            v(3,j)=xz
         endif
!     
         if(filab(18)(1:4).eq.'PHS ') then
            b(1,1)=stn(1,j)*a(1,1)+stn(4,j)*a(1,2)+stn(5,j)*a(1,3)
            b(1,2)=stn(1,j)*a(2,1)+stn(4,j)*a(2,2)+stn(5,j)*a(2,3)
            b(1,3)=stn(1,j)*a(3,1)+stn(4,j)*a(3,2)+stn(5,j)*a(3,3)
            b(2,1)=stn(4,j)*a(1,1)+stn(2,j)*a(1,2)+stn(6,j)*a(1,3)
            b(2,2)=stn(4,j)*a(2,1)+stn(2,j)*a(2,2)+stn(6,j)*a(2,3)
            b(2,3)=stn(4,j)*a(3,1)+stn(2,j)*a(3,2)+stn(6,j)*a(3,3)
            b(3,1)=stn(5,j)*a(1,1)+stn(6,j)*a(1,2)+stn(3,j)*a(1,3)
            b(3,2)=stn(5,j)*a(2,1)+stn(6,j)*a(2,2)+stn(3,j)*a(2,3)
            b(3,3)=stn(5,j)*a(3,1)+stn(6,j)*a(3,2)+stn(3,j)*a(3,3)
!     
            stn(1,j)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
            stn(2,j)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
            stn(3,j)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
            stn(4,j)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
            stn(5,j)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
            stn(6,j)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
         endif
!     
      enddo
!     
      return
      end
