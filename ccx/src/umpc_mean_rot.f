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
      subroutine umpc_mean_rot(x,u,f,a,jdof,n,force,iit,idiscon)
!
!     updates the coefficients in a mean rotation mpc
!
!     INPUT:
!
!     x(3,n)             Carthesian coordinates of the nodes in the
!                        user mpc.
!     u(3,n)             Actual displacements of the nodes in the
!                        user mpc.     
!     jdof               Actual degrees of freedom of the mpc terms
!     n                  number of terms in the user mpc
!     force              Actual value of the mpc force
!     iit                iteration number
!
!     OUTPUT:
!
!     f                  Actual value of the mpc. If the mpc is
!                        exactly satisfied, this value is zero
!     a(n)               coefficients of the linearized mpc
!     jdof               Corrected degrees of freedom of the mpc terms
!     idiscon            0: no discontinuity
!                        1: discontinuity
!                        If a discontinuity arises the previous
!                        results are not extrapolated at the start of
!                        a new increment
!
      implicit none
!
      integer jdof(*),n,nkn,i,j,k,imax,iit,idiscon
!
      real*8 x(3,*),u(3,*),f,a(*),aa(3),cgx(3),cgu(3),pi(3),
     &  xi(3),dd,al,a1,amax,c1,c2,c3,c4,c9,c10,force
!
      nkn=(n-1)/3
      if(3*nkn.ne.n-1) then
         write(*,*)
     &     '*ERROR in meanrotmpc: MPC has wrong number of terms'
         stop
      endif
!
!     normal along the rotation axis
!
      dd=0.d0
      do i=1,3
         aa(i)=x(i,n)
         dd=dd+aa(i)**2
      enddo
      dd=dsqrt(dd)
      if(dd.lt.1.d-10) then
         write(*,*) 
     &     '*ERROR in meanrotmpc: rotation vector has zero length'
         stop
      endif
      do i=1,3
         aa(i)=aa(i)/dd
      enddo
!
!     finding the center of gravity of the position and the
!     displacements of the nodes involved in the MPC
!
      do i=1,3
         cgx(i)=0.d0
         cgu(i)=0.d0
      enddo
!
      do i=1,nkn
c         write(*,*) 'x,u'
c         write(*,101) (x(j,3*i-2),j=1,3),(u(j,3*i-2),j=1,3)
c 101     format(6(1x,e11.4))
         do j=1,3
            cgx(j)=cgx(j)+x(j,3*i-2)
            cgu(j)=cgu(j)+u(j,3*i-2)
         enddo
      enddo
!
      do i=1,3
         cgx(i)=cgx(i)/nkn
         cgu(i)=cgu(i)/nkn
      enddo
c      write(*,*) 'cgx ',(cgx(i),i=1,3)
c      write(*,*) 'cgu ',(cgu(i),i=1,3)
!
!     initializing a
!
      do i=1,n
         a(i)=0.d0
      enddo
!
!     calculating the partial derivatives and storing them in a
!
      f=0.d0
      do i=1,nkn
!
!        relative positions
!
         do j=1,3
            pi(j)=x(j,3*i-2)-cgx(j)
            xi(j)=u(j,3*i-2)-cgu(j)+pi(j)
         enddo
!
         c1=pi(1)*pi(1)+pi(2)*pi(2)+pi(3)*pi(3)
         if(c1.lt.1.d-20) then
            write(*,*) '*WARNING in meanrotmpc: node on rotation axis'
            cycle
         endif
         c3=xi(1)*xi(1)+xi(2)*xi(2)+xi(3)*xi(3)
         c2=dsqrt(c1*c3)
!
         al=(aa(1)*pi(2)*xi(3)+aa(2)*pi(3)*xi(1)+aa(3)*pi(1)*xi(2)
     &     -aa(3)*pi(2)*xi(1)-aa(1)*pi(3)*xi(2)-aa(2)*pi(1)*xi(3))
     &     /c2
!
         f=f+dasin(al)
c         write(*,*) 'f ',dasin(al)
!
         do j=1,3
            if(j.eq.1) then
               c4=aa(2)*pi(3)-aa(3)*pi(2)
            elseif(j.eq.2) then
               c4=aa(3)*pi(1)-aa(1)*pi(3)
            else
               c4=aa(1)*pi(2)-aa(2)*pi(1)
            endif
            c9=(c4/c2-al*xi(j)/c3)/dsqrt(1.d0-al*al)
!
            do k=1,nkn
               if(i.eq.k) then
                  c10=c9*(1.d0-1.d0/real(nkn))
               else
                  c10=-c9/real(nkn)
               endif
               a(k*3-3+j)=a(k*3-3+j)+c10
            enddo
         enddo
      enddo
      a(n)=-nkn
      f=f-nkn*u(1,n)
!
!     assigning the degrees of freedom
!
      do i=1,nkn
         jdof(i*3-2)=1
         jdof(i*3-1)=2
         jdof(i*3)=3
      enddo
      jdof(n)=1
!
!     looking for the maximum tangent to decide which DOF should be
!     taken to be the dependent one
!
      if(dabs(a(1)).lt.1.d-5) then
         amax=0.d0
         do i=1,3
            if(dabs(a(i)).gt.amax) then
               amax=abs(a(i))
               imax=i
            endif
         enddo
c         write(*,*) 'a(1),a(2),a(3) ',a(1),a(2),a(3)
c         write(*,*) 'jdof ',jdof(1),jdof(2),jdof(3)
!
         jdof(1)=imax
         a1=a(1)
         a(1)=a(imax)
         do i=2,3
            if(i.eq.imax) then
               jdof(i)=1
               a(i)=a1
               write(*,*) '*INFO: DOF in umpc_mean_rot changed'
c     stop
            else
               jdof(i)=i
            endif
         enddo
      endif
c      write(*,*) 'a(1),a(2),a(3) ',a(1),a(2),a(3)
c      write(*,*) 'jdof ',jdof(1),jdof(2),jdof(3)
!
      return
      end
