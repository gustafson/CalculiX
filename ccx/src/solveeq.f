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
      subroutine solveeq(adb,aub,adl,addiv,b,sol,aux,icol,irow,jq,
     &  neq,nzs,nzl)
!
!     solving a system of equations by iteratively solving the
!     lumped version
!     The diagonal terms f the original system are stored in adb,
!     the off-diagonal terms in aub
!     Ref: The Finite Element Method for Fluid Dynamics,
!          O.C. Zienkiewicz, R.L. Taylor & P. Nithiarasu
!          6th edition (2006) ISBN 0 7506 6322 7
!          p. 61
!
      implicit none
!
      integer icol(*),irow(*),jq(*),neq,nzs,nzl,i,j,k
!
      real*8 adb(*),aub(*),adl(*),addiv(*),b(*),sol(*),aux(*),p
!
!     lumping the matrix set (adb,aux) and storing the resulting
!     diagonal terms in adl
!
      do i=1,neq
         adl(i)=adb(i)
      enddo
!
      do j=1,neq
         do k=jq(j),jq(j+1)-1
            i=irow(k)
            adl(i)=adl(i)+aub(k)
            adl(j)=adl(j)+aub(k)
         enddo
      enddo
!
!     storing the difference in addiv
!
      do i=1,neq
         addiv(i)=adb(i)-adl(i)
      enddo
!
!     setting the solution to zero
!
      do i=1,neq
         sol(i)=0.d0
      enddo
!
!     iterating three times
!
c      write(*,*) 'solveeq.... ',adl(1),adl(2),adl(3)
      do k=1,20
c      do k=1,5
!
!        multiplying the difference of the original matrix
!        with the lumped matrix with the actual solution 
!
         call op(neq,p,sol,aux,addiv,aub,icol,irow,nzl)
!
         do i=1,neq
            sol(i)=(b(i)-aux(i))/adl(i)
         enddo
c            write(*,*) 'solveeq... ',k,sol(1),sol(2),sol(3)
!
      enddo
!
      return
      end
