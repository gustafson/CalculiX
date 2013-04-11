!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine radresult(ntr,xloadact,ntm,bc,nloadtr,tarea,
     &     tenv,physcon,erad,f,fenv)
!     
      implicit none
!     
      integer i,j,ntm,ntr,nloadtr(*)
!     
      real*8 xloadact(2,*), tarea(*),tenv(*),
     &     erad(*),q,fenv(*),physcon(3),f(ntr,*),bc(ntm,1)
!     
!     calculating the flux and transforming the flux into an
!     equivalent temperature
!     
      write(*,*) ''
      
      do i=1,ntr
         q=bc(i,1)
         do j=1,ntr
            if(i.eq.j)cycle
            q=q-f(i,j)*bc(j,1)
         enddo
         q=q-fenv(i)*physcon(2)*tenv(i)**4
         xloadact(2,nloadtr(i))=
     &        max(tarea(i)**4-q/(erad(i)*physcon(2)),0.d0)
c     write(*,*)  xloadact(2,nloadtr(i))
         xloadact(2,nloadtr(i))=
     &        (xloadact(2,nloadtr(i)))**0.25+physcon(1)
c     write(*,*)  xloadact(2,nloadtr(i)) 
c     write(*,*) i,bc(i,1),q,xloadact(2,nloadtr(i)
      enddo
!     
      return
      end
